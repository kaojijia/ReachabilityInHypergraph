#include "Hypergraph.h"
#include "HypergraphTreeIndex.h"
#include "UWeightedPLL.h"
#include <gtest/gtest.h>
#include <fstream>
#include <iostream>
#include <chrono>
#include <random>
#include <vector>
#include <iomanip>
#include <thread>
#include <atomic>
#include <map>
#include <random>
#include <set>
#include <filesystem>
#include <queue>
#include <unordered_set>
#include <optional>

using namespace std;
using namespace std::chrono;

namespace BenchmarkConfig
{
    std::string g_dataset_file = "";
    int g_num_add_vertex_ops = 1000;
    int g_num_remove_vertex_ops = 1000;
}

static std::string getCurrentTimestamp()
{
    auto now = std::chrono::system_clock::now();
    auto now_time_t = std::chrono::system_clock::to_time_t(now);
    auto now_us = std::chrono::duration_cast<std::chrono::microseconds>(now.time_since_epoch()) % 1000000;
    tm local_tm;
#if defined(_WIN32) || defined(_WIN64)
    localtime_s(&local_tm, &now_time_t);
#else
    localtime_r(&now_time_t, &local_tm);
#endif

    std::stringstream ss;
    ss << std::put_time(&local_tm, "%Y-%m-%d %H:%M:%S");
    ss << "." << std::setw(6) << std::setfill('0') << now_us.count();
    return "[" + ss.str() + "] ";
}

class HypergraphTest : public ::testing::Test
{
protected:
    void SetUp() override {}
    Hypergraph hg;
};

bool bfs_path_exists(Hypergraph &hg, int start_node, int end_node)
{
    if (start_node == end_node)
        return true;
    if (start_node >= hg.numVertices() || end_node >= hg.numVertices() || start_node < 0 || end_node < 0)
        return false;
    if (hg.numVertices() == 0)
        return false;

    std::queue<int> q;
    q.push(start_node);
    std::vector<bool> visited_vertices(hg.numVertices(), false);
    visited_vertices[start_node] = true;

    while (!q.empty())
    {
        int u = q.front();
        q.pop();

        const auto &incident_edges_u = hg.getIncidentHyperedges(u);
        for (int edge_id : incident_edges_u)
        {
            if (edge_id < 0 || edge_id >= hg.numHyperedges())
                continue;
            const auto &hyperedge_obj = hg.getHyperedge(edge_id);

            for (int v_neighbor : hyperedge_obj.vertices)
            {
                if (v_neighbor < 0 || v_neighbor >= hg.numVertices())
                    continue;
                if (v_neighbor == end_node)
                    return true;
                if (!visited_vertices[v_neighbor])
                {
                    visited_vertices[v_neighbor] = true;
                    q.push(v_neighbor);
                }
            }
        }
    }
    return false;
}

TEST_F(HypergraphTest, HyperIndexBenchmark)
{
    std::mt19937 random_generator(32);
    std::string dataset_file_default = "example_hypergraph.txt";
    std::string dataset_file = BenchmarkConfig::g_dataset_file.empty() ? dataset_file_default : BenchmarkConfig::g_dataset_file;
    int num_add_vertex_ops = BenchmarkConfig::g_num_add_vertex_ops;
    int num_remove_vertex_ops = BenchmarkConfig::g_num_remove_vertex_ops;

    std::string cache_file_prefix;
    if (!dataset_file.empty())
    {
        std::filesystem::path dataset_path(dataset_file);
        std::filesystem::path dataset_dir = dataset_path.parent_path();
        std::string dataset_name = dataset_path.stem().string();
        std::filesystem::path cache_dir = dataset_dir / "IndexCache";

        try
        {
            if (!std::filesystem::exists(cache_dir))
            {
                std::filesystem::create_directories(cache_dir);
                std::cout << "创建缓存目录: " << cache_dir << std::endl;
            }
        }
        catch (const std::filesystem::filesystem_error &e)
        {
            std::cerr << "警告: 无法创建缓存目录 " << cache_dir << ": " << e.what() << std::endl;
            cache_dir = std::filesystem::temp_directory_path() / "HypergraphIndexCache";
            std::filesystem::create_directories(cache_dir);
        }

        cache_file_prefix = (cache_dir / (dataset_name + "_")).string();
    }
    else
    {
        std::filesystem::path default_cache_dir = std::filesystem::temp_directory_path() / "HypergraphIndexCache";
        std::filesystem::create_directories(default_cache_dir);
        cache_file_prefix = (default_cache_dir / "default_").string();
    }

    try
    {
        hg = Hypergraph::fromFile(dataset_file);
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error loading hypergraph from file: " << dataset_file << " - " << e.what() << std::endl;
        FAIL() << "无法加载超图: " << dataset_file << " - " << e.what();
        return;
    }

    if (hg.numVertices() == 0 || hg.numHyperedges() == 0)
    {
        GTEST_SKIP() << "超图为空，无法运行基准测试。";
    }

    HypergraphTreeIndex tree_index(hg);

    std::cout << getCurrentTimestamp() << "--- HyperIndex 基准测试 ---" << std::endl;
    std::cout << getCurrentTimestamp() << "数据集: " << dataset_file << std::endl;
    std::cout << getCurrentTimestamp() << "顶点添加操作次数: " << num_add_vertex_ops << std::endl;
    std::cout << getCurrentTimestamp() << "顶点移除操作次数: " << num_remove_vertex_ops << std::endl;

    std::cout << getCurrentTimestamp() << "\n[1. 索引构建]" << std::endl;
    auto build_start_time = high_resolution_clock::now();
    tree_index.buildIndexCacheSizeOnly(cache_file_prefix);
    auto build_end_time = high_resolution_clock::now();
    auto build_duration_ms = duration_cast<milliseconds>(build_end_time - build_start_time);
    std::cout << getCurrentTimestamp() << "索引构建时间: " << build_duration_ms.count() << " ms" << std::endl;

    std::cout << getCurrentTimestamp() << "\n[2. 索引大小]" << std::endl;
    double index_size_mb = tree_index.getMemoryUsageMB();
    std::cout << getCurrentTimestamp() << "索引大小: " << std::fixed << std::setprecision(2) << index_size_mb << " MB" << std::endl;

    std::cout << getCurrentTimestamp() << "\n[3. 动态更新性能]" << std::endl;
    long long total_add_op_time_nanoseconds = 0;
    int add_ops_performed_count = 0;
    long long total_remove_op_time_nanoseconds = 0;
    int remove_ops_performed_count = 0;

    const int MAX_PARAM_GENERATION_ATTEMPTS_PER_OP = 20;
    const int ATTEMPT_MULTIPLIER = 5;

    if (num_add_vertex_ops > 0 && hg.numVertices() > 0 && hg.numHyperedges() > 0)
    {
        int total_add_attempts_overall = 0;
        while (add_ops_performed_count < num_add_vertex_ops && total_add_attempts_overall < num_add_vertex_ops * ATTEMPT_MULTIPLIER)
        {
            total_add_attempts_overall++;
            if (hg.numVertices() == 0 || hg.numHyperedges() == 0)
            {
                break;
            }

            int vertex_id_for_update = -1;
            std::vector<int> edge_list_for_update;
            bool params_generated = false;

            for (int attempt = 0; attempt < MAX_PARAM_GENERATION_ATTEMPTS_PER_OP; ++attempt)
            {
                edge_list_for_update.clear();
                std::uniform_int_distribution<> current_vertex_dist(0, hg.numVertices() - 1);
                vertex_id_for_update = current_vertex_dist(random_generator);

                std::uniform_int_distribution<> edge_list_size_dist(1, std::min(5, (int)hg.numHyperedges()));
                int list_size = edge_list_size_dist(random_generator);
                std::uniform_int_distribution<> current_edge_dist(0, hg.numHyperedges() - 1);
                std::unordered_set<int> chosen_edges_for_this_update;

                for (int j = 0; j < list_size && chosen_edges_for_this_update.size() < hg.numHyperedges(); ++j)
                {
                    int edge_to_select = current_edge_dist(random_generator);
                    if (chosen_edges_for_this_update.find(edge_to_select) == chosen_edges_for_this_update.end())
                    {
                        edge_list_for_update.push_back(edge_to_select);
                        chosen_edges_for_this_update.insert(edge_to_select);
                    }
                }
                if (!edge_list_for_update.empty() && vertex_id_for_update != -1)
                {
                    params_generated = true;
                    break;
                }
            }

            if (!params_generated)
            {
                continue;
            }

            try
            {
                auto update_start_time = high_resolution_clock::now();
                tree_index.applyAddVertexToEdges(vertex_id_for_update, edge_list_for_update);
                auto update_end_time = high_resolution_clock::now();
                total_add_op_time_nanoseconds += duration_cast<nanoseconds>(update_end_time - update_start_time).count();
                add_ops_performed_count++;
            }
            catch (const std::exception &e)
            {
                continue;
            }
        }
    }

    if (num_remove_vertex_ops > 0 && hg.numVertices() > 0 && hg.numHyperedges() > 0)
    {
        int total_remove_attempts_overall = 0;
        while (remove_ops_performed_count < num_remove_vertex_ops && total_remove_attempts_overall < num_remove_vertex_ops * ATTEMPT_MULTIPLIER)
        {
            total_remove_attempts_overall++;
            if (hg.numVertices() == 0 || hg.numHyperedges() == 0)
            {
                break;
            }

            int vertex_id_for_update = -1;
            std::vector<int> edge_list_for_update;
            bool params_generated = false;

            for (int attempt = 0; attempt < MAX_PARAM_GENERATION_ATTEMPTS_PER_OP; ++attempt)
            {
                edge_list_for_update.clear();
                std::uniform_int_distribution<> current_vertex_dist(0, hg.numVertices() - 1);
                int potential_vertex_id = current_vertex_dist(random_generator);

                const auto &incident_edges = hg.getIncidentHyperedges(potential_vertex_id);
                if (incident_edges.empty())
                {
                    continue;
                }
                vertex_id_for_update = potential_vertex_id;

                std::uniform_int_distribution<> edge_list_size_dist(1, std::min(5, (int)incident_edges.size()));
                int list_size = edge_list_size_dist(random_generator);

                std::vector<int> incident_edges_vec(incident_edges.begin(), incident_edges.end());
                std::shuffle(incident_edges_vec.begin(), incident_edges_vec.end(), random_generator);
                for (int j = 0; j < list_size && j < incident_edges_vec.size(); ++j)
                {
                    edge_list_for_update.push_back(incident_edges_vec[j]);
                }

                if (!edge_list_for_update.empty() && vertex_id_for_update != -1)
                {
                    params_generated = true;
                    break;
                }
            }

            if (!params_generated)
            {
                continue;
            }

            try
            {
                auto update_start_time = high_resolution_clock::now();
                tree_index.applyRemoveVertexFromEdges(vertex_id_for_update, edge_list_for_update);
                auto update_end_time = high_resolution_clock::now();
                total_remove_op_time_nanoseconds += duration_cast<nanoseconds>(update_end_time - update_start_time).count();
                remove_ops_performed_count++;
            }
            catch (const std::exception &e)
            {
                continue;
            }
        }
    }

    long long total_update_time_nanoseconds = total_add_op_time_nanoseconds + total_remove_op_time_nanoseconds;
    int total_update_ops_performed = add_ops_performed_count + remove_ops_performed_count;

    if (total_update_ops_performed > 0)
    {
        double avg_update_time_ns = static_cast<double>(total_update_time_nanoseconds) / total_update_ops_performed;
        std::cout << "平均更新操作时间 (" << total_update_ops_performed << " 次总更新, "
                  << add_ops_performed_count << " 添加, " << remove_ops_performed_count << " 移除): "
                  << std::fixed << std::setprecision(2)
                  << avg_update_time_ns << " ns | "
                  << avg_update_time_ns / 1000.0 << " µs | "
                  << avg_update_time_ns / 1000000.0 << " ms" << std::endl;
    }

    std::cout << getCurrentTimestamp() << "\n[4. 操作组合分析]" << std::endl;

    std::optional<double> avg_add_time_ns_opt;
    if (add_ops_performed_count > 0)
    {
        avg_add_time_ns_opt = static_cast<double>(total_add_op_time_nanoseconds) / add_ops_performed_count;
    }

    std::optional<double> avg_remove_time_ns_opt;
    if (remove_ops_performed_count > 0)
    {
        avg_remove_time_ns_opt = static_cast<double>(total_remove_op_time_nanoseconds) / remove_ops_performed_count;
    }

    auto print_projected_time = [](const std::string &label, std::optional<double> time_ns)
    {
        std::cout << std::left << std::setw(40) << label << ": ";
        if (time_ns)
        {
            std::cout << std::fixed << std::setprecision(2)
                      << *time_ns << " ns | "
                      << *time_ns / 1000.0 << " µs | "
                      << *time_ns / 1000000.0 << " ms" << std::endl;
        }
        else
        {
            std::cout << "N/A" << std::endl;
        }
    };

    print_projected_time("  平均时间 (100% 添加)", avg_add_time_ns_opt);
    print_projected_time("  平均时间 (100% 删除)", avg_remove_time_ns_opt);

    std::optional<double> mix_25A_75R_time_ns;
    if (avg_add_time_ns_opt && avg_remove_time_ns_opt)
    {
        mix_25A_75R_time_ns = (0.25 * *avg_add_time_ns_opt) + (0.75 * *avg_remove_time_ns_opt);
    }
    print_projected_time("  平均时间 (25% 添加, 75% 删除)", mix_25A_75R_time_ns);

    std::optional<double> mix_50A_50R_time_ns;
    if (avg_add_time_ns_opt && avg_remove_time_ns_opt)
    {
        mix_50A_50R_time_ns = (0.50 * *avg_add_time_ns_opt) + (0.50 * *avg_remove_time_ns_opt);
    }
    print_projected_time("  平均时间 (50% 添加, 50% 删除)", mix_50A_50R_time_ns);

    std::optional<double> mix_75A_25R_time_ns;
    if (avg_add_time_ns_opt && avg_remove_time_ns_opt)
    {
        mix_75A_25R_time_ns = (0.75 * *avg_add_time_ns_opt) + (0.25 * *avg_remove_time_ns_opt);
    }
    print_projected_time("  平均时间 (75% 添加, 25% 删除)", mix_75A_25R_time_ns);

    std::string actual_cache_file_path = cache_file_prefix + "hypergraph_tree_index_SizeOnly";
    if (std::filesystem::exists(actual_cache_file_path))
    {
        std::filesystem::remove(actual_cache_file_path);
    }

    std::cout << getCurrentTimestamp() << "\n--- HyperIndex 基准测试结束 ---" << std::endl;
}

int main(int argc, char **argv)
{
    for (int i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];
        if (arg == "--dataset" && i + 1 < argc)
        {
            BenchmarkConfig::g_dataset_file = argv[++i];
        }
        else if (arg == "--help")
        {
            std::cout << "Benchmark Options:\n"
                      << "  --dataset <filepath>      Path to the dataset file.\n"
                      << "  --help                    Show this help message.\n";
            return 0;
        }
    }

    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
