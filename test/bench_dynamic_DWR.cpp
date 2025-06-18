#include "testBottleneck.h"
#include "Hypergraph.h"
#include <gtest/gtest.h>
#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <random>
#include <iomanip>
#include <numeric>
#include <algorithm>
#include <optional>
#include <set>
#include <filesystem>

using namespace std;
using namespace std::chrono;

namespace BenchmarkConfigDLCR
{
    std::string g_dataset_file = "";
    int g_num_add_edge_ops = 1000;
    int g_num_remove_edge_ops = 1000;
    int g_k_threshold_conversion = 1;
}

class BottleneckGraphBenchmarkTest : public ::testing::Test
{
protected:
    TestBottleneckGraph *test_graph_ptr = nullptr;
    bottleneckPath::BottleneckGraph *g1_ptr = nullptr;
    Hypergraph hg_source;

    void SetUp() override {}

    void TearDown() override
    {
        delete test_graph_ptr;
        test_graph_ptr = nullptr;
        g1_ptr = nullptr;
    }

    std::vector<std::tuple<int, int, int>> getCurrentEdges()
    {
        std::vector<std::tuple<int, int, int>> edges;
        if (!g1_ptr)
            return edges;
        std::set<std::pair<int, int>> added_undirected_edges;

        for (int u = 1; u <= g1_ptr->n; ++u)
        {
            if (u < g1_ptr->Adj.size())
            {
                for (const auto &edge_pair : g1_ptr->Adj[u])
                {
                    int v = edge_pair.first;
                    int w = edge_pair.second;
                    int u_norm = std::min(u, v);
                    int v_norm = std::max(u, v);
                    if (added_undirected_edges.find({u_norm, v_norm}) == added_undirected_edges.end())
                    {
                        edges.emplace_back(u, v, w);
                        added_undirected_edges.insert({u_norm, v_norm});
                    }
                }
            }
        }
        return edges;
    }
};

TEST_F(BottleneckGraphBenchmarkTest, BottleneckIndexBenchmark)
{

    std::mt19937 random_generator(32);
    std::string dataset_file_default = "example_hypergraph.txt";

    std::string original_hypergraph_file = BenchmarkConfigDLCR::g_dataset_file.empty() ? dataset_file_default : BenchmarkConfigDLCR::g_dataset_file;
    int num_add_edge_ops = BenchmarkConfigDLCR::g_num_add_edge_ops;
    int num_remove_edge_ops = BenchmarkConfigDLCR::g_num_remove_edge_ops;
    int k_threshold_for_conversion = BenchmarkConfigDLCR::g_k_threshold_conversion;

    std::string temp_weighted_graph_file = BenchmarkConfigDLCR::g_dataset_file + "_weighted_graph_k" + std::to_string(k_threshold_for_conversion);

    std::cout << "Loading original hypergraph: " << original_hypergraph_file << std::endl;
    try
    {
        hg_source = Hypergraph::fromFile(original_hypergraph_file);
    }
    catch (const std::exception &e)
    {
        FAIL() << "Cannot load original hypergraph file: " << original_hypergraph_file << " - " << e.what();
        return;
    }
    if (hg_source.numVertices() == 0 || hg_source.numHyperedges() == 0)
    {
        GTEST_SKIP() << "Original hypergraph is empty, cannot run benchmark. File: " << original_hypergraph_file;
        return;
    }
    std::cout << "Original hypergraph loaded: " << hg_source.numVertices() << " vertices, " << hg_source.numHyperedges() << " hyperedges." << std::endl;

    bool weighted_graph_loaded_from_cache = false;
    if (std::filesystem::exists(temp_weighted_graph_file))
    {
        std::cout << "Found existing weighted graph file: " << temp_weighted_graph_file << std::endl;
        std::cout << "Attempting to load directly..." << std::endl;
        try
        {
            test_graph_ptr = new TestBottleneckGraph(temp_weighted_graph_file, true, false);
            if (test_graph_ptr && test_graph_ptr->g1 && test_graph_ptr->g1->n > 0)
            {
                g1_ptr = test_graph_ptr->g1;
                std::cout << "Successfully loaded weighted graph from cache." << std::endl;
                weighted_graph_loaded_from_cache = true;
            }
            else
            {
                std::cout << "Failed to load weighted graph from cache (empty or invalid). Will regenerate from hypergraph." << std::endl;
                delete test_graph_ptr;
                test_graph_ptr = nullptr;
                g1_ptr = nullptr;
            }
        }
        catch (const std::exception &e)
        {
            std::cout << "Error loading weighted graph from cache: " << e.what() << std::endl;
            std::cout << "Will regenerate from hypergraph." << std::endl;
            delete test_graph_ptr;
            test_graph_ptr = nullptr;
            g1_ptr = nullptr;
        }
    }

    if (!weighted_graph_loaded_from_cache)
    {
        std::cout << "Converting hypergraph to weighted graph..." << std::endl;

        std::cout << "--- BottleneckGraph (DLCR-like) Benchmark ---" << std::endl;
        std::cout << "Original hypergraph dataset: " << original_hypergraph_file << std::endl;

        try
        {
            std::cout << "Converting hypergraph to weighted graph" << std::endl;
            auto weightedGraph = hg_source.convertToWeightedGraph();
            std::cout << "Conversion complete. Saving weighted graph to: " << temp_weighted_graph_file << std::endl;
            weightedGraph->saveAdjList(temp_weighted_graph_file);
            std::cout << "Save complete." << std::endl;
        }
        catch (const std::exception &e)
        {
            std::filesystem::remove(temp_weighted_graph_file);
            FAIL() << "Failed to convert hypergraph to weighted graph: " << e.what();
            return;
        }
    }

    if (!weighted_graph_loaded_from_cache)
    {
        try
        {
            test_graph_ptr = new TestBottleneckGraph(temp_weighted_graph_file, true, false);
            if (!test_graph_ptr || !test_graph_ptr->g1)
            {
                std::filesystem::remove(temp_weighted_graph_file);
                FAIL() << "Cannot create TestBottleneckGraph or its internal g1 is empty. File: " << temp_weighted_graph_file;
                return;
            }
            g1_ptr = test_graph_ptr->g1;
        }
        catch (const std::exception &e)
        {
            std::filesystem::remove(temp_weighted_graph_file);
            FAIL() << "Failed to load newly generated weighted graph: " << temp_weighted_graph_file << " - " << e.what();
            return;
        }
    }
    else
    {
        ASSERT_NE(g1_ptr, nullptr) << "g1_ptr should not be null after cache loading.";
        std::cout << "--- BottleneckGraph (DLCR-like) Benchmark (using cached weighted graph) ---" << std::endl;
        std::cout << "Weighted graph file (cached): " << temp_weighted_graph_file << std::endl;
        std::cout << "K threshold used: " << k_threshold_for_conversion << std::endl;
    }

    if (g1_ptr->n == 0)
    {
        std::filesystem::remove(temp_weighted_graph_file);
        GTEST_SKIP() << "Converted weighted graph is empty (0 vertices), cannot run benchmark. File: " << temp_weighted_graph_file;
    }

    std::cout << "Weighted graph statistics (after conversion): " << g1_ptr->n << " vertices, " << g1_ptr->m << " edges" << std::endl;
    std::cout << "Edge addition operations (on weighted graph): " << num_add_edge_ops << std::endl;
    std::cout << "Edge removal operations (on weighted graph): " << num_remove_edge_ops << std::endl;

    std::cout << "\n[1. Index Construction]" << std::endl;
    auto build_start_time = high_resolution_clock::now();
    test_graph_ptr->TestConstruction();
    auto build_end_time = high_resolution_clock::now();
    auto build_duration_ms = duration_cast<milliseconds>(build_end_time - build_start_time);
    std::cout << "Index construction time (via TestConstruction): " << build_duration_ms.count() << " ms" << std::endl;

    std::cout << "\n[2. Index Size]" << std::endl;
    unsigned long long index_size_bytes = g1_ptr->GetIndexSize();
    double index_size_mb = static_cast<double>(index_size_bytes) / (1024.0 * 1024.0);
    std::cout << "Index size: " << std::fixed << std::setprecision(2) << index_size_mb << " MB" << std::endl;

    std::cout << "\n[3. Dynamic Update Performance]" << std::endl;
    long long total_add_op_time_nanoseconds = 0;
    int add_ops_performed_count = 0;
    long long total_remove_op_time_nanoseconds = 0;
    int remove_ops_performed_count = 0;
    long long total_rebuild_time_delete_ns = 0;
    long long total_added_pairs_count = 0;
    long long total_removed_pairs_count = 0;

    const int MAX_PARAM_GENERATION_ATTEMPTS_PER_OP = 20;
    const int ATTEMPT_MULTIPLIER = 5;
    const int MAX_HYPEREDGES_TO_SELECT_PER_OP = 100;

    if (num_add_edge_ops > 0 && g1_ptr->n > 1 && hg_source.numVertices() > 0 && hg_source.numHyperedges() >= 2)
    {
        int total_add_attempts_overall = 0;
        while (add_ops_performed_count < num_add_edge_ops && total_add_attempts_overall < num_add_edge_ops * ATTEMPT_MULTIPLIER)
        {
            total_add_attempts_overall++;

            std::uniform_int_distribution<> hg_vertex_dist(0, hg_source.numVertices() - 1);
            int hg_vertex_to_operate_on = hg_vertex_dist(random_generator);

            int num_hg_edges_to_select = std::min((int)hg_source.numHyperedges(), MAX_HYPEREDGES_TO_SELECT_PER_OP);
            if (num_hg_edges_to_select < 2)
            {
                continue;
            }
            std::uniform_int_distribution<> num_he_dist(2, num_hg_edges_to_select);
            int k_selected_hyperedges = num_he_dist(random_generator);

            std::vector<int> all_hyperedge_ids(hg_source.numHyperedges());
            std::iota(all_hyperedge_ids.begin(), all_hyperedge_ids.end(), 0);
            std::shuffle(all_hyperedge_ids.begin(), all_hyperedge_ids.end(), random_generator);

            std::vector<int> selected_hyperedge_ids;
            for (int i = 0; i < k_selected_hyperedges; ++i)
            {
                selected_hyperedge_ids.push_back(all_hyperedge_ids[i]);
            }

            long long current_op_affected_pairs = 0;
            auto update_start_time = high_resolution_clock::now();
            try
            {
                for (size_t i = 0; i < selected_hyperedge_ids.size(); ++i)
                {
                    for (size_t j = i + 1; j < selected_hyperedge_ids.size(); ++j)
                    {
                        int u_g1 = selected_hyperedge_ids[i] + 1;
                        int v_g1 = selected_hyperedge_ids[j] + 1;
                        if (u_g1 > 0 && u_g1 <= g1_ptr->n && v_g1 > 0 && v_g1 <= g1_ptr->n)
                        {
                            g1_ptr->DynamicAddEdge(u_g1, v_g1, 1);
                            current_op_affected_pairs++;
                        }
                    }
                }
                auto update_end_time = high_resolution_clock::now();
                if (current_op_affected_pairs > 0)
                {
                    total_add_op_time_nanoseconds += duration_cast<nanoseconds>(update_end_time - update_start_time).count();
                    add_ops_performed_count++;
                    total_added_pairs_count += current_op_affected_pairs;
                }
            }
            catch (const std::exception &e)
            {
                continue;
            }
        }
    }

    if (num_remove_edge_ops > 0 && g1_ptr->m > 0 && hg_source.numVertices() > 0 && hg_source.numHyperedges() >= 2)
    {
        g1_ptr->rebuild_time = 0;
        int total_remove_op_attempts = 0;

        while (remove_ops_performed_count < num_remove_edge_ops && total_remove_op_attempts < num_remove_edge_ops * ATTEMPT_MULTIPLIER)
        {
            total_remove_op_attempts++;

            std::uniform_int_distribution<> hg_vertex_dist(0, hg_source.numVertices() - 1);
            int hg_vertex_to_operate_on = hg_vertex_dist(random_generator);

            std::vector<int> incident_hyperedges = hg_source.getIncidentHyperedges(hg_vertex_to_operate_on);

            if (incident_hyperedges.size() < 2)
            {
                continue;
            }

            int num_hg_edges_to_select = std::min((int)incident_hyperedges.size(), MAX_HYPEREDGES_TO_SELECT_PER_OP);
            std::uniform_int_distribution<> num_he_dist(2, num_hg_edges_to_select);
            int k_selected_hyperedges = num_he_dist(random_generator);

            std::shuffle(incident_hyperedges.begin(), incident_hyperedges.end(), random_generator);
            std::vector<int> selected_hyperedge_ids_for_removal;
            for (int i = 0; i < k_selected_hyperedges; ++i)
            {
                selected_hyperedge_ids_for_removal.push_back(incident_hyperedges[i]);
            }

            long long current_op_affected_pairs = 0;
            auto op_start_time = high_resolution_clock::now();
            try
            {
                for (size_t i = 0; i < selected_hyperedge_ids_for_removal.size(); ++i)
                {
                    for (size_t j = i + 1; j < selected_hyperedge_ids_for_removal.size(); ++j)
                    {
                        int u_g1 = selected_hyperedge_ids_for_removal[i] + 1;
                        int v_g1 = selected_hyperedge_ids_for_removal[j] + 1;
                        if (u_g1 > 0 && u_g1 <= g1_ptr->n && v_g1 > 0 && v_g1 <= g1_ptr->n)
                        {
                            g1_ptr->DynamicDeleteEdge(u_g1, v_g1, 1);
                            current_op_affected_pairs++;
                        }
                    }
                }
                auto op_end_time = high_resolution_clock::now();
                if (current_op_affected_pairs > 0)
                {
                    total_remove_op_time_nanoseconds += duration_cast<nanoseconds>(op_end_time - op_start_time).count();
                    remove_ops_performed_count++;
                    total_removed_pairs_count += current_op_affected_pairs;
                }
            }
            catch (const std::exception &e)
            {
                continue;
            }
        }
        total_rebuild_time_delete_ns = g1_ptr->rebuild_time;
    }

    long long total_update_time_nanoseconds = total_add_op_time_nanoseconds + total_remove_op_time_nanoseconds;
    int total_update_ops_performed = add_ops_performed_count + remove_ops_performed_count;

    if (total_update_ops_performed > 0)
    {
        double avg_update_time_ns = static_cast<double>(total_update_time_nanoseconds) / total_update_ops_performed;
        std::cout << "Average update operation time (" << total_update_ops_performed << " total updates, "
                  << add_ops_performed_count << " additions, " << remove_ops_performed_count << " removals): "
                  << std::fixed << std::setprecision(2)
                  << avg_update_time_ns << " ns | "
                  << avg_update_time_ns / 1000.0 << " µs | "
                  << avg_update_time_ns / 1000000.0 << " ms" << std::endl;
        if (remove_ops_performed_count > 0)
        {
            double avg_rebuild_per_delete_op_ns = static_cast<double>(total_rebuild_time_delete_ns) / remove_ops_performed_count;
            std::cout << "  Average internal rebuild time (during delete operations): "
                      << std::fixed << std::setprecision(2)
                      << avg_rebuild_per_delete_op_ns << " ns | "
                      << avg_rebuild_per_delete_op_ns / 1000.0 << " µs | "
                      << avg_rebuild_per_delete_op_ns / 1000000.0 << " ms" << std::endl;
            if (total_removed_pairs_count > 0)
            {
                double avg_removed_pairs_per_op = static_cast<double>(total_removed_pairs_count) / remove_ops_performed_count;
                std::cout << "  Average hyperedge pairs affected per removal operation: "
                          << std::fixed << std::setprecision(2)
                          << avg_removed_pairs_per_op << std::endl;
            }
        }
        if (add_ops_performed_count > 0 && total_added_pairs_count > 0)
        {
            double avg_added_pairs_per_op = static_cast<double>(total_added_pairs_count) / add_ops_performed_count;
            std::cout << "  Average hyperedge pairs affected per addition operation: "
                      << std::fixed << std::setprecision(2)
                      << avg_added_pairs_per_op << std::endl;
        }
    }

    std::cout << "\n[4. Operation Mix Analysis]" << std::endl;
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

    auto print_projected_time_dlcr = [](const std::string &label, std::optional<double> time_ns)
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

    print_projected_time_dlcr("  Average time (100% additions)", avg_add_time_ns_opt);
    print_projected_time_dlcr("  Average time (100% deletions)", avg_remove_time_ns_opt);

    std::optional<double> mix_25A_75R_time_ns;
    if (avg_add_time_ns_opt && avg_remove_time_ns_opt)
    {
        mix_25A_75R_time_ns = (0.25 * *avg_add_time_ns_opt) + (0.75 * *avg_remove_time_ns_opt);
    }
    print_projected_time_dlcr("  Average time (25% add, 75% delete)", mix_25A_75R_time_ns);

    std::optional<double> mix_50A_50R_time_ns;
    if (avg_add_time_ns_opt && avg_remove_time_ns_opt)
    {
        mix_50A_50R_time_ns = (0.50 * *avg_add_time_ns_opt) + (0.50 * *avg_remove_time_ns_opt);
    }
    print_projected_time_dlcr("  Average time (50% add, 50% delete)", mix_50A_50R_time_ns);

    std::optional<double> mix_75A_25R_time_ns;
    if (avg_add_time_ns_opt && avg_remove_time_ns_opt)
    {
        mix_75A_25R_time_ns = (0.75 * *avg_add_time_ns_opt) + (0.25 * *avg_remove_time_ns_opt);
    }
    print_projected_time_dlcr("  Average time (75% add, 25% delete)", mix_75A_25R_time_ns);

    std::cout << "\n--- BottleneckGraph (DLCR-like) Benchmark End ---" << std::endl;
}

int main(int argc, char **argv)
{
    for (int i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];
        if (arg == "--dataset" && i + 1 < argc)
        {
            BenchmarkConfigDLCR::g_dataset_file = argv[++i];
        }
        else if (arg == "--k" && i + 1 < argc)
        {
            BenchmarkConfigDLCR::g_k_threshold_conversion = std::stoi(argv[++i]);
        }
        else if (arg == "--help")
        {
            std::cout << "BottleneckGraph Benchmark Options:\n"
                      << "  --dataset <filepath>      Path to the original hypergraph dataset file.\n"
                      << "  --k <threshold>           K threshold for hypergraph to weighted graph conversion (default: 1).\n"
                      << "  --help                    Show this help message.\n";
            return 0;
        }
    }

    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
