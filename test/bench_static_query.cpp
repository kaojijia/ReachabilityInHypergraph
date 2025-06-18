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
#include <set>
#include <filesystem>

using namespace std;
using namespace std::chrono;

class HypergraphTest : public ::testing::Test
{
public:
    static std::string hypergraph_file;

protected:
    void SetUp() override {}
};

std::string HypergraphTest::hypergraph_file = "";

TEST_F(HypergraphTest, PerformanceComparison)
{
    std::filesystem::path input_path(hypergraph_file);
    string cache_path = input_path.parent_path().string();
    if (!cache_path.empty() && cache_path.back() != '/')
    {
        cache_path += "/";
    }

    cout << "Cache path set to: " << cache_path << endl;

    Hypergraph hg;
    try
    {
        hg = Hypergraph::fromFile(hypergraph_file);
    }
    catch (const std::exception &e)
    {
        cerr << "Error loading hypergraph: " << e.what() << endl;
        FAIL() << "Failed to load hypergraph file: " << hypergraph_file;
        return;
    }

    cout << "Hypergraph loaded: " << hg.numVertices() << " vertices, "
         << hg.numHyperedges() << " hyperedges." << endl;

    if (hg.numVertices() == 0 || hg.numHyperedges() == 0)
    {
        cout << "Hypergraph is empty, skipping performance test." << endl;
        GTEST_SKIP() << "Skipping performance test on empty hypergraph.";
        return;
    }

    cout << "Building indices concurrently (Layered DS, UWeightedPLL, TreeIndex)..." << endl;

    std::chrono::milliseconds build_duration_baseline(0);
    std::chrono::milliseconds build_duration_pll(0);
    std::chrono::milliseconds build_duration_tree(0);
    std::atomic<bool> baseline_build_error(false);
    std::atomic<bool> pll_build_error(false);
    std::atomic<bool> tree_build_error(false);

    std::unique_ptr<HypergraphTreeIndex> tree_index;
    try
    {
        tree_index = std::make_unique<HypergraphTreeIndex>(hg);
    }
    catch (const std::exception &e)
    {
        cerr << "Error creating HypergraphTreeIndex object: " << e.what() << endl;
        FAIL() << "Failed to create HypergraphTreeIndex object.";
        return;
    }

    std::thread baseline_thread([&hg, &build_duration_baseline, &baseline_build_error, &cache_path]()
                                {
        auto start = high_resolution_clock::now();
        try {
             hg.offline_industry_baseline();
        } catch (const std::exception& e) {
             cerr << "Error in baseline build thread: " << e.what() << endl;
             baseline_build_error = true;
        }
        auto end = high_resolution_clock::now();
        build_duration_baseline = duration_cast<milliseconds>(end - start);
        cout << "Layered DS build completed." << endl; });

    std::thread pll_thread([&hg, &build_duration_pll, &pll_build_error, &cache_path]()
                           {
        auto start = high_resolution_clock::now();
         try {
            hg.offline_industry_pll();
         } catch (const std::exception& e) {
             cerr << "Error in PLL build thread: " << e.what() << endl;
             pll_build_error = true;
         }
        auto end = high_resolution_clock::now();
        build_duration_pll = duration_cast<milliseconds>(end - start);
        cout<< "PLL build completed." << endl; });

    HypergraphTreeIndex *tree_index_ptr = tree_index.get();
    std::thread tree_thread([tree_index_ptr, &build_duration_tree, &tree_build_error, &cache_path]()
                            {
         if (!tree_index_ptr) {
             cerr << "Error: tree_index_ptr is null in thread." << endl;
             tree_build_error = true;
             return;
         }
         auto start = high_resolution_clock::now();
         try {
             tree_index_ptr->buildIndexCacheSizeOnly();
         } catch (const std::exception& e) {
             cerr << "Error in TreeIndex build thread: " << e.what() << endl;
             tree_build_error = true;
         }
         auto end = high_resolution_clock::now();
         build_duration_tree = duration_cast<milliseconds>(end - start);
         cout << "TreeIndex build completed." << endl; });

    baseline_thread.join();
    pll_thread.join();
    tree_thread.join();

    if (baseline_build_error || pll_build_error || tree_build_error)
    {
        FAIL() << "Error occurred during concurrent index building. Check logs.";
        return;
    }

    cout << "Index building finished." << endl;
    cout << " - Layered DS build time:          " << build_duration_baseline.count() << " ms." << endl;
    cout << " - UWeightedPLL build time:        " << build_duration_pll.count() << " ms." << endl;
    cout << " - HypergraphTreeIndex build time: " << build_duration_tree.count() << " ms." << endl;

    cout << "\n--- Estimated Memory Usage ---" << endl;
    cout << fixed << setprecision(2);
    cout << "2. Layered DS Index (Total):        " << setw(8) << hg.getWeightedGraphsMemoryUsageMB() << " MB" << endl;
    cout << "3. UWeightedPLL Index:            " << setw(8) << hg.getPllMemoryUsageMB() << " MB" << endl;

    double tree_index_mem_mb = 0.0;
    if (tree_index)
    {
        tree_index_mem_mb = tree_index->getMemoryUsageMB();
    }
    cout << "4. HypergraphTreeIndex:           " << setw(8) << tree_index_mem_mb << " MB" << endl;

    const int num_queries_per_k = 5000;
    const int target_k_values[] = {2, 4, 6, 8, 10};
    const int num_k_values = sizeof(target_k_values) / sizeof(target_k_values[0]);
    const int total_queries = num_k_values * num_queries_per_k;

    std::vector<std::tuple<int, int, int, int>> queries;
    queries.reserve(total_queries);
    std::mt19937 rng(std::random_device{}());
    std::uniform_int_distribution<int> vertex_dist(0, hg.numVertices() - 1);

    cout << "\nGenerating " << total_queries << " random queries ("
         << num_queries_per_k << " per k value)..." << endl;

    for (int target_k : target_k_values)
    {
        for (int i = 0; i < num_queries_per_k; ++i)
        {
            int u = vertex_dist(rng);
            int v = vertex_dist(rng);
            queries.emplace_back(u, v, target_k, target_k);
        }
    }

    cout << "Generated " << queries.size() << " queries." << endl;

    cout << "\nRunning " << queries.size() << " queries..." << endl;
    std::map<std::string, std::map<int, std::pair<long long, int>>> timings;
    std::vector<std::string> method_names = {"BFS", "LayeredDS", "UWeightedPLL", "TreeIndex"};

    for (const auto &method : method_names)
    {
        for (int k : target_k_values)
        {
            timings[method][k] = {0, 0};
        }
    }

    cout << "  Executing BFS..." << endl;
    for (const auto &q_tuple : queries)
    {
        int u, v, k_query, k_group;
        std::tie(u, v, k_query, k_group) = q_tuple;
        auto start_ns = high_resolution_clock::now();
        try
        {
            hg.isReachableBidirectionalBFSByEdge(u, v, k_query);
        }
        catch (...)
        {
        }
        auto end_ns = high_resolution_clock::now();
        long long delta_ns = duration_cast<nanoseconds>(end_ns - start_ns).count();
        timings["BFS"][k_group].first += delta_ns;
        timings["BFS"][k_group].second++;
    }

    cout << "  Executing LayeredDS..." << endl;
    for (const auto &q_tuple : queries)
    {
        int u, v, k_query, k_group;
        std::tie(u, v, k_query, k_group) = q_tuple;
        auto start_ns = high_resolution_clock::now();
        try
        {
            hg.isReachableViaWeightedGraphByEdge(u, v, k_query);
        }
        catch (...)
        {
        }
        auto end_ns = high_resolution_clock::now();
        long long delta_ns = duration_cast<nanoseconds>(end_ns - start_ns).count();
        timings["LayeredDS"][k_group].first += delta_ns;
        timings["LayeredDS"][k_group].second++;
    }

    cout << "  Executing UWeightedPLL..." << endl;
    for (const auto &q_tuple : queries)
    {
        int u, v, k_query, k_group;
        std::tie(u, v, k_query, k_group) = q_tuple;
        auto start_ns = high_resolution_clock::now();
        try
        {
            hg.isReachableViaUWeightedPLLByEdge(u, v, k_query);
        }
        catch (...)
        {
        }
        auto end_ns = high_resolution_clock::now();
        long long delta_ns = duration_cast<nanoseconds>(end_ns - start_ns).count();
        timings["UWeightedPLL"][k_group].first += delta_ns;
        timings["UWeightedPLL"][k_group].second++;
    }

    cout << "  Executing TreeIndex..." << endl;
    for (const auto &q_tuple : queries)
    {
        int u, v, k_query, k_group;
        std::tie(u, v, k_query, k_group) = q_tuple;
        auto start_ns = high_resolution_clock::now();
        try
        {
            tree_index->queryByEdgeId(u, v, k_query);
        }
        catch (...)
        {
        }
        auto end_ns = high_resolution_clock::now();
        long long delta_ns = duration_cast<nanoseconds>(end_ns - start_ns).count();
        timings["TreeIndex"][k_group].first += delta_ns;
        timings["TreeIndex"][k_group].second++;
    }

    cout << "\n--- Query Performance Results ---" << endl;
    cout << fixed << setprecision(3);
    cout << "  k   | Method          | Avg. Time (µs) | Count |" << endl;
    cout << "------|-----------------|----------------|-------|" << endl;

    for (int k : target_k_values)
    {
        for (const auto &method : method_names)
        {
            long long total_ns = timings[method][k].first;
            int count = timings[method][k].second;
            double avg_us = (count == 0) ? 0.0 : (double)total_ns / count / 1000.0;

            cout << "  " << left << setw(3) << k
                 << " | " << left << setw(15) << method
                 << " | " << right << setw(14) << avg_us
                 << " | " << right << setw(5) << count << endl;
        }
        if (k != target_k_values[num_k_values - 1])
        {
            cout << "------|-----------------|----------------|-------|" << endl;
        }
    }
    cout << "\n--- Performance Comparison Test End ---" << endl;
}

int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    for (int i = 1; i < argc; i++)
    {
        if (std::string(argv[i]) == "--dataset" && i + 1 < argc)
        {
            HypergraphTest::hypergraph_file = argv[i + 1];
            cout << "Input file set to: " << HypergraphTest::hypergraph_file << endl;
            break;
        }
    }
    return RUN_ALL_TESTS();
}
