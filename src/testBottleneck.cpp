#include "testBottleneck.h"
#include <iostream>
#include <algorithm>
#include <chrono>
TestBottleneckGraph::TestBottleneckGraph(const std::string &filePath, bool useOrder, bool loadBinary)
    : filePath(filePath), useOrder(useOrder), loadBinary(loadBinary), g1(nullptr)
{
    timer.StartTimer("total_bottleneck_test_time");
    g1 = new bottleneckPath::BottleneckGraph(filePath, useOrder, loadBinary);
}

TestBottleneckGraph::~TestBottleneckGraph()
{
    delete g1;
    timer.EndTimerAndPrint("total_bottleneck_test_time");
}

void TestBottleneckGraph::TestConstruction()
{
    printf("\n=========== Start TestBottleneckGraph Construction ===========\n");
    timer.StartTimer("BottleneckGraph_Initialization");
    if (!g1)
    {
        printf("Failed to create BottleneckGraph object.\n");
        printf("=========== End TestBottleneckGraph Construction ===========\n");
        return;
    }
    timer.EndTimerAndPrint("BottleneckGraph_Initialization");
    printf("Graph initialization: OK\n");

    timer.StartTimer("BottleneckGraph_IndexConstruction");
    g1->ConstructIndex();
    timer.EndTimerAndPrint("BottleneckGraph_IndexConstruction");
    printf("Graph index construction: OK\n");

    g1->PrintStat();
    printf("=========== End TestBottleneckGraph Construction ===========\n");
}

void TestBottleneckGraph::TestQueryCorrectness(int num_queries)
{
    if (!g1)
    {
        printf("Graph not constructed. Call TestConstruction first.\n");
        TestConstruction();
    }
    if (!g1)
    {
        printf("Failed to construct graph for TestQueryCorrectness.\n");
        return;
    }
    if (g1->n == 0)
    {
        printf("Failed to construct graph or graph is empty for TestQueryCorrectness.\n");
        return;
    }
    g1->QueryBFS(3, 6, 1);

    printf("\n=========== Start TestBottleneckGraph Query Correctness (%d queries) ===========\n", num_queries);

    std::default_random_engine e(time(nullptr));

    std::uniform_int_distribution<int> vertex_dist(1, g1->n);

    int min_w = 1, max_w = 100;

    std::uniform_int_distribution<int> weight_dist(min_w, max_w);

    int correct_count = 0;
    long long total_query_time_ns = 0;
    long long total_bfs_time_ns = 0;

    for (int i = 0; i < num_queries; ++i)
    {
        int s = vertex_dist(e);
        int t = vertex_dist(e);
        int query_w = weight_dist(e);

        timer.StartTimer("Bottleneck_Query_2Hop");
        bool res_2hop = g1->Query(s, t, query_w);
        total_query_time_ns += timer.EndTimer("Bottleneck_Query_2Hop");

        timer.StartTimer("Bottleneck_Query_BFS");
        bool res_bfs = g1->QueryBFS(s, t, query_w);
        total_bfs_time_ns += timer.EndTimer("Bottleneck_Query_BFS");

        if (res_2hop == res_bfs)
        {
            correct_count++;
        }
        else
        {
            printf("Mismatch! s=%d, t=%d, query_w=%d. 2Hop: %d, BFS: %d\n",
                   s, t, query_w, res_2hop, res_bfs);
        }
    }

    double avg_2hop_time = num_queries > 0 ? (double)total_query_time_ns / num_queries : 0.0;
    double avg_bfs_time = num_queries > 0 ? (double)total_bfs_time_ns / num_queries : 0.0;

    printf("=========== Query Correctness Summary ===========\n");
    printf("总查询数: %d\n", num_queries);
    printf("2-Hop与单向BFS结果一致数: %d\n", correct_count);
    printf("一致率: %.2f%%\n", num_queries > 0 ? 100.0 * correct_count / num_queries : 0.0);
    printf("平均2-Hop查询耗时: %.2f ns\n", avg_2hop_time);
    printf("平均单向BFS查询耗时: %.2f ns\n", avg_bfs_time);
    if (correct_count == num_queries)
    {
        printf("所有查询结果完全一致，正确性验证通过。\n");
    }
    else
    {
        printf("存在不一致的查询结果，请检查实现。\n");
    }
    printf("=========== End of Query Correctness Summary ===========\n");
    printf("=========== End TestBottleneckGraph Query Correctness ===========\n");
}

void TestBottleneckGraph::TestDynamicAddAndQuery(int add_u, int add_v, int add_w, int num_queries_after_add)
{
    if (!g1)
    {
        printf("Graph not constructed. Call TestConstruction first.\n");
        TestConstruction();
    }
    if (!g1)
    {
        printf("Failed to construct graph for TestDynamicAddAndQuery.\n");
        return;
    }

    printf("\n=========== Start TestBottleneckGraph Dynamic Add Edge (%d, %d, %d) ===========\n", add_u, add_v, add_w);

    printf("--- Labels before adding edge (%d, %d, %d) ---\n", add_u, add_v, add_w);
    if (add_u > 0 && add_u <= g1->n)
        g1->PrintNodeLabel(add_u);
    if (add_v > 0 && add_v <= g1->n && add_u != add_v)
        g1->PrintNodeLabel(add_v);

    timer.StartTimer("BottleneckGraph_DynamicAddEdge");
    g1->DynamicAddEdge(add_u, add_v, add_w);
    timer.EndTimerAndPrint("BottleneckGraph_DynamicAddEdge");
    printf("Edge (%d, %d, %d) added.\n", add_u, add_v, add_w);

    g1->PrintStat();

    printf("--- Labels after adding edge (%d, %d, %d) ---\n", add_u, add_v, add_w);
    if (add_u > 0 && add_u <= g1->n)
        g1->PrintNodeLabel(add_u);
    if (add_v > 0 && add_v <= g1->n && add_u != add_v)
        g1->PrintNodeLabel(add_v);

    printf("--- Running general query correctness test after edge addition ---\n");
    TestQueryCorrectness(num_queries_after_add > 0 ? num_queries_after_add : DEFAULT_BOTTLENECK_TEST_NUM);

    printf("=========== End TestBottleneckGraph Dynamic Add Edge ===========\n");
}

void TestBottleneckGraph::TestBatchDynamicAddAndQuery(const std::vector<std::tuple<int, int, int>> &edges_to_add, int num_queries_after_batch_add)
{
    if (!g1 || g1->n == 0)
    {
        std::cerr << "Graph not constructed or empty. Run TestConstruction first." << std::endl;
        return;
    }
    if (edges_to_add.empty())
    {
        std::cout << "No edges provided for batch dynamic add. Skipping." << std::endl;
        TestQueryCorrectness(num_queries_after_batch_add);
        return;
    }

    std::cout << "--- Bottleneck Batch Dynamic Add Edges & Query Test (" << edges_to_add.size() << " edges) ---" << std::endl;

    timer.StartTimer("BottleneckGraph_BatchDynamicAddEdges");
    for (const auto &edge_tuple : edges_to_add)
    {
        int u = std::get<0>(edge_tuple);
        int v = std::get<1>(edge_tuple);
        int w = std::get<2>(edge_tuple);
        std::cout << "Adding edge: (" << u << ", " << v << ", w: " << w << ")" << std::endl;
        g1->DynamicAddEdge(u, v, w);
    }
    timer.EndTimer("BottleneckGraph_BatchDynamicAddEdges");

    g1->PrintStat();

    TestQueryCorrectness(num_queries_after_batch_add);
    std::cout << "--- End Bottleneck Batch Dynamic Add Edges & Query Test ---" << std::endl;
}

void TestBottleneckGraph::TestDynamicDeleteAndQuery(int u, int v, int weight_change)
{
    if (!g1)
    {
        printf("Error: Graph not initialized. Call TestConstruction first.\n");
        return;
    }

    printf("\n--- Testing Dynamic Delete/Reduce for Edge (%d,%d) with weight change %d ---\n", u, v, weight_change);

    int original_weight = 0;
    bool edge_exists = false;
    for (const auto &edge_pair : g1->Adj[u])
    {
        if (edge_pair.first == v)
        {
            edge_exists = true;
            original_weight = edge_pair.second;
            break;
        }
    }

    if (!edge_exists)
    {
        printf("Error: Edge (%d,%d) does not exist in the graph. Cannot delete/reduce.\n", u, v);
        return;
    }

    int expected_new_weight = std::max(0, original_weight + weight_change);
    bool complete_removal = (expected_new_weight <= 0);

    std::vector<std::tuple<int, int, int, bool>> pre_delete_queries;
    std::default_random_engine e(time(nullptr));
    std::uniform_int_distribution<int> vertex_dist(1, g1->n);
    std::uniform_int_distribution<int> weight_dist(1, 100);

    for (int i = 0; i < 20; ++i)
    {
        int s = vertex_dist(e);
        int t = vertex_dist(e);
        int w = weight_dist(e);
        bool result_bfs = g1->QueryBFS(s, t, w);
        pre_delete_queries.emplace_back(s, t, w, result_bfs);
    }

    pre_delete_queries.emplace_back(u, v, original_weight, true);
    if (!complete_removal)
    {
        pre_delete_queries.emplace_back(u, v, expected_new_weight, true);
    }
    pre_delete_queries.emplace_back(u, v, original_weight + 1, false);

    auto start_time = std::chrono::high_resolution_clock::now();
    g1->DynamicDeleteEdge(u, v, -weight_change);
    auto end_time = std::chrono::high_resolution_clock::now();
    auto elapsed_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count();

    printf("Edge %s in %lld nanoseconds (%.6f ms)\n",
           complete_removal ? "deleted" : "weight reduced",
           elapsed_ns, elapsed_ns / 1000000.0);

    bool all_correct = true;
    int errors = 0;

    long long bfs_ns_total = 0, pl_ns_total = 0;
    int bfs_count = 0, pl_count = 0;

    for (const auto &[s, t, w, expected_result_before] : pre_delete_queries)
    {
        bool expected_result_after = expected_result_before;

        if ((s == u && t == v) || (s == v && t == u))
        {

            if (w > expected_new_weight && expected_result_before)
            {
                expected_result_after = false;
            }
        }

        start_time = std::chrono::high_resolution_clock::now();
        bool pl_result = g1->Query(s, t, w);
        end_time = std::chrono::high_resolution_clock::now();
        pl_ns_total += std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count();
        pl_count++;

        start_time = std::chrono::high_resolution_clock::now();
        bool bfs_result = g1->QueryBFS(s, t, w);
        end_time = std::chrono::high_resolution_clock::now();
        bfs_ns_total += std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count();
        bfs_count++;

        if (pl_result != bfs_result)
        {
            printf("Query (%d,%d,%d) - PLL and BFS results differ: PLL=%s, BFS=%s\n",
                   s, t, w, pl_result ? "true" : "false", bfs_result ? "true" : "false");
            all_correct = false;
            errors++;
        }

        if (bfs_result != expected_result_after)
        {

            if (!((s == u && t == v) || (s == v && t == u)))
            {
                printf("Query (%d,%d,%d) - Expected %s after modification, but got %s\n",
                       s, t, w, expected_result_after ? "true" : "false",
                       bfs_result ? "true" : "false");
            }
        }
    }

    if (complete_removal)
    {

        bool direct_query = g1->Query(u, v, 1);
        bool direct_bfs = g1->QueryBFS(u, v, 1);
        if (direct_query != direct_bfs)
        {
            printf("Direct edge query (%d,%d,1) after removal: PLL=%s, BFS=%s\n",
                   u, v, direct_query ? "true" : "false", direct_bfs ? "true" : "false");
            all_correct = false;
            errors++;
        }
    }
    else
    {

        bool below_query = g1->Query(u, v, expected_new_weight);
        bool below_bfs = g1->QueryBFS(u, v, expected_new_weight);
        bool above_query = g1->Query(u, v, original_weight);
        bool above_bfs = g1->QueryBFS(u, v, original_weight);

        if (below_query != below_bfs || !below_bfs)
        {
            printf("Query below (%d,%d,%d) after reduction: PLL=%s, BFS=%s, expected true\n",
                   u, v, expected_new_weight, below_query ? "true" : "false",
                   below_bfs ? "true" : "false");
            all_correct = false;
            errors++;
        }

        if (above_query != above_bfs || above_bfs)
        {
            printf("Query above (%d,%d,%d) after reduction: PLL=%s, BFS=%s, expected false\n",
                   u, v, original_weight, above_query ? "true" : "false",
                   above_bfs ? "true" : "false");
            all_correct = false;
            errors++;
        }
    }

    if (all_correct)
    {
        printf("All test queries after edge modification are consistent between PLL and BFS! ✓\n");
    }
    else
    {
        printf("%d errors found in query verification. Please check the implementation.\n", errors);
    }

    printf("PLL Query average: %lld nanoseconds (%.6f ms) over %d queries\n",
           pl_count > 0 ? pl_ns_total / pl_count : 0,
           pl_count > 0 ? (pl_ns_total / pl_count) / 1000000.0 : 0,
           pl_count);

    printf("BFS Query average: %lld nanoseconds (%.6f ms) over %d queries\n",
           bfs_count > 0 ? bfs_ns_total / bfs_count : 0,
           bfs_count > 0 ? (bfs_ns_total / bfs_count) / 1000000.0 : 0,
           bfs_count);

    printf("PLL/BFS speedup ratio: %.2fx\n",
           (bfs_count > 0 && pl_count > 0) ? (double)(bfs_ns_total / bfs_count) / (pl_ns_total / pl_count) : 0);

    printf("--- End of Dynamic Delete Test ---\n\n");
}

void TestBottleneckGraph::TestBatchDynamicDeleteAndQuery(
    const std::vector<std::tuple<int, int, int>> &edges_to_modify, int num_queries)
{

    if (!g1)
    {
        printf("Error: Graph not initialized. Call TestConstruction first.\n");
        return;
    }

    printf("\n--- Testing Batch Dynamic Delete/Reduce for %zu edges ---\n", edges_to_modify.size());

    int original_edge_count = g1->m;
    printf("Graph before modifications: %d nodes, %lld edges\n", g1->n, g1->m);

    auto start_time = std::chrono::high_resolution_clock::now();

    int complete_deletions = 0;
    int weight_reductions = 0;

    for (const auto &[u, v, weight_change] : edges_to_modify)
    {
        if (weight_change >= 0)
        {
            printf("Warning: Positive weight change %d for edge (%d,%d) - expected negative for deletion.\n",
                   weight_change, u, v);
            continue;
        }

        int current_weight = 0;
        bool edge_exists = false;
        for (const auto &edge_pair : g1->Adj[u])
        {
            if (edge_pair.first == v)
            {
                edge_exists = true;
                current_weight = edge_pair.second;
                break;
            }
        }

        if (!edge_exists)
        {
            printf("Warning: Edge (%d,%d) does not exist, skipping.\n", u, v);
            continue;
        }

        int actual_weight = -weight_change;
        if (actual_weight >= current_weight)
        {
            complete_deletions++;
        }
        else
        {
            weight_reductions++;
        }

        g1->DynamicDeleteEdge(u, v, actual_weight);
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    auto elapsed_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count();
    auto avg_ns_per_op = (complete_deletions + weight_reductions > 0) ? elapsed_ns / (complete_deletions + weight_reductions) : 0;

    printf("Batch modification completed in %lld nanoseconds (%.6f ms)\n",
           elapsed_ns, elapsed_ns / 1000000.0);
    printf("Average time per edge modification: %lld nanoseconds (%.6f ms)\n",
           avg_ns_per_op, avg_ns_per_op / 1000000.0);
    printf("Modifications: %d complete deletions, %d weight reductions.\n",
           complete_deletions, weight_reductions);
    printf("Graph after modifications: %d nodes, %lld edges\n", g1->n, g1->m);

    unsigned long long label_count = 0;
    for (int i = 1; i <= g1->n; ++i)
    {
        if (i < g1->Label.size())
        {
            label_count += g1->Label[i].size();
        }
    }
    printf("Total label entries after modifications: %llu\n", label_count);

    printf("\nPerforming %d random queries to verify index correctness...\n", num_queries);

    std::default_random_engine e(time(nullptr));
    std::uniform_int_distribution<int> vertex_dist(1, g1->n);
    std::uniform_int_distribution<int> weight_dist(1, 100);

    int errors = 0;
    long long bfs_ns_total = 0, pl_ns_total = 0;

    start_time = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < num_queries; ++i)
    {
        int s = vertex_dist(e);
        int t = vertex_dist(e);
        int w = weight_dist(e);

        auto query_start = std::chrono::high_resolution_clock::now();
        bool pl_result = g1->Query(s, t, w);
        auto query_end = std::chrono::high_resolution_clock::now();
        pl_ns_total += std::chrono::duration_cast<std::chrono::nanoseconds>(query_end - query_start).count();

        query_start = std::chrono::high_resolution_clock::now();
        bool bfs_result = g1->QueryBFS(s, t, w);
        query_end = std::chrono::high_resolution_clock::now();
        bfs_ns_total += std::chrono::duration_cast<std::chrono::nanoseconds>(query_end - query_start).count();

        if (pl_result != bfs_result)
        {
            printf("Query (%d,%d,%d) - ERROR: PLL=%s, BFS=%s\n",
                   s, t, w, pl_result ? "true" : "false", bfs_result ? "true" : "false");
            errors++;

            if (errors > 20)
            {
                printf("Too many errors, stopping verification.\n");
                break;
            }
        }
    }

    end_time = std::chrono::high_resolution_clock::now();
    auto query_elapsed_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count();

    if (errors == 0)
    {
        printf("All %d queries correctly match BFS results! ✓\n", num_queries);
    }
    else
    {
        printf("%d errors found in %d queries. Please check the implementation.\n", errors, num_queries);
    }

    double avg_pl_ns = num_queries > 0 ? (double)pl_ns_total / num_queries : 0;
    double avg_bfs_ns = num_queries > 0 ? (double)bfs_ns_total / num_queries : 0;
    double speedup = avg_bfs_ns > 0 ? avg_bfs_ns / avg_pl_ns : 0;

    printf("Query verification completed in %lld nanoseconds (%.6f ms)\n",
           query_elapsed_ns, query_elapsed_ns / 1000000.0);
    printf("PLL Query average: %.2f nanoseconds (%.6f ms)\n", avg_pl_ns, avg_pl_ns / 1000000.0);
    printf("BFS Query average: %.2f nanoseconds (%.6f ms)\n", avg_bfs_ns, avg_bfs_ns / 1000000.0);
    printf("PLL/BFS speedup ratio: %.2fx\n", speedup);

    double pl_qps = 1000000000.0 / (avg_pl_ns > 0 ? avg_pl_ns : 1);
    printf("PLL Query QPS: %.2f queries/second\n", pl_qps);

    printf("--- End of Batch Dynamic Delete Test ---\n\n");
}