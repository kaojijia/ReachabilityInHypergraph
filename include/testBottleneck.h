#pragma once

#include "bottleneckGraph.h"
#include "timer.h"
#include <string>
#include <vector>
#include <tuple>
#include <random>
#include <ctime>

const int DEFAULT_BOTTLENECK_TEST_NUM = 100;

class TestBottleneckGraph
{
public:
    std::string filePath;
    bool useOrder;
    bool loadBinary;

    bottleneckPath::BottleneckGraph *g1;

    cx::Timer timer;

    TestBottleneckGraph(const std::string &filePath, bool useOrder = true, bool loadBinary = false);
    ~TestBottleneckGraph();

    void TestConstruction();
    void TestQueryCorrectness(int num_queries = DEFAULT_BOTTLENECK_TEST_NUM);

    void TestBatchDynamicAddAndQuery(const std::vector<std::tuple<int, int, int>> &edges_to_add, int num_queries_after_batch_add = DEFAULT_BOTTLENECK_TEST_NUM);
    void TestDynamicAddAndQuery(int add_u, int add_v, int add_w, int num_queries_after_add);
    void TestDynamicDeleteAndQuery(int u, int v, int weight_change);
    void TestBatchDynamicDeleteAndQuery(const std::vector<std::tuple<int, int, int>> &edges_to_modify, int num_queries);
};
