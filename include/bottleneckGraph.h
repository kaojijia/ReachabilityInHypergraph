#pragma once

#include <cstdio>
#include <cstdlib>
#include <string>
#include <algorithm>
#include <set>
#include <utility>
#include <vector>
#include <random>
#include <iostream>
#include <queue>
#include <map>
#include <limits>
#include "timer.h"

namespace bottleneckPath
{

    const int INF_WEIGHT = std::numeric_limits<int>::max();
    const int INVALID_PREDECESSOR = 0;

    struct EdgeNode
    {
        int s;
        int t;
        int weight;

        int getOtherEndpoint(int u) const
        {
            if (s == u)
                return t;
            if (t == u)
                return s;
            return -1;
        }
    };

    struct BottleneckLabelEntry
    {
        int hub_id;
        int bottleneck_weight;
        int predecessor_on_path_to_hub;

        bool operator<(const BottleneckLabelEntry &other) const
        {
            if (hub_id != other.hub_id)
            {
                return hub_id < other.hub_id;
            }
            if (bottleneck_weight != other.bottleneck_weight)
            {
                return bottleneck_weight > other.bottleneck_weight;
            }
            return predecessor_on_path_to_hub < other.predecessor_on_path_to_hub;
        }
    };

    struct DegreeNode
    {
        int id;
        int num;

        DegreeNode(int i = -1, int n = 0) : id(i), num(n) {}
    };

    typedef std::map<int, BottleneckLabelEntry> HubBottleneckMap;

    class BottleneckGraph
    {
    public:
        long long rebuild_time;
        cx::Timer timer;

        int n;
        long long m;

        std::vector<std::vector<std::pair<int, int>>> Adj;

        std::vector<DegreeNode> degreeList;
        std::vector<int> rankList;

        std::string index_file_path;

        std::vector<HubBottleneckMap> Label;

        BottleneckGraph(const std::string &filePath, bool useOrder, bool loadBinary);
        ~BottleneckGraph();

        void ConstructIndex();

        bool SaveIndexToFile() const;

        bool LoadIndexFromFile();

        bool Query(int s, int t, int query_weight);
        bool QueryBFS(int s, int t, int query_weight);

        void DynamicAddEdge(int u, int v, int weight);
        void DynamicDeleteEdge(int u, int v, int weight);
        void DynamicBatchAdd(const std::vector<std::tuple<int, int, int>> &addedEdges);
        void DynamicBatchDelete(const std::vector<std::tuple<int, int, int>> &deletedEdges);

        unsigned long long GetIndexSize();
        void PrintStat();
        void PrintNodeLabel(int node_id);

        void SpreadHubLabelUpdate(int start_node, int hub_id, int current_bottleneck, int predecessor);

        void RunComprehensiveLabelUpdate(int u, int v, int current_edge_weight);

    private:
        void ReadGraph(const std::string &filePath);
        void CalculateDegreesAndRank(bool useOrder);

        static bool CmpDegreeNode(const DegreeNode &a, const DegreeNode &b)
        {
            return a.num > b.num;
        }

        void BFSForIndex(int hub_node_id);

        void BFSForIndexOptimized(int hub_node_id);

        bool TryInsertLabel(HubBottleneckMap &labels, int hub_id, int bottleneck_w, int predecessor_node);

        void RunIncrementalBFSUpdate(int source_node, int first_target_node, int bn_to_first_target, int pred_of_first_target);

        void GenerateIndexFilePath(const std::string &graph_file_path);
    };

}
