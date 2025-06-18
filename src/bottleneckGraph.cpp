#include "bottleneckGraph.h"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <vector>
#include <queue>
#include <set>
#include <map>
#include <chrono>

namespace bottleneckPath
{

    const long long DYNAMIC_EDGE_ADDITION_BUFFER_VECTOR = 100000;

    BottleneckGraph::BottleneckGraph(const std::string &filePath, bool useOrder, bool loadBinary) : n(0), m(0)
    {

        if (loadBinary)
        {
            std::cerr << "Binary load not implemented for BottleneckGraph. Loading from text." << std::endl;
        }
        ReadGraph(filePath);
        CalculateDegreesAndRank(useOrder);
        rebuild_time = 0;
        Label.resize(n + 1);

        GenerateIndexFilePath(filePath);
    }

    BottleneckGraph::~BottleneckGraph()
    {
    }

    void BottleneckGraph::ReadGraph(const std::string &filePath)
    {
        std::ifstream infile(filePath);
        if (!infile.is_open())
        {
            std::cerr << "Error opening file: " << filePath << std::endl;
            exit(EXIT_FAILURE);
        }

        long long M_from_file;
        infile >> this->n >> M_from_file;

        Adj.assign(n + 1, std::vector<std::pair<int, int>>());

        std::map<std::pair<int, int>, int> unique_undirected_edges_weights;
        long long actual_edge_count = 0;

        for (long long i = 0; i < M_from_file; ++i)
        {
            int u, v, w;
            infile >> u >> v >> w;
            if (u > n || v > n || u < 0 || v < 0)
            {
                std::cerr << "Warning: Vertex out of bounds " << u << " or " << v << " (n=" << n << "). Skipping edge." << std::endl;
                continue;
            }
            if (u == v)
            {
                std::cerr << "Warning: Self-loop for vertex " << u << " with weight " << w << ". Skipping." << std::endl;
                continue;
            }

            int u_norm = std::min(u, v);
            int v_norm = std::max(u, v);

            auto it = unique_undirected_edges_weights.find({u_norm, v_norm});
            if (it == unique_undirected_edges_weights.end())
            {

                Adj[u].push_back({v, w});
                Adj[v].push_back({u, w});
                unique_undirected_edges_weights[{u_norm, v_norm}] = w;
                actual_edge_count++;
            }
            else
            {

#ifdef DEBUG
                std::cerr << "Warning: Duplicate undirected edge definition for (" << u_norm << ", " << v_norm
                          << ") in input file. Original weight: " << it->second
                          << ", Duplicate attempt with (" << u << "," << v << " w:" << w << "). Using first encountered." << std::endl;
#endif
            }
        }
        this->m = actual_edge_count;
        infile.close();
    }

    void BottleneckGraph::CalculateDegreesAndRank(bool useOrder)
    {
        degreeList.resize(n + 1);
        rankList.resize(n + 1, n);

        for (int i = 1; i <= n; ++i)
        {
            degreeList[i].id = i;
            degreeList[i].num = Adj[i].size();
        }

        if (useOrder)
        {
            std::vector<DegreeNode> sorted_degrees;
            for (int i = 1; i <= n; ++i)
                sorted_degrees.push_back(degreeList[i]);

            std::sort(sorted_degrees.begin(), sorted_degrees.end(), CmpDegreeNode);

            for (int i = 0; i < n; ++i)
            {
                rankList[sorted_degrees[i].id] = i;
            }
        }
        else
        {
            for (int i = 1; i <= n; ++i)
            {
                rankList[i] = i - 1;
            }
        }
    }

    bool BottleneckGraph::TryInsertLabel(HubBottleneckMap &current_labels, int hub_id, int bottleneck_w, int predecessor_node)
    {
        auto it = current_labels.find(hub_id);
        if (it != current_labels.end())
        {
            if (it->second.bottleneck_weight >= bottleneck_w)
            {
                return false;
            }
        }
        current_labels[hub_id] = {hub_id, bottleneck_w, predecessor_node};
        return true;
    }

    void BottleneckGraph::BFSForIndex(int hub_node_id)
    {
        std::queue<std::pair<int, int>> q;
        std::vector<int> visited_bottlenecks(n + 1, 0);
        std::vector<int> path_predecessor_nodes(n + 1, INVALID_PREDECESSOR);

        q.push({hub_node_id, INF_WEIGHT});
        visited_bottlenecks[hub_node_id] = INF_WEIGHT;

        while (!q.empty())
        {
            int u = q.front().first;
            int bn_hub_to_u = q.front().second;
            q.pop();

            if (rankList[u] > rankList[hub_node_id])
            {
                TryInsertLabel(Label[u], hub_node_id, bn_hub_to_u, path_predecessor_nodes[u]);
            }

            for (const auto &edge_pair : Adj[u])
            {
                int v = edge_pair.first;
                int weight = edge_pair.second;

                int new_bn_hub_to_v = std::min(bn_hub_to_u, weight);

                if (new_bn_hub_to_v > visited_bottlenecks[v])
                {
                    visited_bottlenecks[v] = new_bn_hub_to_v;
                    path_predecessor_nodes[v] = u;
                    q.push({v, new_bn_hub_to_v});
                }
            }
        }
    }

    void BottleneckGraph::BFSForIndexOptimized(int hub_node_id)
    {
        std::priority_queue<std::pair<int, int>> pq;
        std::vector<int> best_bottleneck(n + 1, 0);
        std::vector<int> path_predecessor_nodes(n + 1, INVALID_PREDECESSOR);

        pq.push({INF_WEIGHT, hub_node_id});
        best_bottleneck[hub_node_id] = INF_WEIGHT;

        while (!pq.empty())
        {
            auto [current_bn, u] = pq.top();
            pq.pop();

            if (current_bn < best_bottleneck[u])
                continue;

            if (rankList[u] > rankList[hub_node_id])
            {
                TryInsertLabel(Label[u], hub_node_id, current_bn, path_predecessor_nodes[u]);
            }

            for (const auto &[v, weight] : Adj[u])
            {
                int new_bn = std::min(current_bn, weight);
                if (new_bn > best_bottleneck[v])
                {
                    best_bottleneck[v] = new_bn;
                    path_predecessor_nodes[v] = u;
                    pq.push({new_bn, v});
                }
            }
        }
    }

    void BottleneckGraph::ConstructIndex()
    {

        if (LoadIndexFromFile())
        {
            std::cout << "成功从缓存文件加载索引: " << index_file_path << std::endl;
            return;
        }
        std::cout << "未找到索引缓存文件或加载失败，开始构建新索引..." << std::endl;

        for (int i = 1; i <= n; ++i)
        {
            if (i < Label.size())
                TryInsertLabel(Label[i], i, INF_WEIGHT, INVALID_PREDECESSOR);
        }

        std::vector<DegreeNode> sorted_degree_nodes;
        for (int i = 1; i <= n; ++i)
            sorted_degree_nodes.push_back(degreeList[i]);
        std::sort(sorted_degree_nodes.begin(), sorted_degree_nodes.end(), CmpDegreeNode);

        for (const auto &ranked_node_info : sorted_degree_nodes)
        {
            int hub_candidate_id = ranked_node_info.id;
            if (hub_candidate_id > n || hub_candidate_id <= 0)
                continue;

            BFSForIndexOptimized(hub_candidate_id);
        }

        if (SaveIndexToFile())
        {
            std::cout << "索引已保存到缓存文件: " << index_file_path << std::endl;
        }
        else
        {
            std::cerr << "警告: 索引保存失败" << std::endl;
        }
    }

    bool BottleneckGraph::SaveIndexToFile() const
    {
        try
        {
            std::ofstream out(index_file_path, std::ios::binary);
            if (!out.is_open())
            {
                std::cerr << "无法创建索引文件: " << index_file_path << std::endl;
                return false;
            }

            const char *magic = "BTLNCK_IDX";
            out.write(magic, 10);

            out.write(reinterpret_cast<const char *>(&n), sizeof(n));
            out.write(reinterpret_cast<const char *>(&m), sizeof(m));

            for (int i = 0; i <= n; ++i)
            {
                out.write(reinterpret_cast<const char *>(&rankList[i]), sizeof(int));
            }

            for (int i = 0; i <= n; ++i)
            {
                out.write(reinterpret_cast<const char *>(&degreeList[i].id), sizeof(int));
                out.write(reinterpret_cast<const char *>(&degreeList[i].num), sizeof(int));
            }

            for (int i = 1; i <= n; ++i)
            {
                if (i >= Label.size())
                {
                    int label_size = 0;
                    out.write(reinterpret_cast<const char *>(&label_size), sizeof(int));
                    continue;
                }

                int label_size = Label[i].size();
                out.write(reinterpret_cast<const char *>(&label_size), sizeof(int));

                for (const auto &label_pair : Label[i])
                {
                    const BottleneckLabelEntry &entry = label_pair.second;
                    out.write(reinterpret_cast<const char *>(&entry.hub_id), sizeof(int));
                    out.write(reinterpret_cast<const char *>(&entry.bottleneck_weight), sizeof(int));
                    out.write(reinterpret_cast<const char *>(&entry.predecessor_on_path_to_hub), sizeof(int));
                }
            }

            out.close();
            return true;
        }
        catch (const std::exception &e)
        {
            std::cerr << "保存索引时发生异常: " << e.what() << std::endl;
            return false;
        }
    }

    bool BottleneckGraph::LoadIndexFromFile()
    {
        try
        {
            std::ifstream in(index_file_path, std::ios::binary);
            if (!in.is_open())
            {
                return false;
            }

            char magic[11] = {0};
            in.read(magic, 10);
            if (std::string(magic) != "BTLNCK_IDX")
            {
                std::cerr << "索引文件格式无效: " << index_file_path << std::endl;
                in.close();
                return false;
            }

            int file_n, file_m;
            in.read(reinterpret_cast<char *>(&file_n), sizeof(int));
            in.read(reinterpret_cast<char *>(&file_m), sizeof(long long));

            if (file_n != n || file_m != m)
            {
                std::cerr << "索引文件与当前图不匹配 - 文件: n=" << file_n << ", m=" << file_m
                          << " 当前: n=" << n << ", m=" << m << std::endl;
                in.close();
                return false;
            }

            rankList.resize(n + 1);
            for (int i = 0; i <= n; ++i)
            {
                in.read(reinterpret_cast<char *>(&rankList[i]), sizeof(int));
            }

            degreeList.resize(n + 1);
            for (int i = 0; i <= n; ++i)
            {
                in.read(reinterpret_cast<char *>(&degreeList[i].id), sizeof(int));
                in.read(reinterpret_cast<char *>(&degreeList[i].num), sizeof(int));
            }

            Label.clear();
            Label.resize(n + 1);

            for (int i = 1; i <= n; ++i)
            {
                int label_size;
                in.read(reinterpret_cast<char *>(&label_size), sizeof(int));

                for (int j = 0; j < label_size; ++j)
                {
                    BottleneckLabelEntry entry;
                    in.read(reinterpret_cast<char *>(&entry.hub_id), sizeof(int));
                    in.read(reinterpret_cast<char *>(&entry.bottleneck_weight), sizeof(int));
                    in.read(reinterpret_cast<char *>(&entry.predecessor_on_path_to_hub), sizeof(int));

                    Label[i][entry.hub_id] = entry;
                }
            }

            in.close();

            if (in.bad())
            {
                std::cerr << "读取索引文件时发生IO错误" << std::endl;
                return false;
            }

            return true;
        }
        catch (const std::exception &e)
        {
            std::cerr << "加载索引时发生异常: " << e.what() << std::endl;
            return false;
        }
    }

    void BottleneckGraph::GenerateIndexFilePath(const std::string &graph_file_path)
    {

        index_file_path = graph_file_path + ".bottleneck_index";
    }

    bool BottleneckGraph::Query(int s, int t, int query_weight)
    {
        if (s > n || t > n || s <= 0 || t <= 0)
            return false;
        if (s == t)
            return INF_WEIGHT >= query_weight;

        for (const auto &s_label_pair : Label[s])
        {
            int hub_id = s_label_pair.first;
            const BottleneckLabelEntry &s_entry = s_label_pair.second;

            auto t_label_it = Label[t].find(hub_id);
            if (t_label_it != Label[t].end())
            {
                const BottleneckLabelEntry &t_entry = t_label_it->second;

                int path_bottleneck = std::min(s_entry.bottleneck_weight, t_entry.bottleneck_weight);
                if (path_bottleneck >= query_weight)
                {
                    return true;
                }
            }
        }
        return false;
    }

    bool BottleneckGraph::QueryBFS(int s, int t, int query_weight)
    {
        if (s > n || t > n || s <= 0 || t <= 0)
            return false;
        if (s == t)
            return true;

        std::queue<int> q_bfs;
        std::vector<bool> visited(n + 1, false);

        q_bfs.push(s);
        visited[s] = true;

        while (!q_bfs.empty())
        {
            int u = q_bfs.front();
            q_bfs.pop();

            for (const auto &edge_pair : Adj[u])
            {
                int v = edge_pair.first;
                int weight = edge_pair.second;

                if (weight < query_weight)
                    continue;

                if (v == t)
                    return true;

                if (!visited[v])
                {
                    visited[v] = true;
                    q_bfs.push(v);
                }
            }
        }

        return false;
    }

    void BottleneckGraph::RunIncrementalBFSUpdate(int source_node, int first_target_node, int bn_to_first_target, int pred_of_first_target_from_source)
    {
        std::queue<std::tuple<int, int, int>> q;
        std::vector<int> max_bn_found_in_this_bfs(n + 1, 0);

        if (source_node <= 0 || source_node > n || first_target_node <= 0 || first_target_node > n)
            return;

        if (rankList[source_node] > rankList[first_target_node])
        {
            TryInsertLabel(Label[source_node], first_target_node, bn_to_first_target, pred_of_first_target_from_source);
        }

        q.push({first_target_node, bn_to_first_target, pred_of_first_target_from_source});
        max_bn_found_in_this_bfs[first_target_node] = bn_to_first_target;

        while (!q.empty())
        {
            auto [curr_node, bn_source_to_curr, pred_of_curr_node_on_path_from_source] = q.front();
            q.pop();

            if (curr_node <= 0 || curr_node > n)
                continue;

            for (const auto &edge_pair : Adj[curr_node])
            {
                int neighbor = edge_pair.first;
                int edge_weight = edge_pair.second;

                if (neighbor <= 0 || neighbor > n)
                    continue;

                if (curr_node == first_target_node && neighbor == source_node && pred_of_curr_node_on_path_from_source == source_node)
                {
                    continue;
                }

                int bn_source_to_neighbor = std::min(bn_source_to_curr, edge_weight);

                if (rankList[source_node] > rankList[neighbor])
                {
                    TryInsertLabel(Label[source_node], neighbor, bn_source_to_neighbor, curr_node);
                }

                if (bn_source_to_neighbor > max_bn_found_in_this_bfs[neighbor])
                {
                    max_bn_found_in_this_bfs[neighbor] = bn_source_to_neighbor;
                    q.push({neighbor, bn_source_to_neighbor, curr_node});
                }
            }
        }
    }

    void BottleneckGraph::DynamicDeleteEdge(int u, int v, int weight)
    {
        if (u < 0 || u > n || v < 0 || v > n)
        {
#ifdef DEBUG
            std::cerr << "Error in DynamicDeleteEdge: Vertex " << u << " or " << v
                      << " is out of bounds (n=" << n << "). Cannot delete edge." << std::endl;
#endif
            return;
        }
        if (u == v)
        {
            std::cerr << "Error in DynamicDeleteEdge: Self-loops (" << u << "," << v << ") are not processed." << std::endl;
            return;
        }

        int old_weight = 0;
        bool edge_exists = false;
        for (const auto &edge_pair : Adj[u])
        {
            if (edge_pair.first == v)
            {
                edge_exists = true;
                old_weight = edge_pair.second;
                break;
            }
        }

        if (!edge_exists)
        {
#ifdef DEBUG
            std::cerr << "Warning: Edge (" << u << "," << v << ") does not exist. Nothing to delete." << std::endl;
#endif
            return;
        }

        int new_weight = old_weight - weight;
        bool complete_removal = (new_weight <= 0);

        std::set<int> affected_candidates;
        std::queue<std::pair<int, int>> bfs_q;
        std::vector<bool> visited_bfs(n + 1, false);
        const int MAX_DEPTH = 15;

        auto rebuild_start = std::chrono::high_resolution_clock::now();

        bfs_q.push({u, 0});
        visited_bfs[u] = true;
        affected_candidates.insert(u);

        while (!bfs_q.empty())
        {
            auto [curr, depth] = bfs_q.front();
            bfs_q.pop();

            if (depth >= MAX_DEPTH)
                continue;

            for (const auto &edge_pair : Adj[curr])
            {
                int neighbor = edge_pair.first;
                if (!visited_bfs[neighbor])
                {
                    visited_bfs[neighbor] = true;
                    affected_candidates.insert(neighbor);
                    bfs_q.push({neighbor, depth + 1});
                }
            }
        }

        std::fill(visited_bfs.begin(), visited_bfs.end(), false);
        bfs_q.push({v, 0});
        visited_bfs[v] = true;
        affected_candidates.insert(v);

        while (!bfs_q.empty())
        {
            auto [curr, depth] = bfs_q.front();
            bfs_q.pop();

            if (depth >= MAX_DEPTH)
                continue;

            for (const auto &edge_pair : Adj[curr])
            {
                int neighbor = edge_pair.first;
                if (!visited_bfs[neighbor])
                {
                    visited_bfs[neighbor] = true;
                    affected_candidates.insert(neighbor);
                    bfs_q.push({neighbor, depth + 1});
                }
            }
        }

        if (complete_removal)
        {

#ifdef DEBUG
            printf("Completely removing edge (%d, %d) with weight %d\n", u, v, old_weight);
#endif

            for (auto it = Adj[u].begin(); it != Adj[u].end();)
            {
                if (it->first == v)
                {
                    it = Adj[u].erase(it);
                }
                else
                {
                    ++it;
                }
            }

            for (auto it = Adj[v].begin(); it != Adj[v].end();)
            {
                if (it->first == u)
                {
                    it = Adj[v].erase(it);
                }
                else
                {
                    ++it;
                }
            }

            m--;
        }
        else
        {

#ifdef DEBUG
            printf("Reducing edge (%d, %d) weight from %d to %d\n", u, v, old_weight, new_weight);
#endif
            for (auto &edge_pair : Adj[u])
            {
                if (edge_pair.first == v)
                {
                    edge_pair.second = new_weight;
                    break;
                }
            }

            for (auto &edge_pair : Adj[v])
            {
                if (edge_pair.first == u)
                {
                    edge_pair.second = new_weight;
                    break;
                }
            }
        }

        for (int node_id : affected_candidates)
        {
            if (node_id < 1 || node_id >= Label.size())
                continue;

            Label[node_id].clear();

            if (node_id <= n)
            {
                TryInsertLabel(Label[node_id], node_id, INF_WEIGHT, INVALID_PREDECESSOR);
            }
        }

        std::vector<DegreeNode> sorted_affected_nodes;
        for (int node_id : affected_candidates)
        {
            if (node_id <= 0 || node_id > n)
                continue;
            sorted_affected_nodes.push_back(degreeList[node_id]);
        }

        std::sort(sorted_affected_nodes.begin(), sorted_affected_nodes.end(), CmpDegreeNode);

        for (const auto &node_info : sorted_affected_nodes)
        {
            BFSForIndex(node_info.id);
        }

        auto rebuild_end = std::chrono::high_resolution_clock::now();
        rebuild_time += std::chrono::duration_cast<std::chrono::nanoseconds>(rebuild_end - rebuild_start).count();
    }

    void BottleneckGraph::DynamicBatchAdd(const std::vector<std::tuple<int, int, int>> &addedEdges)
    {
        std::cout << "DynamicBatchAdd called with " << addedEdges.size() << " edges - Not Implemented" << std::endl;
    }

    void BottleneckGraph::DynamicBatchDelete(const std::vector<std::tuple<int, int, int>> &deletedEdges)
    {
        std::cout << "DynamicBatchDelete called with " << deletedEdges.size() << " edges - Not Implemented" << std::endl;
    }

    unsigned long long BottleneckGraph::GetIndexSize()
    {
        unsigned long long count = 0;
        for (int i = 1; i <= n; ++i)
        {
            if (i < Label.size())
                count += Label[i].size();
        }
        return count * sizeof(BottleneckLabelEntry);
    }

    void BottleneckGraph::PrintStat()
    {
        std::cout << "BottleneckGraph Stats (Undirected):" << std::endl;
        std::cout << "Nodes (n): " << n << std::endl;
        std::cout << "Edges (m): " << m << std::endl;

        std::ifstream check(index_file_path);
        if (check.good())
        {
            check.seekg(0, std::ios::end);
            long long file_size = check.tellg();
            std::cout << "Index File Size: " << file_size << " bytes ("
                      << file_size / (1024.0 * 1024.0) << " MB)" << std::endl;
            check.close();
        }
        else
        {
            std::cout << "Index File: Not found" << std::endl;
        }

        std::cout << "Index Size (approx bytes): " << GetIndexSize() << std::endl;
        unsigned long long total_labels = 0;
        for (int i = 1; i <= n; ++i)
        {
            if (i < Label.size())
                total_labels += Label[i].size();
        }
        std::cout << "Total label entries: " << total_labels << std::endl;
    }

    void BottleneckGraph::PrintNodeLabel(int node_id)
    {
        if (node_id <= 0 || node_id > n)
        {
            printf("Node ID %d is out of bounds.\n", node_id);
            return;
        }
        printf("Labels for Node %d (Rank: %d):\n", node_id, rankList[node_id]);
        if (Label[node_id].empty())
        {
            printf("  No labels.\n");
            return;
        }
        for (const auto &pair : Label[node_id])
        {
            const BottleneckLabelEntry &entry = pair.second;
            printf("  Hub: %d (Rank: %d), Bottleneck: %d", entry.hub_id, rankList[entry.hub_id], entry.bottleneck_weight);
            if (entry.predecessor_on_path_to_hub != INVALID_PREDECESSOR)
            {
                int pred_node = entry.predecessor_on_path_to_hub;
                int edge_w_to_hub = -1;
                bool found_edge = false;
                if (entry.hub_id > 0 && entry.hub_id <= n)
                {
                    for (const auto &edge_pair : Adj[entry.hub_id])
                    {
                        if (edge_pair.first == pred_node)
                        {
                            edge_w_to_hub = edge_pair.second;
                            found_edge = true;
                            break;
                        }
                    }
                }
                if (found_edge)
                {
                    printf(", PathToHubVia: %d (Edge: (%d,%d) w:%d)\n", pred_node, pred_node, entry.hub_id, edge_w_to_hub);
                }
                else
                {
                    printf(", PathToHubVia: %d (Edge weight not found in Adj[%d])\n", pred_node, entry.hub_id);
                }
            }
            else
            {
                printf(", PathToHubVia: null (self-loop or initial)\n");
            }
        }
        printf("--- End Labels for Node %d ---\n", node_id);
    }

    void BottleneckGraph::SpreadHubLabelUpdate(int start_node, int hub_id, int current_bottleneck, int predecessor)
    {
        std::queue<std::tuple<int, int, int>> q;
        std::vector<int> max_bn_found(n + 1, 0);

        q.push({start_node, current_bottleneck, predecessor});
        max_bn_found[start_node] = current_bottleneck;

        while (!q.empty())
        {
            auto [curr, bn_hub_to_curr, pred] = q.front();
            q.pop();

            if (curr <= 0 || curr > n)
                continue;

            for (const auto &edge_pair : Adj[curr])
            {
                int neighbor = edge_pair.first;
                int edge_weight = edge_pair.second;

                if (neighbor <= 0 || neighbor > n)
                    continue;

                int bn_hub_to_neighbor = std::min(bn_hub_to_curr, edge_weight);

                if (bn_hub_to_neighbor > max_bn_found[neighbor] &&
                    rankList[neighbor] > rankList[hub_id])
                {
                    max_bn_found[neighbor] = bn_hub_to_neighbor;
                    if (TryInsertLabel(Label[neighbor], hub_id, bn_hub_to_neighbor, curr))
                    {
                        q.push({neighbor, bn_hub_to_neighbor, curr});
                    }
                }
            }
        }
    }

    void BottleneckGraph::RunComprehensiveLabelUpdate(int u, int v, int current_edge_weight)
    {

        std::set<int> affected_candidates;
        std::queue<std::pair<int, int>> bfs_q;
        std::vector<bool> visited_bfs(n + 1, false);
        const int MAX_DEPTH = 10;

        bfs_q.push({u, 0});
        visited_bfs[u] = true;
        affected_candidates.insert(u);

        while (!bfs_q.empty())
        {
            auto [curr, depth] = bfs_q.front();
            bfs_q.pop();

            if (depth >= MAX_DEPTH)
                continue;

            for (const auto &edge_pair : Adj[curr])
            {
                int neighbor = edge_pair.first;
                if (!visited_bfs[neighbor])
                {
                    visited_bfs[neighbor] = true;
                    affected_candidates.insert(neighbor);
                    bfs_q.push({neighbor, depth + 1});
                }
            }
        }

        std::fill(visited_bfs.begin(), visited_bfs.end(), false);
        bfs_q.push({v, 0});
        visited_bfs[v] = true;
        affected_candidates.insert(v);

        while (!bfs_q.empty())
        {
            auto [curr, depth] = bfs_q.front();
            bfs_q.pop();

            if (depth >= MAX_DEPTH)
                continue;

            for (const auto &edge_pair : Adj[curr])
            {
                int neighbor = edge_pair.first;
                if (!visited_bfs[neighbor])
                {
                    visited_bfs[neighbor] = true;
                    affected_candidates.insert(neighbor);
                    bfs_q.push({neighbor, depth + 1});
                }
            }
        }

        for (int cand_node : affected_candidates)
        {
            if (cand_node <= 0 || cand_node > n)
                continue;
            BFSForIndex(cand_node);
        }
    }

    void BottleneckGraph::DynamicAddEdge(int u, int v, int weight)
    {
        if (u <= 0 || u > n || v <= 0 || v > n)
        {
            std::cerr << "Error in DynamicAddEdge: Vertex " << u << " or " << v
                      << " is out of bounds (n=" << n << "). Edge not added." << std::endl;
            return;
        }
        if (u == v)
        {
            std::cerr << "Error in DynamicAddEdge: Self-loops (" << u << "," << v << ") are not processed for dynamic addition." << std::endl;
            return;
        }

        int old_weight = 0;
        bool edge_exists = false;

        for (const auto &edge_pair : Adj[u])
        {
            if (edge_pair.first == v)
            {
                edge_exists = true;
                old_weight = edge_pair.second;
                break;
            }
        }

        int current_edge_weight;

        if (edge_exists)
        {
            current_edge_weight = old_weight + weight;
#ifdef DEBUG
            printf("Attempting to update edge (%d, %d) from old weight %d to new weight %d.\n", u, v, old_weight, current_edge_weight);
#endif
            bool updated_u = false, updated_v = false;
            for (auto &edge_pair : Adj[u])
            {
                if (edge_pair.first == v)
                {
                    edge_pair.second = current_edge_weight;
                    updated_u = true;
                    break;
                }
            }
            for (auto &edge_pair : Adj[v])
            {
                if (edge_pair.first == u)
                {
                    edge_pair.second = current_edge_weight;
                    updated_v = true;
                    break;
                }
            }
            if (!updated_u || !updated_v)
            {
                std::cerr << "Error: Inconsistency found while trying to update edge weights for (" << u << "," << v << "). Aborting update." << std::endl;
                return;
            }
        }
        else
        {
            current_edge_weight = weight;
            Adj[u].push_back({v, current_edge_weight});
            Adj[v].push_back({u, current_edge_weight});
            m++;
#ifdef DEBUG
            printf("Dynamically added new edge (%d, %d, %d)\n", u, v, current_edge_weight);
#endif
        }

        /*
        if (rankList[u] > rankList[v]) {
            TryInsertLabel(Label[u], v, current_edge_weight, u);
        }

        RunIncrementalBFSUpdate(u, v, current_edge_weight, u);
        RunIncrementalBFSUpdate(v, u, current_edge_weight, v);
        */

        RunComprehensiveLabelUpdate(u, v, current_edge_weight);
    }

}
