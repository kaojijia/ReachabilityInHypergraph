
#ifndef WEIGHTED_GRAPH_H
#define WEIGHTED_GRAPH_H

#include <vector>
#include <limits>
#include <stdexcept>
#include <memory>
#include <queue>
#include <utility>
#include <numeric>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iostream>

using namespace std;

class GraphDisjointSets
{
public:
    GraphDisjointSets(size_t n) : parent(n), rank(n, 0)
    {
        for (size_t i = 0; i < n; i++)
            parent[i] = i;
    }

    int find(int x)
    {
        if (parent[x] != x)
            parent[x] = find(parent[x]);
        return parent[x];
    }

    void merge(int x, int y)
    {
        int root_x = find(x);
        int root_y = find(y);

        if (root_x == root_y)
            return;

        if (rank[root_x] < rank[root_y])
            parent[root_x] = root_y;
        else if (rank[root_x] > rank[root_y])
            parent[root_y] = root_x;
        else
        {
            parent[root_y] = root_x;
            rank[root_x]++;
        }
    }

    std::vector<int> parent;
    std::vector<int> rank;

private:
};

class WeightedGraph
{
public:
    WeightedGraph(size_t num_vertices = 0, size_t min_weight = 1, size_t estimated_edges_per_vertex = 10)
        : adj_list(num_vertices), min_weight(min_weight)
    {

        for (auto &neighbors : adj_list)
        {
            neighbors.reserve(estimated_edges_per_vertex);
        }

        labels.resize(num_vertices);
    }

    void addVertices(size_t count, size_t estimated_edges_per_vertex = 10)
    {
        size_t original_size = adj_list.size();
        adj_list.resize(original_size + count);

        for (size_t i = original_size; i < adj_list.size(); i++)
        {
            adj_list[i].reserve(estimated_edges_per_vertex);
        }
    }

    int addVertex(size_t estimated_edges = 10)
    {
        adj_list.push_back(std::vector<std::pair<int, int>>());
        adj_list.back().reserve(estimated_edges);
        return adj_list.size() - 1;
    }

    void addEdges(const std::vector<std::tuple<int, int, int>> &edges)
    {
        for (const auto &[u, v, weight] : edges)
        {
            if (u < 0 || u >= static_cast<int>(adj_list.size()) ||
                v < 0 || v >= static_cast<int>(adj_list.size()))
                throw std::invalid_argument("Vertex does not exist");

            adj_list[u].emplace_back(v, weight);
            adj_list[v].emplace_back(u, weight);
        }
    }

    void addEdge(int u, int v, int weight)
    {
        if (u < 0 || u >= static_cast<int>(adj_list.size()) ||
            v < 0 || v >= static_cast<int>(adj_list.size()))
            throw std::invalid_argument("Vertex does not exist");

        adj_list[u].emplace_back(v, weight);
        adj_list[v].emplace_back(u, weight);
    }

    size_t numVertices() const
    {
        return adj_list.size();
    }

    size_t numEdges() const
    {
        size_t edge_count = 0;
        for (size_t u = 0; u < adj_list.size(); ++u)
        {
            for (const auto &edge : adj_list[u])
            {
                if (u < static_cast<size_t>(edge.first))
                {
                    edge_count++;
                }
            }
        }
        return edge_count;
    }

    const std::vector<std::pair<int, int>> &getNeighbors(int vertex) const
    {
        if (vertex < 0 || vertex >= static_cast<int>(adj_list.size()))
            throw std::invalid_argument("Vertex does not exist");
        return adj_list[vertex];
    }

    void offline_industry()
    {

        ds = std::make_unique<GraphDisjointSets>(adj_list.size());

        for (size_t u = 0; u < adj_list.size(); u++)
        {
            for (const auto &[v, weight] : adj_list[u])
            {

                if (u < static_cast<size_t>(v) && static_cast<size_t>(weight) >= min_weight)
                {
                    ds->merge(u, v);
                }
            }
        }
    }

    bool disjointSet_reachability_query(int source, int target) const
    {
        if (source < 0 || source >= static_cast<int>(adj_list.size()) ||
            target < 0 || target >= static_cast<int>(adj_list.size()))
            throw std::invalid_argument("Vertex does not exist");

        if (source == target)
            return true;

        if (!ds)
        {

            throw std::runtime_error("Disjoint set not built. Call offline_industry first.");
        }

        return ds->find(source) == ds->find(target);
    }

    bool landmark_reachability_query(int u, int v) const
    {
        if (u < 0 || u >= (int)numVertices() ||
            v < 0 || v >= (int)numVertices())
            throw std::invalid_argument("Vertex does not exist");

        if (u == v)
            return true;

        return intersect(labels[u], labels[v]);
    }

    double getAdjListMemoryUsageMB() const
    {
        size_t memory_bytes = 0;
        memory_bytes += sizeof(adj_list);
        for (const auto &neighbors : adj_list)
        {
            memory_bytes += sizeof(neighbors);

            memory_bytes += neighbors.capacity() * sizeof(std::pair<int, int>);
        }
        return static_cast<double>(memory_bytes) / (1024.0 * 1024.0);
    }

    double getDsMemoryUsageMB() const
    {
        size_t memory_bytes = 0;
        if (ds)
        {
            memory_bytes += sizeof(*ds);
            memory_bytes += ds->parent.capacity() * sizeof(int);
            memory_bytes += ds->rank.capacity() * sizeof(int);
        }
        return static_cast<double>(memory_bytes) / (1024.0 * 1024.0);
    }

    double getMemoryUsageMB() const
    {
        return getAdjListMemoryUsageMB() + getDsMemoryUsageMB();
    }

    friend class Hypergraph;

    friend class WeightedPrunedLandmarkIndex;

    bool saveAdjList(const std::string &filename) const
    {
        std::ofstream file(filename);
        if (!file.is_open())
        {
            std::cerr << "Error: Cannot open file for writing WeightedGraph AdjList: " << filename << std::endl;
            return false;
        }

        file << adj_list.size() << " " << numEdges() << "\n";

        for (size_t u = 0; u < adj_list.size(); ++u)
        {

            for (const auto &edge_pair : adj_list[u])
            {
                int v = edge_pair.first;
                int weight = edge_pair.second;

                if (u < static_cast<size_t>(v))
                {
                    file << u << " " << v << " " << weight << "\n";
                }
            }
        }
        return file.good();
    }

    bool loadAdjList(const std::string &filename)
    {
        std::ifstream file(filename);
        if (!file.is_open())
        {

            return false;
        }
        size_t num_vertices_expected;
        file >> num_vertices_expected;
        if (file.fail() || num_vertices_expected != adj_list.size())
        {
            std::cerr << "Error: WeightedGraph AdjList cache size mismatch or read error." << std::endl;
            return false;
        }

        for (auto &neighbors : adj_list)
        {
            neighbors.clear();
        }

        std::string line;
        std::getline(file, line);
        int u, neighbor, weight;
        while (std::getline(file, line))
        {
            std::istringstream iss(line);
            iss >> u;
            if (iss.fail() || u < 0 || u >= num_vertices_expected)
            {
                std::cerr << "Warning: Skipping invalid line in WeightedGraph AdjList cache: " << line << std::endl;
                continue;
            }
            while (iss >> neighbor >> weight)
            {
                if (neighbor >= 0 && neighbor < num_vertices_expected)
                {
                    adj_list[u].emplace_back(neighbor, weight);
                }
                else
                {
                    std::cerr << "Warning: Invalid neighbor ID " << neighbor << " for vertex " << u << " in AdjList cache." << std::endl;
                }
            }
        }
        if (file.bad())
        {
            std::cerr << "Error: File read error occurred during WeightedGraph AdjList cache loading." << std::endl;
            return false;
        }
        return true;
    }

    bool saveDisjointSets(const std::string &filename) const
    {
        if (!ds)
        {
            std::cerr << "Error: Cannot save WeightedGraph DisjointSets, it's not built." << std::endl;
            return false;
        }
        std::ofstream file(filename);
        if (!file.is_open())
        {
            std::cerr << "Error: Cannot open file for writing WeightedGraph DisjointSets: " << filename << std::endl;
            return false;
        }
        file << ds->parent.size() << "\n";
        for (size_t i = 0; i < ds->parent.size(); ++i)
        {
            file << ds->parent[i] << " " << ds->rank[i] << "\n";
        }
        return file.good();
    }

    bool loadDisjointSets(const std::string &filename)
    {
        std::ifstream file(filename);
        if (!file.is_open())
        {

            return false;
        }
        size_t n_expected;
        file >> n_expected;
        if (file.fail() || n_expected != adj_list.size())
        {
            std::cerr << "Error: WeightedGraph DisjointSets cache size mismatch or read error." << std::endl;
            return false;
        }

        ds = std::make_unique<GraphDisjointSets>(n_expected);
        for (size_t i = 0; i < n_expected; ++i)
        {
            if (!(file >> ds->parent[i] >> ds->rank[i]))
            {
                std::cerr << "Error: Failed to read parent/rank data from WeightedGraph DisjointSets cache for index " << i << "." << std::endl;
                ds.reset();
                return false;
            }
        }
        if (file.bad())
        {
            std::cerr << "Error: File read error occurred during WeightedGraph DisjointSets cache loading." << std::endl;
            return false;
        }

        int check_eof;
        file >> check_eof;
        if (!file.eof())
        {
            std::cerr << "Warning: Extra data found at the end of WeightedGraph DisjointSets cache file: " << filename << std::endl;
        }

        return true;
    }

private:
    std::vector<int> getConnectedComponent(int vertex) const
    {
        if (vertex < 0 || vertex >= static_cast<int>(adj_list.size()))
            throw std::invalid_argument("Vertex does not exist");

        int root = ds->find(vertex);
        std::vector<int> component;

        for (size_t i = 0; i < adj_list.size(); i++)
        {
            if (ds->find(i) == root)
            {
                component.push_back(i);
            }
        }

        return component;
    }

    std::vector<std::vector<int>> getAllConnectedComponents() const
    {

        std::vector<int> root_to_index(adj_list.size(), -1);
        std::vector<std::vector<int>> components;

        std::vector<int> component_size(adj_list.size(), 0);
        for (size_t i = 0; i < adj_list.size(); i++)
        {
            int root = ds->find(i);
            component_size[root]++;
        }

        int component_count = 0;
        for (size_t i = 0; i < adj_list.size(); i++)
        {
            if (component_size[i] > 0)
            {
                root_to_index[i] = component_count++;
                components.push_back(std::vector<int>());
                components.back().reserve(component_size[i]);
            }
        }

        for (size_t i = 0; i < adj_list.size(); i++)
        {
            int root = ds->find(i);
            components[root_to_index[root]].push_back(i);
        }

        return components;
    }

    void buildPrunedLandmarkIndex()
    {

        int n = numVertices();

        for (auto &label : labels)
        {
            label.clear();
        }

        std::vector<std::pair<int, int>> degrees(n);
        for (int i = 0; i < n; ++i)
        {
            degrees[i] = {static_cast<int>(adj_list[i].size()), i};
        }

        std::sort(degrees.begin(), degrees.end(), [](const std::pair<int, int> &a, const std::pair<int, int> &b)
                  { return a.first > b.first; });

        std::vector<int> verticesOrder(n);
        for (int i = 0; i < n; ++i)
        {
            verticesOrder[i] = degrees[i].second;
        }

        for (int landmark : verticesOrder)
        {
            prunedBFS(landmark);
        }

        for (int i = 0; i < n; ++i)
        {
            insertSorted(labels[i], i);
        }
    }

    void prunedBFS(int curLandmark)
    {

        std::vector<bool> visited(numVertices(), false);
        visited[curLandmark] = true;

        std::queue<int> q;
        q.push(curLandmark);

        while (!q.empty())
        {
            int current = q.front();
            q.pop();

            if (hopQuery(curLandmark, current) && current != curLandmark)
            {
                continue;
            }

            if (current != curLandmark)
            {
                insertSorted(labels[current], curLandmark);
            }

            const auto &neighbors = adj_list[current];
            for (auto &[nbr, w] : neighbors)
            {
                if (w >= this->min_weight && !visited[nbr])
                {
                    visited[nbr] = true;
                    q.push(nbr);
                }
            }
        }
    }

    bool hopQuery(int curLandmark, int node) const
    {

        if (node == curLandmark)
        {
            return true;
        }

        return intersect(labels[curLandmark], labels[node]);
    }

    void insertSorted(std::vector<int> &vec, int val) const
    {
        auto it = std::lower_bound(vec.begin(), vec.end(), val);
        if (it == vec.end() || *it != val)
        {
            vec.insert(it, val);
        }
    }

    bool intersect(const std::vector<int> &vec1, const std::vector<int> &vec2) const
    {
        int i = 0, j = 0;
        while (i < (int)vec1.size() && j < (int)vec2.size())
        {
            if (vec1[i] == vec2[j])
            {
                return true;
            }
            else if (vec1[i] < vec2[j])
            {
                i++;
            }
            else
            {
                j++;
            }
        }
        return false;
    }

    std::vector<std::vector<std::pair<int, int>>> adj_list;

    mutable std::unique_ptr<GraphDisjointSets> ds;

    std::vector<std::vector<int>> labels;

    size_t min_weight;
};

#endif