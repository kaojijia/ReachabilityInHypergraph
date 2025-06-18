#ifndef HYPERGRAPH_H
#define HYPERGRAPH_H
#include <thread>
#include <mutex>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <string>
#include <fstream>
#include <sstream>
#include <unordered_set>
#include <queue>
#include <memory>
#include <numeric>
#include "WeightedGraph.h"
#include "UWeightedPLL.h"

using namespace std;

class DisjointSets
{
public:
    DisjointSets(size_t n) : parent(n), rank(n, 0)
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

    bool saveToFile(const std::string &filename) const
    {
        std::ofstream file(filename);
        if (!file.is_open())
        {
            std::cerr << "Error: Cannot open file for writing DisjointSets: " << filename << std::endl;
            return false;
        }
        file << parent.size() << "\n";
        for (size_t i = 0; i < parent.size(); ++i)
        {
            file << parent[i] << " " << rank[i] << "\n";
        }
        return file.good();
    }

    bool loadFromFile(const std::string &filename)
    {
        std::ifstream file(filename);
        if (!file.is_open())
        {

            return false;
        }
        size_t n_expected;
        file >> n_expected;
        if (file.fail())
        {
            std::cerr << "Error: Failed to read size from DisjointSets cache: " << filename << std::endl;
            return false;
        }

        parent.resize(n_expected);
        rank.resize(n_expected);

        for (size_t i = 0; i < n_expected; ++i)
        {
            if (!(file >> parent[i] >> rank[i]))
            {
                std::cerr << "Error: Failed to read parent/rank data from DisjointSets cache for index " << i << " in file: " << filename << std::endl;

                parent.clear();
                rank.clear();
                return false;
            }
        }
        if (file.bad())
        {
            std::cerr << "Error: File read error occurred during DisjointSets cache loading." << std::endl;
            return false;
        }

        int check_eof;
        file >> check_eof;
        if (!file.eof())
        {
            std::cerr << "Warning: Extra data found at the end of DisjointSets cache file: " << filename << std::endl;
        }

        return true;
    }

    std::vector<int> parent;
    std::vector<int> rank;

private:
};

struct Hyperedge
{
    std::vector<int> vertices;

    Hyperedge()
    {

        vertices.reserve(16);
    }

    Hyperedge(size_t expected_size)
    {
        vertices.reserve(expected_size);
    }

    Hyperedge(const std::vector<int> &vertices_)
        : vertices(vertices_) {}

    size_t size() const
    {
        return vertices.size();
    }

    std::vector<int> intersection(const Hyperedge &other) const
    {
        std::vector<int> result;
        result.reserve(std::min(vertices.size(), other.vertices.size()));

        std::vector<int> sorted_vertices = vertices;
        std::vector<int> other_sorted = other.vertices;
        std::sort(sorted_vertices.begin(), sorted_vertices.end());
        std::sort(other_sorted.begin(), other_sorted.end());

        std::set_intersection(
            sorted_vertices.begin(), sorted_vertices.end(),
            other_sorted.begin(), other_sorted.end(),
            std::back_inserter(result));
        return result;
    }
};

class Hypergraph
{
public:
    Hypergraph(size_t expected_vertices = 100, size_t expected_edges = 100)
        : ds_valid(false), pll_graph(nullptr), pll(nullptr)
    {
        vertices.reserve(expected_vertices);
        hyperedges.reserve(expected_edges);

        vertex_to_edges.reserve(expected_vertices);
    }

    std::vector<int> getHyperedgeIntersection(int edgeId1, int edgeId2) const
    {
        return getHyperedge(edgeId1).intersection(getHyperedge(edgeId2));
    }

    DisjointSets buildDisjointSets() const
    {
        DisjointSets ds(vertices.size());

        for (const auto &edge : hyperedges)
        {
            if (edge.vertices.empty())
                continue;

            int first = edge.vertices[0];
            for (size_t i = 1; i < edge.vertices.size(); i++)
            {
                ds.merge(first, edge.vertices[i]);
            }
        }

        return ds;
    }

    void addVertexToEdges(int vertexId, const std::vector<int> &edgeList)
    {

        if (vertexId < 0)
        {
            throw std::invalid_argument("Vertex ID cannot be negative");
        }

        if (vertexId >= static_cast<int>(vertex_to_edges.size()) || vertex_to_edges[vertexId].empty())
        {

            addVertexWithId(vertexId);
        }

        ds_valid = false;
        graphs_built = false;
        all_intersections.clear();
        pll.reset();
        pll_graph.reset();

        for (int edgeId : edgeList)
        {

            if (edgeId < 0 || edgeId >= static_cast<int>(hyperedges.size()))
            {
                std::cerr << "Warning: Invalid edge ID " << edgeId << " in edgeList. Skipping." << std::endl;
                continue;
            }

            auto &edge_vertices = hyperedges[edgeId].vertices;
            if (std::find(edge_vertices.begin(), edge_vertices.end(), vertexId) == edge_vertices.end())
            {
                edge_vertices.push_back(vertexId);
            }

            auto &incident_edges = vertex_to_edges[vertexId];
            if (std::find(incident_edges.begin(), incident_edges.end(), edgeId) == incident_edges.end())
            {
                incident_edges.push_back(edgeId);
            }
        }
    }

    void removeVertexFromEdges(int vertexId, const std::vector<int> &edgeList)
    {

        if (vertexId < 0 || vertexId >= static_cast<int>(vertex_to_edges.size()) || vertex_to_edges[vertexId].empty())
        {

            std::cerr << "Warning: Vertex " << vertexId << " does not exist or has no incident edges. Cannot remove." << std::endl;
            return;
        }

        ds_valid = false;
        graphs_built = false;
        all_intersections.clear();
        pll.reset();
        pll_graph.reset();

        bool vertex_removed_from_any = false;
        for (int edgeId : edgeList)
        {

            if (edgeId < 0 || edgeId >= static_cast<int>(hyperedges.size()))
            {
                std::cerr << "Warning: Invalid edge ID " << edgeId << " in edgeList. Skipping." << std::endl;
                continue;
            }

            auto &edge_vertices = hyperedges[edgeId].vertices;
            auto it = std::find(edge_vertices.begin(), edge_vertices.end(), vertexId);
            if (it != edge_vertices.end())
            {
                edge_vertices.erase(it);
                vertex_removed_from_any = true;
            }
            else
            {
                std::cerr << "Warning: Vertex " << vertexId << " not found in edge " << edgeId << ". Skipping removal from this edge." << std::endl;
            }

            auto &incident_edges = vertex_to_edges[vertexId];
            auto edge_it = std::find(incident_edges.begin(), incident_edges.end(), edgeId);
            if (edge_it != incident_edges.end())
            {
                incident_edges.erase(edge_it);
            }
        }

        if (vertex_removed_from_any)
        {
        }
        else
        {
            std::cout << "Info: Vertex " << vertexId << " was not found in any of the specified valid edges." << std::endl;
        }
    }

    int addVertex()
    {
        ds_valid = false;
        graphs_built = false;
        int vertexId = vertices.size();
        vertices.push_back(vertexId);

        vertex_to_edges.push_back(std::vector<int>());
        return vertexId;
    }

    int addVertices(size_t count)
    {
        ds_valid = false;
        graphs_built = false;
        int firstId = vertices.size();
        vertices.resize(firstId + count);
        vertex_to_edges.resize(firstId + count);

        for (size_t i = 0; i < count; i++)
        {
            vertices[firstId + i] = firstId + i;
        }
        return vertices.size();
    }

    void addVertexWithId(int vertexId)
    {

        ds_valid = false;
        graphs_built = false;

        if (vertexId < 0)
        {
            throw std::invalid_argument("Vertex ID must be non-negative");
        }
        if (vertexId >= static_cast<int>(vertices.size()))
        {

            size_t new_size = std::max(static_cast<size_t>(vertexId + 1), vertices.size() + (vertices.size() / 2) + 1);
            vertices.resize(new_size);
            vertex_to_edges.resize(new_size);

            for (int i = vertices.size() - 1; i >= vertexId && vertices[i] == 0; --i)
            {
                if (i < vertex_to_edges.size())
                {
                    vertices[i] = i;
                }
                else
                {

                    throw std::runtime_error("Internal error: vertex_to_edges size mismatch after resize");
                }
            }

            if (vertexId < vertices.size() && vertices[vertexId] == 0)
            {
                vertices[vertexId] = vertexId;
            }
        }
    }

    void addHyperedgeWithId(int edgeId, const std::vector<int> &vertexList)
    {

        ds_valid = false;
        graphs_built = false;
        all_intersections.clear();

        if (edgeId < 0)
        {
            throw std::invalid_argument("Edge ID must be non-negative");
        }
        if (edgeId >= static_cast<int>(hyperedges.size()))
        {

            hyperedges.resize(edgeId + 1);
        }

        for (int v : vertexList)
        {
            if (v < 0)
            {
                throw std::invalid_argument("Vertex ID in edge cannot be negative");
            }
            if (v >= static_cast<int>(vertex_to_edges.size()))
            {

                addVertexWithId(v);
            }
        }

        hyperedges[edgeId] = Hyperedge(vertexList);

        for (int v : vertexList)
        {

            if (v >= 0 && v < static_cast<int>(vertex_to_edges.size()))
            {

                if (std::find(vertex_to_edges[v].begin(), vertex_to_edges[v].end(), edgeId) == vertex_to_edges[v].end())
                {
                    vertex_to_edges[v].push_back(edgeId);
                }
            }
            else
            {

                throw std::runtime_error("Internal error: Vertex ID out of bounds after addVertexWithId");
            }
        }
    }

    int addHyperedge(const std::vector<int> &vertexList)
    {

        ds_valid = false;
        graphs_built = false;

        int maxVertex = vertices.size() - 1;
        for (int v : vertexList)
        {

            if (v < 0 || v > maxVertex /* || vertex_to_edges[v].empty() if empty means deleted */)
                throw std::invalid_argument("Vertex id does not exist or is invalid");
        }

        int newEdgeId = hyperedges.size();
        hyperedges.emplace_back(Hyperedge(vertexList));

        for (int v : vertexList)
        {
            if (std::find(vertex_to_edges[v].begin(), vertex_to_edges[v].end(), newEdgeId) == vertex_to_edges[v].end())
            {
                vertex_to_edges[v].push_back(newEdgeId);
            }
        }

        all_intersections.clear();

        return newEdgeId;
    }

    void reserveHyperedges(size_t count)
    {
        hyperedges.reserve(hyperedges.size() + count);
    }

    const std::vector<int> &getIncidentHyperedges(int vertexId) const
    {
        if (vertexId < 0 || vertexId >= static_cast<int>(vertices.size()))
            throw std::invalid_argument("Vertex id does not exist");
        return vertex_to_edges[vertexId];
    }

    const Hyperedge &getHyperedge(int edgeId) const
    {
        if (edgeId < 0 || edgeId >= static_cast<int>(hyperedges.size()))
            throw std::invalid_argument("Hyperedge id does not exist");
        return hyperedges[edgeId];
    }

    size_t numVertices() const
    {
        return vertices.size();
    }

    size_t numHyperedges() const
    {
        return hyperedges.size();
    }

    static Hypergraph fromFile(const std::string &filename)
    {
        std::ifstream file(filename);
        if (!file.is_open())
        {
            throw std::runtime_error("Cannot open file: " + filename);
        }

        std::string line;
        int max_vertex_id = -1;
        std::vector<std::vector<int>> edge_data;

        while (std::getline(file, line))
        {
            std::istringstream iss(line);
            int vertex_id;
            std::vector<int> edge;

            while (iss >> vertex_id)
            {
                edge.push_back(vertex_id);
                max_vertex_id = std::max(max_vertex_id, vertex_id);
            }

            if (!edge.empty())
            {
                edge_data.push_back(std::move(edge));
            }
        }

        Hypergraph hg(max_vertex_id + 1, edge_data.size());
        hg.addVertices(max_vertex_id + 1);

        for (const auto &edge : edge_data)
        {
            hg.addHyperedge(edge);
        }

        return hg;
    }

    void offline_industry()
    {
        void offline_industry_baseline();
        void offline_industry_pll();
    }
    void offline_industry_pll(string filename = "")
    {

        pll_graph.release();
        pll_graph = std::make_unique<WeightedGraph>(hyperedges.size());
        for (size_t i = 0; i < hyperedges.size(); i++)
        {
            if (hyperedges[i].size() == 0)
                continue;
            for (size_t j = i + 1; j < hyperedges.size(); j++)
            {
                if (hyperedges[j].size() == 0)
                    continue;
                auto intersection = getHyperedgeIntersection(i, j);
                int size = intersection.size();
                if (size > 0)
                {
                    pll_graph->addEdge(i, j, size);
                }
            }
        }
        pll = std::make_unique<WeightedPrunedLandmarkIndex>(*pll_graph);
        pll->offline_industry(filename);
        pll_total_label_size = pll->getTotalLabelSize();

        pll_memory_bytes = static_cast<size_t>((pll_graph->getAdjListMemoryUsageMB() + pll->getMemoryUsageMB()) * 1024.0 * 1024.0);
    }

    void offline_industry_baseline(const std::string &cache_prefix = "")
    {
        bool hg_ds_loaded = false;
        bool layered_graphs_loaded = false;

        if (!cache_prefix.empty())
        {
            std::string hg_ds_file = cache_prefix + "_hg_ds.idx";
            try
            {

                auto temp_ds = std::make_unique<DisjointSets>(0);
                if (temp_ds->loadFromFile(hg_ds_file))
                {

                    if (temp_ds->parent.size() == vertices.size())
                    {
                        ds = std::move(temp_ds);
                        ds_valid = true;
                        ds_nodes = ds->parent.size();
                        ds_memory_bytes = sizeof(*ds) + ds->parent.capacity() * sizeof(int) + ds->rank.capacity() * sizeof(int);
                        hg_ds_loaded = true;
                        std::cout << "Hypergraph Disjoint Set loaded from cache: " << hg_ds_file << std::endl;
                    }
                    else
                    {
                        std::cerr << "Warning: Hypergraph DS cache size mismatch (" << temp_ds->parent.size()
                                  << ") vs current vertex count (" << vertices.size() << "). Rebuilding DS." << std::endl;
                        ds.reset();
                        ds_valid = false;
                    }
                }
                else
                {
                    std::cout << "Info: Hypergraph DS cache not found or invalid at '" << hg_ds_file << "'. Building DS." << std::endl;
                }
            }
            catch (const std::exception &e)
            {
                std::cerr << "Warning: Failed to load Hypergraph DS cache '" << hg_ds_file << "'. Building DS. Error: " << e.what() << std::endl;
                ds.reset();
                ds_valid = false;
            }
        }

        if (!hg_ds_loaded)
        {
            std::cout << "Building Hypergraph Disjoint Set..." << std::endl;
            buildHypergraphDS();

            if (!cache_prefix.empty())
            {
                std::string hg_ds_file = cache_prefix + "_hg_ds.idx";
                try
                {
                    std::filesystem::path p(hg_ds_file);
                    if (p.has_parent_path())
                    {
                        std::filesystem::create_directories(p.parent_path());
                    }
                    if (ds && ds->saveToFile(hg_ds_file))
                    {
                        std::cout << "Hypergraph Disjoint Set saved to cache: " << hg_ds_file << std::endl;
                    }
                    else if (ds)
                    {
                    }
                }
                catch (const std::exception &e)
                {
                    std::cerr << "Warning: Failed to save Hypergraph DS cache '" << hg_ds_file << "'. Error: " << e.what() << std::endl;
                }
            }
        }

        if (!cache_prefix.empty() && !graphs_built)
        {
            bool loaded_all_layers = true;
            weighted_graphs.clear();
            weighted_graphs.resize(MAX_INTERSECTION_SIZE + 1);
            weighted_graphs_adj_list_memory_bytes.assign(MAX_INTERSECTION_SIZE + 1, 0);
            weighted_graphs_ds_memory_bytes.assign(MAX_INTERSECTION_SIZE + 1, 0);
            weighted_graphs_total_nodes = 0;
            weighted_graphs_total_edges = 0;

            for (int k = 1; k <= MAX_INTERSECTION_SIZE; ++k)
            {
                std::string adj_file = cache_prefix + "_lds_k" + std::to_string(k) + "_adj.idx";
                std::string ds_file = cache_prefix + "_lds_k" + std::to_string(k) + "_ds.idx";

                weighted_graphs[k] = std::make_unique<WeightedGraph>(hyperedges.size(), k);

                if (!weighted_graphs[k]->loadAdjList(adj_file) || !weighted_graphs[k]->loadDisjointSets(ds_file))
                {
                    std::cout << "Info: Layered DS cache not found or invalid for k=" << k << ". Building required." << std::endl;
                    loaded_all_layers = false;
                    weighted_graphs.clear();

                    weighted_graphs_adj_list_memory_bytes.assign(MAX_INTERSECTION_SIZE + 1, 0);
                    weighted_graphs_ds_memory_bytes.assign(MAX_INTERSECTION_SIZE + 1, 0);
                    weighted_graphs_total_nodes = 0;
                    weighted_graphs_total_edges = 0;
                    break;
                }

                weighted_graphs_total_nodes += weighted_graphs[k]->numVertices();
                weighted_graphs_total_edges += weighted_graphs[k]->numEdges();
                weighted_graphs_adj_list_memory_bytes[k] = static_cast<size_t>(weighted_graphs[k]->getAdjListMemoryUsageMB() * 1024.0 * 1024.0);
                weighted_graphs_ds_memory_bytes[k] = static_cast<size_t>(weighted_graphs[k]->getDsMemoryUsageMB() * 1024.0 * 1024.0);
            }
            if (loaded_all_layers)
            {
                std::cout << "All Layered DS levels loaded from cache prefix: " << cache_prefix << std::endl;
                graphs_built = true;
                layered_graphs_loaded = true;
            }
        }

        if (!layered_graphs_loaded && !graphs_built)
        {
            std::cout << "Building Layered DS..." << std::endl;

            if (all_intersections.empty())
            {
                calculateAllIntersectionsParallel();
            }

            weighted_graphs.clear();
            weighted_graphs.resize(MAX_INTERSECTION_SIZE + 1);
            weighted_graphs_adj_list_memory_bytes.assign(MAX_INTERSECTION_SIZE + 1, 0);
            weighted_graphs_ds_memory_bytes.assign(MAX_INTERSECTION_SIZE + 1, 0);
            weighted_graphs_total_nodes = 0;
            weighted_graphs_total_edges = 0;

            for (int min_size = 1; min_size <= MAX_INTERSECTION_SIZE; min_size++)
            {

                weighted_graphs[min_size] = std::make_unique<WeightedGraph>(hyperedges.size(), min_size);

                for (const auto &intersection_info : all_intersections)
                {
                    if (std::get<2>(intersection_info) >= min_size)
                    {
                        weighted_graphs[min_size]->addEdge(std::get<0>(intersection_info), std::get<1>(intersection_info), std::get<2>(intersection_info));
                    }
                }

                weighted_graphs[min_size]->offline_industry();

                weighted_graphs_total_nodes += weighted_graphs[min_size]->numVertices();
                weighted_graphs_total_edges += weighted_graphs[min_size]->numEdges();
                weighted_graphs_adj_list_memory_bytes[min_size] = static_cast<size_t>(weighted_graphs[min_size]->getAdjListMemoryUsageMB() * 1024.0 * 1024.0);
                weighted_graphs_ds_memory_bytes[min_size] = static_cast<size_t>(weighted_graphs[min_size]->getDsMemoryUsageMB() * 1024.0 * 1024.0);

                if (!cache_prefix.empty())
                {
                    std::string adj_file = cache_prefix + "_lds_k" + std::to_string(min_size) + "_adj.idx";
                    std::string ds_file = cache_prefix + "_lds_k" + std::to_string(min_size) + "_ds.idx";
                    try
                    {
                        std::filesystem::path p(adj_file);
                        if (p.has_parent_path())
                        {
                            std::filesystem::create_directories(p.parent_path());
                        }

                        if (!weighted_graphs[min_size]->saveAdjList(adj_file))
                        {
                            std::cerr << "Warning: Failed to save Layered DS adjacency list cache for k=" << min_size << std::endl;
                        }
                        if (!weighted_graphs[min_size]->saveDisjointSets(ds_file))
                        {
                            std::cerr << "Warning: Failed to save Layered DS disjoint set cache for k=" << min_size << std::endl;
                        }
                    }
                    catch (const std::exception &e)
                    {
                        std::cerr << "Warning: Failed to save Layered DS cache for k=" << min_size << ". Error: " << e.what() << std::endl;
                    }
                }
            }
            graphs_built = true;
            if (!cache_prefix.empty())
            {
                std::cout << "Layered DS saved to cache prefix: " << cache_prefix << std::endl;
            }
        }
    }

    bool isReachable(int u, int v) const
    {
        if (u < 0 || u >= static_cast<int>(vertices.size()) ||
            v < 0 || v >= static_cast<int>(vertices.size()))
            throw std::invalid_argument("Vertex id does not exist");

        if (ds_valid && ds)
        {
            return ds->find(u) == ds->find(v);
        }
        else
        {

            DisjointSets temp_ds = buildDisjointSets();
            return temp_ds.find(u) == temp_ds.find(v);
        }
    }

    bool isReachableBidirectionalBFS(int source, int target)
    {
        if (source < 0 || source >= static_cast<int>(vertices.size()) ||
            target < 0 || target >= static_cast<int>(vertices.size()))
            throw std::invalid_argument("Vertex id does not exist");

        if (source == target)
            return true;

        std::vector<bool> visited_forward(vertices.size(), false);
        std::vector<bool> visited_backward(vertices.size(), false);

        std::queue<int> q_forward;
        std::queue<int> q_backward;

        q_forward.push(source);
        q_backward.push(target);
        visited_forward[source] = true;
        visited_backward[target] = true;

        while (!q_forward.empty() && !q_backward.empty())
        {

            int forward_size = q_forward.size();
            for (int i = 0; i < forward_size; i++)
            {
                int current = q_forward.front();
                q_forward.pop();

                for (int edgeId : vertex_to_edges[current])
                {
                    const auto &edge = hyperedges[edgeId];

                    for (int neighbor : edge.vertices)
                    {

                        if (visited_backward[neighbor])
                        {
                            return true;
                        }

                        if (!visited_forward[neighbor])
                        {
                            visited_forward[neighbor] = true;
                            q_forward.push(neighbor);
                        }
                    }
                }
            }

            int backward_size = q_backward.size();
            for (int i = 0; i < backward_size; i++)
            {
                int current = q_backward.front();
                q_backward.pop();

                for (int edgeId : vertex_to_edges[current])
                {
                    const auto &edge = hyperedges[edgeId];

                    for (int neighbor : edge.vertices)
                    {

                        if (visited_forward[neighbor])
                        {
                            return true;
                        }

                        if (!visited_backward[neighbor])
                        {
                            visited_backward[neighbor] = true;
                            q_backward.push(neighbor);
                        }
                    }
                }
            }
        }

        return false;
    }
    bool isReachableBidirectionalBFS(int source, int target, int minIntersectionSize = 0) const
    {
        if (source < 0 || source >= static_cast<int>(vertices.size()) ||
            target < 0 || target >= static_cast<int>(vertices.size()))
            throw std::invalid_argument("Vertex id does not exist");

        if (source == target)
            return true;

        const auto &source_edges = vertex_to_edges[source];
        const auto &target_edges = vertex_to_edges[target];

        for (int edge_id : source_edges)
        {
            if (std::find(target_edges.begin(), target_edges.end(), edge_id) != target_edges.end())
            {
                return true;
            }
        }

        if (minIntersectionSize <= 0)
            return isReachableBidirectionalBFS(source, target);

        std::vector<int> pred_edge_forward(vertices.size(), -1);
        std::vector<int> pred_edge_backward(vertices.size(), -1);
        std::vector<int> distance_forward(vertices.size(), -1);
        std::vector<int> distance_backward(vertices.size(), -1);

        std::queue<int> q_forward;
        std::queue<int> q_backward;

        q_forward.push(source);
        q_backward.push(target);
        distance_forward[source] = 0;
        distance_backward[target] = 0;

        std::vector<std::tuple<int, int, int>> meeting_points;

        while (!q_forward.empty() && !q_backward.empty())
        {

            int forward_size = q_forward.size();
            for (int i = 0; i < forward_size; i++)
            {
                int current = q_forward.front();
                q_forward.pop();

                for (int edge_id : vertex_to_edges[current])
                {
                    const auto &edge = hyperedges[edge_id];

                    if (pred_edge_forward[current] == -1 ||
                        getHyperedgeIntersection(pred_edge_forward[current], edge_id).size() >=
                            static_cast<size_t>(minIntersectionSize))
                    {

                        for (int next : edge.vertices)
                        {

                            if (distance_forward[next] == -1)
                            {
                                distance_forward[next] = distance_forward[current] + 1;
                                pred_edge_forward[next] = edge_id;
                                q_forward.push(next);

                                if (distance_backward[next] != -1)
                                {
                                    meeting_points.emplace_back(next, edge_id, pred_edge_backward[next]);
                                }
                            }
                        }
                    }
                }
            }

            int backward_size = q_backward.size();
            for (int i = 0; i < backward_size; i++)
            {
                int current = q_backward.front();
                q_backward.pop();

                for (int edge_id : vertex_to_edges[current])
                {
                    const auto &edge = hyperedges[edge_id];

                    if (pred_edge_backward[current] == -1 ||
                        getHyperedgeIntersection(pred_edge_backward[current], edge_id).size() >=
                            static_cast<size_t>(minIntersectionSize))
                    {

                        for (int next : edge.vertices)
                        {

                            if (distance_backward[next] == -1)
                            {
                                distance_backward[next] = distance_backward[current] + 1;
                                pred_edge_backward[next] = edge_id;
                                q_backward.push(next);

                                if (distance_forward[next] != -1)
                                {
                                    meeting_points.emplace_back(next, pred_edge_forward[next], edge_id);
                                }
                            }
                        }
                    }
                }
            }

            for (const auto &[meet_point, forward_edge, backward_edge] : meeting_points)
            {

                if (forward_edge == -1 || backward_edge == -1)
                {
                    continue;
                }

                auto intersection = getHyperedgeIntersection(forward_edge, backward_edge);
                if (intersection.size() >= static_cast<size_t>(minIntersectionSize))
                {
                    return true;
                }
            }

            meeting_points.clear();
        }

        return false;
    }

    bool isReachableBidirectionalBFSByEdge(int sourceEdge, int targetEdge, int minIntersectionSize = 0) const
    {

        if (sourceEdge < 0 || sourceEdge >= static_cast<int>(hyperedges.size()) ||
            targetEdge < 0 || targetEdge >= static_cast<int>(hyperedges.size()))
            return false;

        if (sourceEdge == targetEdge)
            return true;

        auto direct_intersection = getHyperedgeIntersection(sourceEdge, targetEdge);
        if (direct_intersection.size() >= static_cast<size_t>(minIntersectionSize))
        {
            return true;
        }

        const auto &source_vertices = hyperedges[sourceEdge].vertices;
        const auto &target_vertices = hyperedges[targetEdge].vertices;

        if (source_vertices.empty() || target_vertices.empty())
        {
            return false;
        }

        if (minIntersectionSize <= 0)
        {

            return isReachableBidirectionalBFS(source_vertices[0], target_vertices[0]);
        }

        std::vector<int> pred_edge_forward(vertices.size(), -1);
        std::vector<int> pred_edge_backward(vertices.size(), -1);
        std::vector<int> distance_forward(vertices.size(), -1);
        std::vector<int> distance_backward(vertices.size(), -1);

        std::queue<int> q_forward;
        std::queue<int> q_backward;

        for (int vertex : source_vertices)
        {
            q_forward.push(vertex);
            distance_forward[vertex] = 0;
            pred_edge_forward[vertex] = sourceEdge;
        }

        for (int vertex : target_vertices)
        {
            q_backward.push(vertex);
            distance_backward[vertex] = 0;
            pred_edge_backward[vertex] = targetEdge;
        }

        std::vector<std::tuple<int, int, int>> meeting_points;

        while (!q_forward.empty() && !q_backward.empty())
        {

            int forward_size = q_forward.size();
            for (int i = 0; i < forward_size; i++)
            {
                int current = q_forward.front();
                q_forward.pop();

                for (int edge_id : vertex_to_edges[current])
                {
                    const auto &edge = hyperedges[edge_id];

                    if (pred_edge_forward[current] == -1 ||
                        getHyperedgeIntersection(pred_edge_forward[current], edge_id).size() >=
                            static_cast<size_t>(minIntersectionSize))
                    {

                        for (int next : edge.vertices)
                        {

                            if (distance_forward[next] == -1)
                            {
                                distance_forward[next] = distance_forward[current] + 1;
                                pred_edge_forward[next] = edge_id;
                                q_forward.push(next);

                                if (distance_backward[next] != -1)
                                {
                                    meeting_points.emplace_back(next, edge_id, pred_edge_backward[next]);
                                }
                            }
                        }
                    }
                }
            }

            int backward_size = q_backward.size();
            for (int i = 0; i < backward_size; i++)
            {
                int current = q_backward.front();
                q_backward.pop();

                for (int edge_id : vertex_to_edges[current])
                {
                    const auto &edge = hyperedges[edge_id];

                    if (pred_edge_backward[current] == -1 ||
                        getHyperedgeIntersection(pred_edge_backward[current], edge_id).size() >=
                            static_cast<size_t>(minIntersectionSize))
                    {

                        for (int next : edge.vertices)
                        {

                            if (distance_backward[next] == -1)
                            {
                                distance_backward[next] = distance_backward[current] + 1;
                                pred_edge_backward[next] = edge_id;
                                q_backward.push(next);

                                if (distance_forward[next] != -1)
                                {
                                    meeting_points.emplace_back(next, pred_edge_forward[next], edge_id);
                                }
                            }
                        }
                    }
                }
            }

            for (const auto &[meet_point, forward_edge, backward_edge] : meeting_points)
            {

                auto intersection = getHyperedgeIntersection(forward_edge, backward_edge);
                if (intersection.size() >= static_cast<size_t>(minIntersectionSize))
                {
                    return true;
                }
            }

            meeting_points.clear();
        }

        return false;
    }

    bool isReachableViaWeightedGraph(int sourceVertex, int targetVertex, int minIntersectionSize = 1)
    {
        if (sourceVertex < 0 || sourceVertex >= static_cast<int>(vertices.size()) ||
            targetVertex < 0 || targetVertex >= static_cast<int>(vertices.size()))
            throw std::invalid_argument("Vertex id does not exist");

        if (sourceVertex == targetVertex)
            return true;

        if (!graphs_built)
        {
            const_cast<Hypergraph *>(this)->offline_industry_baseline();
        }

        int effective_min_size = minIntersectionSize;
        if (effective_min_size <= 0)
        {

            effective_min_size = 1;
        }
        if (effective_min_size > MAX_INTERSECTION_SIZE)
        {

            effective_min_size = MAX_INTERSECTION_SIZE;
        }

        if (effective_min_size < 1 || effective_min_size > MAX_INTERSECTION_SIZE || !weighted_graphs[effective_min_size])
        {

            if (!weighted_graphs[effective_min_size])
            {
                throw std::runtime_error("Requested WeightedGraph layer " + std::to_string(effective_min_size) + " is not built.");
            }
        }

        const auto &source_edges = vertex_to_edges[sourceVertex];
        const auto &target_edges = vertex_to_edges[targetVertex];

        if (source_edges.empty() || target_edges.empty())
        {
            return false;
        }

        for (int source_edge : source_edges)
        {
            if (std::find(target_edges.begin(), target_edges.end(), source_edge) != target_edges.end())
            {
                return true;
            }
        }

        for (int source_edge : source_edges)
        {
            for (int target_edge : target_edges)
            {

                if (weighted_graphs[effective_min_size]->disjointSet_reachability_query(source_edge, target_edge))
                {
                    return true;
                }
            }
        }

        return false;
    }

    bool isReachableViaWeightedGraphByEdge(int sourceEdge, int targetEdge, int minIntersectionSize = 1)
    {

        if (sourceEdge < 0 || sourceEdge >= static_cast<int>(hyperedges.size()) ||
            targetEdge < 0 || targetEdge >= static_cast<int>(hyperedges.size()))
            return false;

        if (sourceEdge == targetEdge)
            return true;

        if (!graphs_built)
        {
            const_cast<Hypergraph *>(this)->offline_industry_baseline();
        }

        int effective_min_size = minIntersectionSize;
        if (effective_min_size <= 0)
        {

            effective_min_size = 1;
        }
        if (effective_min_size > MAX_INTERSECTION_SIZE)
        {

            effective_min_size = MAX_INTERSECTION_SIZE;
        }

        if (effective_min_size < 1 || effective_min_size > MAX_INTERSECTION_SIZE || !weighted_graphs[effective_min_size])
        {

            if (!weighted_graphs[effective_min_size])
            {
                throw std::runtime_error("Requested WeightedGraph layer " + std::to_string(effective_min_size) + " is not built.");
            }
        }

        if (weighted_graphs[effective_min_size]->disjointSet_reachability_query(sourceEdge, targetEdge))
            return true;
        return false;
    }

    bool isReachableViaUWeightedPLL(int sourceVertex, int targetVertex, int minIntersectionSize = 0) const
    {
        if (sourceVertex < 0 || sourceVertex >= static_cast<int>(vertices.size()) ||
            targetVertex < 0 || targetVertex >= static_cast<int>(vertices.size()))
            throw std::invalid_argument("Vertex id does not exist");

        if (sourceVertex == targetVertex)
            return true;

        if (!graphs_built || !pll)
        {
            const_cast<Hypergraph *>(this)->offline_industry();
            if (!pll)
            {
                throw std::runtime_error("PLL index failed to build.");
            }
        }

        const auto &source_edges = vertex_to_edges[sourceVertex];
        const auto &target_edges = vertex_to_edges[targetVertex];

        if (source_edges.empty() || target_edges.empty())
        {
            return false;
        }

        for (int source_edge : source_edges)
        {
            if (std::find(target_edges.begin(), target_edges.end(), source_edge) != target_edges.end())
            {
                return true;
            }
        }

        for (int source_edge : source_edges)
        {
            for (int target_edge : target_edges)
            {
                if (pll->reachability_query(source_edge, target_edge, minIntersectionSize))
                {
                    return true;
                }
            }
        }

        return false;
    }

    bool isReachableViaUWeightedPLLByEdge(int sourceEdge, int targetEdge, int minIntersectionSize = 0) const
    {

        if (sourceEdge == targetEdge)
            return true;

        if (!graphs_built || !pll)
        {
            const_cast<Hypergraph *>(this)->offline_industry();
            if (!pll)
            {
                throw std::runtime_error("PLL index failed to build.");
            }
        }

        if (pll->reachability_query(sourceEdge, targetEdge, minIntersectionSize))
        {
            return true;
        }

        return false;
    }

    std::vector<std::vector<int>> getConnectedComponents()
    {

        DisjointSets *working_ds = nullptr;
        DisjointSets temp_ds(vertices.size());

        if (ds_valid && ds)
        {

            working_ds = ds.get();
        }
        else
        {

            temp_ds = buildDisjointSets();
            working_ds = &temp_ds;
        }

        std::vector<int> component_size(vertices.size(), 0);
        for (size_t i = 0; i < vertices.size(); i++)
        {
            int root = working_ds->find(i);
            component_size[root]++;
        }

        int num_components = 0;
        for (size_t i = 0; i < vertices.size(); i++)
        {
            if (component_size[i] > 0)
            {
                num_components++;
            }
        }

        std::vector<std::vector<int>> result(num_components);
        std::vector<int> component_index(vertices.size(), -1);

        int curr_idx = 0;
        for (size_t i = 0; i < vertices.size(); i++)
        {
            if (component_size[i] > 0)
            {
                component_index[i] = curr_idx;
                result[curr_idx].reserve(component_size[i]);
                curr_idx++;
            }
        }

        for (size_t i = 0; i < vertices.size(); i++)
        {
            int root = working_ds->find(i);
            result[component_index[root]].push_back(i);
        }

        return result;
    }

    std::unique_ptr<WeightedGraph> convertToWeightedGraph() const
    {

        auto graph = std::make_unique<WeightedGraph>(hyperedges.size());

        for (size_t i = 0; i < hyperedges.size(); i++)
        {
            for (size_t j = i + 1; j < hyperedges.size(); j++)
            {

                auto intersection = getHyperedgeIntersection(i, j);
                if (!intersection.empty())
                {

                    graph->addEdge(i, j, intersection.size());
                }
            }
        }

        return graph;
    }

    std::vector<std::pair<int, int>> getVertexToGraphEdges(int vertexId, const WeightedGraph &graph) const
    {
        if (vertexId < 0 || vertexId >= static_cast<int>(vertices.size()))
            throw std::invalid_argument("Vertex id does not exist");

        std::vector<std::pair<int, int>> result;

        const auto &incident_edges = vertex_to_edges[vertexId];

        for (size_t i = 0; i < incident_edges.size(); i++)
        {
            for (size_t j = i + 1; j < incident_edges.size(); j++)
            {
                int edge1 = incident_edges[i];
                int edge2 = incident_edges[j];

                result.emplace_back(edge1, edge2);
            }
        }

        return result;
    }

    int getHyperedgeToGraphVertex(int edgeId) const
    {
        if (edgeId < 0 || edgeId >= static_cast<int>(hyperedges.size()))
            throw std::invalid_argument("Hyperedge id does not exist");

        return edgeId;
    }

    size_t getDsNodes() const
    {
        return ds_nodes;
    }

    size_t getPllTotalLabelSize() const
    {
        return pll_total_label_size;
    }

    size_t getWeightedGraphsTotalNodes() const
    {
        return weighted_graphs_total_nodes;
    }

    size_t getWeightedGraphsTotalEdges() const
    {
        return weighted_graphs_total_edges;
    }

    double getDsMemoryUsageMB() const
    {
        return static_cast<double>(ds_memory_bytes) / (1024.0 * 1024.0);
    }

    double getWeightedGraphAdjListMemoryUsageMB(int layer) const
    {
        if (layer < 0 || layer > MAX_INTERSECTION_SIZE || layer >= weighted_graphs_adj_list_memory_bytes.size())
        {
            return 0.0;
        }
        return static_cast<double>(weighted_graphs_adj_list_memory_bytes[layer]) / (1024.0 * 1024.0);
    }

    double getWeightedGraphDsMemoryUsageMB(int layer) const
    {
        if (layer < 0 || layer > MAX_INTERSECTION_SIZE || layer >= weighted_graphs_ds_memory_bytes.size())
        {
            return 0.0;
        }
        return static_cast<double>(weighted_graphs_ds_memory_bytes[layer]) / (1024.0 * 1024.0);
    }

    double getWeightedGraphsMemoryUsageMB() const
    {
        size_t total_bytes = 0;
        for (size_t bytes : weighted_graphs_adj_list_memory_bytes)
        {
            total_bytes += bytes;
        }
        for (size_t bytes : weighted_graphs_ds_memory_bytes)
        {
            total_bytes += bytes;
        }
        return static_cast<double>(total_bytes) / (1024.0 * 1024.0);
    }

    double getPllMemoryUsageMB() const
    {
        return static_cast<double>(pll_memory_bytes) / (1024.0 * 1024.0);
    }

    void deleteVertex(int vertexId)
    {

        if (vertexId < 0 || vertexId >= static_cast<int>(vertex_to_edges.size()))
        {
            throw std::invalid_argument("Vertex id is out of bounds for deletion");
        }

        bool already_deleted = true;
        if (vertexId < vertex_to_edges.size() && !vertex_to_edges[vertexId].empty())
        {
            already_deleted = false;
        }

        if (vertex_to_edges[vertexId].empty() && /* 检查 vertices[vertexId] 是否也表示不存在 */ true)
        {

            return;
        }

        ds_valid = false;
        graphs_built = false;

        pll.reset();
        pll_graph.reset();

        all_intersections.clear();

        std::vector<int> incident_edges_copy = vertex_to_edges[vertexId];
        for (int edgeId : incident_edges_copy)
        {
            if (edgeId >= 0 && edgeId < static_cast<int>(hyperedges.size()))
            {
                auto &edge_verts = hyperedges[edgeId].vertices;

                edge_verts.erase(std::remove(edge_verts.begin(), edge_verts.end(), vertexId), edge_verts.end());
            }
        }

        vertex_to_edges[vertexId].clear();
    }

    void buildHypergraphDS()
    {
        ds = std::make_unique<DisjointSets>(vertices.size());
        for (const auto &edge : hyperedges)
        {
            if (edge.vertices.size() > 1)
            {
                int first = edge.vertices[0];
                for (size_t i = 1; i < edge.vertices.size(); i++)
                {
                    ds->merge(first, edge.vertices[i]);
                }
            }
        }
        ds_valid = true;
        ds_nodes = ds->parent.size();
        ds_memory_bytes = sizeof(*ds) + ds->parent.capacity() * sizeof(int) + ds->rank.capacity() * sizeof(int);
    }

    void calculateAllIntersectionsParallel()
    {

        all_intersections.clear();

        std::cout << "Calculating hyperedge intersections..." << std::endl;
        std::mutex intersections_mutex;
        size_t num_edges = hyperedges.size();
        unsigned int num_threads = std::thread::hardware_concurrency();
        if (num_threads == 0)
            num_threads = 1;
        std::vector<std::thread> threads(num_threads);
        all_intersections.reserve(num_edges);

        auto calculate_intersections_task =
            [&](size_t start_idx, size_t end_idx)
        {
            std::vector<std::tuple<size_t, size_t, int>> local_intersections;
            local_intersections.reserve(num_edges);
            for (size_t i = start_idx; i < end_idx; ++i)
            {
                if (hyperedges[i].size() == 0)
                    continue;
                for (size_t j = i + 1; j < num_edges; ++j)
                {
                    if (hyperedges[j].size() == 0)
                        continue;
                    auto intersection = getHyperedgeIntersection(i, j);
                    int size = intersection.size();
                    if (size > 0)
                    {
                        local_intersections.emplace_back(i, j, size);
                    }
                }
            }
            std::lock_guard<std::mutex> lock(intersections_mutex);
            all_intersections.insert(all_intersections.end(),
                                     local_intersections.begin(),
                                     local_intersections.end());
        };

        size_t chunk_size = (num_edges + num_threads - 1) / num_threads;
        size_t current_start = 0;
        for (unsigned int t = 0; t < num_threads; ++t)
        {
            size_t current_end = std::min(current_start + chunk_size, num_edges);
            if (current_start >= current_end)
                break;
            threads[t] = std::thread(calculate_intersections_task, current_start, current_end);
            current_start = current_end;
        }

        for (unsigned int t = 0; t < threads.size(); ++t)
        {
            if (threads[t].joinable())
            {
                threads[t].join();
            }
        }
        std::cout << "Intersection calculation complete. Found " << all_intersections.size() << " intersections." << std::endl;
    }

    static const int MAX_INTERSECTION_SIZE = 10;
    std::unique_ptr<WeightedGraph> pll_graph;
    std::vector<std::tuple<size_t, size_t, int>> all_intersections;

    std::vector<int> vertices;
    std::vector<Hyperedge> hyperedges;
    std::vector<std::vector<int>> vertex_to_edges;

    std::unique_ptr<DisjointSets> ds;

    bool ds_valid;

    std::vector<std::unique_ptr<WeightedGraph>> weighted_graphs;

    std::unique_ptr<WeightedPrunedLandmarkIndex> pll;

    bool graphs_built = false;

    size_t ds_nodes = 0;
    size_t pll_total_label_size = 0;
    size_t weighted_graphs_total_nodes = 0;
    size_t weighted_graphs_total_edges = 0;

    size_t ds_memory_bytes = 0;
    std::vector<size_t> weighted_graphs_adj_list_memory_bytes;
    std::vector<size_t> weighted_graphs_ds_memory_bytes;
    size_t pll_memory_bytes = 0;

private:
};

#endif