#pragma once
#include "WeightedGraph.h"
#include <vector>
#include <queue>
#include <algorithm>
#include <numeric>
#include <limits>
#include <iostream>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>

class WeightedPrunedLandmarkIndex
{
public:
    WeightedPrunedLandmarkIndex(const WeightedGraph &graph, int weightThreshold = 0)
        : g(graph), weightThreshold(weightThreshold)
    {
        label.resize(g.numVertices());
    }

    void offline_industry(std::string cache_path = "");

    bool reachability_query(int u, int v, int queryThreshold) const;

    size_t getTotalLabelSize() const
    {
        size_t total_size = 0;
        for (const auto &vec : label)
        {
            total_size += vec.size();
        }
        return total_size;
    }

    const std::vector<std::vector<std::pair<int, int>>> &getLabels() const
    {
        return label;
    }

    double getMemoryUsageMB() const
    {
        size_t memory_bytes = 0;
        memory_bytes += sizeof(label);
        for (const auto &vec : label)
        {
            memory_bytes += sizeof(vec);

            memory_bytes += vec.capacity() * sizeof(std::pair<int, int>);
        }
        return static_cast<double>(memory_bytes) / (1024.0 * 1024.0);
    }

    int weightThreshold;

    std::vector<std::vector<std::pair<int, int>>> label;

private:
    const WeightedGraph &g;

    void insertOrUpdateLabel(std::vector<std::pair<int, int>> &vec, int landmark, int bw) const;

    bool hopQuery(int curLandmark, int node, int candidateBW) const;

    void prunedBFS(int curLandmark);

    bool intersectWithThreshold(const std::vector<std::pair<int, int>> &a,
                                const std::vector<std::pair<int, int>> &b,
                                int threshold) const;

    bool saveLabels(const std::string &filename) const
    {
        std::ofstream file(filename);
        if (!file.is_open())
        {
            std::cerr << "Error: Cannot open file for writing PLL labels: " << filename << std::endl;
            return false;
        }
        file << label.size() << "\n";
        for (size_t u = 0; u < label.size(); ++u)
        {
            file << u << " " << label[u].size();
            for (const auto &entry : label[u])
            {
                file << " " << entry.first << " " << entry.second;
            }
            file << "\n";
        }
        return file.good();
    }

    bool loadLabels(const std::string &filename)
    {
        std::ifstream file(filename);
        if (!file.is_open())
        {

            return false;
        }

        size_t num_vertices_expected;
        file >> num_vertices_expected;
        if (file.fail() || num_vertices_expected != g.numVertices())
        {
            std::cerr << "Error: PLL cache file vertex count mismatch or read error." << std::endl;
            return false;
        }

        label.assign(num_vertices_expected, {});

        std::string line;
        std::getline(file, line);
        int u, count, landmark, bottleneck;
        while (std::getline(file, line))
        {
            std::istringstream iss(line);
            iss >> u >> count;
            if (iss.fail() || u < 0 || u >= num_vertices_expected)
            {
                std::cerr << "Warning: Skipping invalid line in PLL cache: " << line << std::endl;
                continue;
            }
            label[u].reserve(count);
            bool read_error = false;
            for (int i = 0; i < count; ++i)
            {
                if (iss >> landmark >> bottleneck)
                {
                    label[u].emplace_back(landmark, bottleneck);
                }
                else
                {
                    std::cerr << "Error: Incorrect number of label entries for vertex " << u << " in PLL cache." << std::endl;
                    label[u].clear();
                    read_error = true;
                    break;
                }
            }
            if (read_error)
                continue;

            std::sort(label[u].begin(), label[u].end());
        }

        if (file.bad())
        {
            std::cerr << "Error: File read error occurred during PLL cache loading." << std::endl;
            return false;
        }

        return true;
    }
};

void WeightedPrunedLandmarkIndex::offline_industry(std::string cache_path /* = "" */)
{

    bool loaded_from_cache = false;
    if (!cache_path.empty())
    {
        cache_path = cache_path + "_pll.idx";
        try
        {
            loaded_from_cache = loadLabels(cache_path);
            if (loaded_from_cache)
            {
                std::cout << "PLL index loaded from cache: " << cache_path << std::endl;
                return;
            }
            else
            {
                std::cout << "Info: PLL cache not found or invalid at '" << cache_path << "'. Building index." << std::endl;
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Warning: Failed to load PLL cache '" << cache_path << "'. Building index. Error: " << e.what() << std::endl;
        }
    }

    std::cout << "Building PLL index..." << std::endl;
    int n = g.numVertices();
    if (n == 0)
    {
        std::cout << "Warning: Graph has no vertices. PLL index is empty." << std::endl;
        return;
    }

    label.assign(n, {});

    std::vector<std::pair<int, int>> deg(n);
    for (int i = 0; i < n; ++i)
    {
        deg[i] = {static_cast<int>(g.getNeighbors(i).size()), i};
    }
    std::sort(deg.begin(), deg.end(), [](const auto &a, const auto &b)
              { return a.first > b.first; });
    std::vector<int> order(n);
    for (int i = 0; i < n; ++i)
    {
        order[i] = deg[i].second;
    }

    for (int landmark : order)
    {
        prunedBFS(landmark);
    }

    for (int i = 0; i < n; ++i)
    {

        insertOrUpdateLabel(label[i], i, std::numeric_limits<int>::max());
    }
    std::cout << "PLL index build complete." << std::endl;

    if (!loaded_from_cache && !cache_path.empty())
    {
        try
        {

            std::filesystem::path p(cache_path);
            if (p.has_parent_path())
            {
                std::filesystem::create_directories(p.parent_path());
            }
            if (saveLabels(cache_path))
            {
                std::cout << "PLL index saved to cache: " << cache_path << std::endl;
            }
            else
            {
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Warning: Failed to save PLL cache '" << cache_path << "'. Error: " << e.what() << std::endl;
        }
    }
}

void WeightedPrunedLandmarkIndex::insertOrUpdateLabel(std::vector<std::pair<int, int>> &vec,
                                                      int landmark, int bw) const
{
    auto it = std::lower_bound(vec.begin(), vec.end(), landmark,
                               [](const std::pair<int, int> &p, int val)
                               {
                                   return p.first < val;
                               });
    if (it != vec.end() && it->first == landmark)
    {

        if (bw > it->second)
        {
            it->second = bw;
        }
    }
    else
    {
        vec.insert(it, {landmark, bw});
    }
}

bool WeightedPrunedLandmarkIndex::hopQuery(int curLandmark, int node, int candidateBW) const
{

    if (node == curLandmark)
        return true;

    const auto &vec1 = label[curLandmark];
    const auto &vec2 = label[node];
    int i = 0, j = 0;
    while (i < (int)vec1.size() && j < (int)vec2.size())
    {
        if (vec1[i].first == vec2[j].first)
        {
            int bw = std::min(vec1[i].second, vec2[j].second);
            if (bw >= candidateBW)
                return true;
            ++i;
            ++j;
        }
        else if (vec1[i].first < vec2[j].first)
        {
            ++i;
        }
        else
        {
            ++j;
        }
    }
    return false;
}

bool WeightedPrunedLandmarkIndex::intersectWithThreshold(
    const std::vector<std::pair<int, int>> &a,
    const std::vector<std::pair<int, int>> &b,
    int threshold) const
{
    int i = 0, j = 0;
    while (i < (int)a.size() && j < (int)b.size())
    {
        if (a[i].first == b[j].first)
        {
            int bw = std::min(a[i].second, b[j].second);
            if (bw >= threshold)
                return true;
            ++i;
            ++j;
        }
        else if (a[i].first < b[j].first)
        {
            ++i;
        }
        else
        {
            ++j;
        }
    }
    return false;
}

void WeightedPrunedLandmarkIndex::prunedBFS(int curLandmark)
{
    int n = g.numVertices();

    std::vector<bool> visited(n, false);

    std::queue<std::pair<int, int>> q;
    visited[curLandmark] = true;
    q.push({curLandmark, std::numeric_limits<int>::max()});

    while (!q.empty())
    {
        auto [current, bw] = q.front();
        q.pop();

        if (current != curLandmark && hopQuery(curLandmark, current, bw))
            continue;
        if (current != curLandmark)
            insertOrUpdateLabel(label[current], curLandmark, bw);

        for (const auto &edge : g.getNeighbors(current))
        {
            int v = edge.first;
            int w = edge.second;
            if (!visited[v])
            {
                visited[v] = true;
                int new_bw = std::min(bw, w);
                q.push({v, new_bw});
            }
        }
    }
}

bool WeightedPrunedLandmarkIndex::reachability_query(int u, int v, int queryThreshold) const
{
    if (u < 0 || u >= (int)g.numVertices() || v < 0 || v >= (int)g.numVertices())
        throw std::invalid_argument("Vertex id out of range");
    if (u == v)
        return true;

    return intersectWithThreshold(label[u], label[v], queryThreshold);
}
