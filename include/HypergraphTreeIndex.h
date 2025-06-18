#ifndef HYPERGRAPH_TREE_INDEX_H
#define HYPERGRAPH_TREE_INDEX_H

#include "Hypergraph.h"
#include <vector>
#include <memory>
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <tuple>
#include <fstream>
#include <string>
#include <iomanip>
#include <sstream>
#include <thread>
#include <mutex>
#include <iostream>
#include <chrono>
#include <filesystem>
#include <numeric>
#include <set>
#include <unordered_set>
#include <atomic>
#include <map>
#include <random>
#include <unordered_map>

class HypergraphTreeIndex;

struct TreeNode
{
    bool is_leaf = false;
    int intersection_size = 0;
    std::vector<int> intersection_vertices;
    int hyperedge_id = -1;
    std::vector<std::shared_ptr<TreeNode>> children;
    std::weak_ptr<TreeNode> parent;
    int node_id = -1;
    int depth = 0;

    TreeNode(int id) : node_id(id) {}
};

using TreeNodePtr = std::shared_ptr<TreeNode>;

class HypergraphTreeIndex
{
public:
    HypergraphTreeIndex(Hypergraph &hg) : hypergraph_(hg), next_node_id_(0)
    {
        size_t num_hyperedges = hypergraph_.numHyperedges();

        nodes_.reserve(num_hyperedges * 2);
        up_.reserve(num_hyperedges * 2);
        hyperedge_to_leaf_.resize(num_hyperedges, nullptr);

        if (max_intersection_k_ >= 1)
        {
            edge_neighbors_by_size_.resize(num_hyperedges);
            for (size_t i = 0; i < num_hyperedges; ++i)
            {

                edge_neighbors_by_size_[i].resize(max_intersection_k_ + 10);
            }
        }
    }

    bool querySizeOnly(int u, int v, int k)
    {
        if (k < 1)
            return true;
        if (u == v)
            return true;

        const auto &edges_u = hypergraph_.getIncidentHyperedges(u);
        const auto &edges_v = hypergraph_.getIncidentHyperedges(v);

        if (edges_u.empty() || edges_v.empty())
            return false;

        std::vector<int> sorted_edges_v = edges_v;
        std::sort(sorted_edges_v.begin(), sorted_edges_v.end());
        for (int edge_id_u : edges_u)
        {
            if (std::binary_search(sorted_edges_v.begin(), sorted_edges_v.end(), edge_id_u))
            {
                return true;
            }
        }

        for (int edge_id_u : edges_u)
        {
            if (edge_id_u < 0 || edge_id_u >= hyperedge_to_leaf_.size() || !hyperedge_to_leaf_[edge_id_u])
                continue;
            TreeNodePtr leaf_u = hyperedge_to_leaf_[edge_id_u];

            for (int edge_id_v : edges_v)
            {
                if (edge_id_v < 0 || edge_id_v >= hyperedge_to_leaf_.size() || !hyperedge_to_leaf_[edge_id_v])
                    continue;
                TreeNodePtr leaf_v = hyperedge_to_leaf_[edge_id_v];

                TreeNodePtr lca_node = find_lca(leaf_u, leaf_v);

                if (lca_node && lca_node->intersection_size >= k)
                {
                    return true;
                }
            }
        }
        return false;
    }

    bool queryByEdgeId(int edge_id_u, int edge_id_v, int k)
    {
        if (k < 1)
            return true;
        if (edge_id_u == edge_id_v)
            return true;

        if (edge_id_u < 0 || edge_id_u >= hyperedge_to_leaf_.size() || !hyperedge_to_leaf_[edge_id_u] ||
            edge_id_v < 0 || edge_id_v >= hyperedge_to_leaf_.size() || !hyperedge_to_leaf_[edge_id_v])
        {

            return false;
        }

        TreeNodePtr leaf_u = hyperedge_to_leaf_[edge_id_u];
        TreeNodePtr leaf_v = hyperedge_to_leaf_[edge_id_v];

        TreeNodePtr lca_node = find_lca(leaf_u, leaf_v);

        if (lca_node && lca_node->intersection_size >= k)
        {
            return true;
        }

        return false;
    }

    void buildIndexCacheSizeOnly(std::string cache_path = "")
    {

        bool loaded_from_cache = false;
        if (!cache_path.empty())
        {
            cache_path = cache_path + "hypergraph_tree_index_SizeOnly";
            try
            {
                loaded_from_cache = loadIndex(cache_path);
                if (loaded_from_cache)
                {
                    std::cout << "HypergraphTreeIndex (SizeOnly) loaded from cache: " << cache_path << std::endl;
#ifdef DEBUG
                    std::cout << "Recalculating neighbor information (SizeOnly) after loading from cache..." << std::endl;
#endif
                    recalculateAndPopulateNeighborsSizeOnly();
#ifdef DEBUG
                    std::cout << "Neighbor information recalculated." << std::endl;
#endif

                    return;
                }
                else
                {
                    std::cout << "Info: HypergraphTreeIndex (SizeOnly) cache not found or invalid at '" << cache_path << "'. Building index." << std::endl;
                }
            }
            catch (const std::exception &e)
            {
                std::cerr << "Warning: Failed to load HypergraphTreeIndex (SizeOnly) cache '" << cache_path << "'. Building index. Error: " << e.what() << std::endl;
            }
        }

        std::cout << "Building HypergraphTreeIndex (SizeOnly)..." << std::endl;
        int num_hyperedges = hypergraph_.numHyperedges();
        if (num_hyperedges == 0)
        {
            std::cout << "Warning: Hypergraph has no hyperedges. Index is empty." << std::endl;
            return;
        }

        std::vector<TreeNodePtr> current_leaf_nodes;
        current_leaf_nodes.resize(num_hyperedges);
        for (int i = 0; i < num_hyperedges; ++i)
        {
            if (hypergraph_.getHyperedge(i).size() == 0)
                continue;
            TreeNodePtr leaf = std::make_shared<TreeNode>(next_node_id_++);
            leaf->is_leaf = true;
            leaf->hyperedge_id = i;
            nodes_.push_back(leaf);
            if (i < hyperedge_to_leaf_.size())
            {
                hyperedge_to_leaf_[i] = leaf;
            }
            else
            {
                throw std::out_of_range("Hyperedge ID out of range for hyperedge_to_leaf_ vector during build (SizeOnly)");
            }
            current_leaf_nodes.push_back(leaf);
        }

        if (current_leaf_nodes.empty())
        {
            std::cout << "Warning: No valid hyperedges found. Index is empty." << std::endl;
            return;
        }

        std::vector<std::tuple<int, int, int>> merge_candidates_size_only;
        calculateAndPrepareMergeCandidatesSizeOnly(current_leaf_nodes, merge_candidates_size_only);

        sortMergeCandidatesSizeOnly(merge_candidates_size_only);

        populateEdgeNeighborsSizeOnly(merge_candidates_size_only);

        initializeDSU();

        mergeNodesSizeOnly(merge_candidates_size_only);

        std::vector<TreeNodePtr> actual_roots;
        std::vector<bool> is_root_added(next_node_id_, false);

        for (int i = 0; i < next_node_id_; ++i)
        {

            if (i < dsu_parent_.size() && dsu_parent_[i] == i)
            {

                if (i < nodes_.size() && nodes_[i])
                {
                    int root_id = i;
                    if (root_id < is_root_added.size() && !is_root_added[root_id])
                    {
                        actual_roots.push_back(nodes_[root_id]);
                        is_root_added[root_id] = true;
                    }
                }
                else
                {
#ifdef DEBUG
                    std::cerr << "Warning: DSU root ID " << i << " does not correspond to a valid node." << std::endl;
#endif
                }
            }
        }

        roots_.clear();

        if (!actual_roots.empty())
        {

            int virtual_root_id = next_node_id_++;
            TreeNodePtr virtual_root = std::make_shared<TreeNode>(virtual_root_id);
            virtual_root->is_leaf = false;
            virtual_root->intersection_size = 0;

            if (virtual_root_id >= nodes_.size())
                nodes_.resize(virtual_root_id + 1);
            nodes_[virtual_root_id] = virtual_root;
            if (virtual_root_id >= dsu_parent_.size())
                dsu_parent_.resize(virtual_root_id + 1);
            dsu_parent_[virtual_root_id] = virtual_root_id;

            virtual_root->children.reserve(actual_roots.size());
            for (const auto &actual_root_node : actual_roots)
            {
                if (actual_root_node)
                {
                    actual_root_node->parent = virtual_root;
                    virtual_root->children.push_back(actual_root_node);

                    if (actual_root_node->node_id >= 0 && actual_root_node->node_id < dsu_parent_.size())
                    {
                        dsu_parent_[actual_root_node->node_id] = virtual_root_id;
                    }
                }
            }

            roots_.push_back(virtual_root);
        }

        if (roots_.empty() && next_node_id_ > 0)
        {
#ifdef DEBUG
            std::cerr << "Warning: No root nodes found after build and virtual root creation. Index might be invalid." << std::endl;
#endif
        }

        precomputeAllLCA();

        std::cout << "HypergraphTreeIndex (SizeOnly) build complete." << std::endl;

        if (!cache_path.empty())
        {
            try
            {
                std::filesystem::path p(cache_path);
                if (p.has_parent_path())
                {
                    std::filesystem::create_directories(p.parent_path());
                }
                if (saveIndex(cache_path))
                {
                    std::cout << "HypergraphTreeIndex (SizeOnly) saved to cache: " << cache_path << std::endl;
                }
                else
                {
                    std::cerr << "Warning: Failed to save HypergraphTreeIndex (SizeOnly) cache to '" << cache_path << "'." << std::endl;
                }
            }
            catch (const std::exception &e)
            {
                std::cerr << "Warning: Failed to save HypergraphTreeIndex (SizeOnly) cache '" << cache_path << "'. Error: " << e.what() << std::endl;
            }
        }
    }

    void applyAddVertexToEdges(int vertexId, const std::vector<int> &edgeList)
    {

        hypergraph_.addVertexToEdges(vertexId, edgeList);
#ifdef DEBUG
        std::cout << "Hypergraph updated. Updating Tree index incrementally for vertex " << vertexId << " addition..." << std::endl;
#endif

        bool structure_changed = false;

        std::set<int> unique_edges(edgeList.begin(), edgeList.end());
        std::vector<int> edges_vec(unique_edges.begin(), unique_edges.end());

        for (size_t i = 0; i < edges_vec.size(); ++i)
        {
            for (size_t j = i + 1; j < edges_vec.size(); ++j)
            {
                int edge_u_id = edges_vec[i];
                int edge_v_id = edges_vec[j];
                bool current_pair_structure_changed = false;

#ifdef DEBUG
                std::cout << "Processing pair (" << edge_u_id << ", " << edge_v_id << ")" << std::endl;
#endif

                if (edge_u_id < 0 || edge_u_id >= hyperedge_to_leaf_.size() || !hyperedge_to_leaf_[edge_u_id] ||
                    edge_v_id < 0 || edge_v_id >= hyperedge_to_leaf_.size() || !hyperedge_to_leaf_[edge_v_id])
                {
#ifdef DEBUG
                    std::cerr << "Warning: Invalid edge ID or missing leaf for pair (" << edge_u_id << ", " << edge_v_id << "). Skipping." << std::endl;
#endif
                    continue;
                }
                TreeNodePtr leaf_u = hyperedge_to_leaf_[edge_u_id];
                TreeNodePtr leaf_v = hyperedge_to_leaf_[edge_v_id];

                TreeNodePtr lca_node = find_lca(leaf_u, leaf_v);
                if (!lca_node)
                {
#ifdef DEBUG
                    std::cerr << "Warning: LCA not found for leaves of edges (" << edge_u_id << ", " << edge_v_id << "). Skipping." << std::endl;
#endif
                    continue;
                }
#ifdef DEBUG
                std::cout << "  LCA found: " << lca_node->node_id << " (size: " << lca_node->intersection_size << ")" << std::endl;
#endif

                int new_intersection_size = hypergraph_.getHyperedgeIntersection(edge_u_id, edge_v_id).size();

                if (lca_node->intersection_size > new_intersection_size)
                {
#ifdef DEBUG
                    std::cout << "  Skipping update for pair (" << edge_u_id << ", " << edge_v_id
                              << "): Existing LCA size (" << lca_node->intersection_size
                              << ") > New intersection size (" << new_intersection_size << ")." << std::endl;
#endif
                    continue;
                }

                TreeNodePtr subtree_u = findChildAncestor(leaf_u, lca_node);
                TreeNodePtr subtree_v = findChildAncestor(leaf_v, lca_node);

                if (!subtree_u || !subtree_v)
                {
#ifdef DEBUG
                    std::cerr << "Warning: Could not find child ancestors under LCA " << lca_node->node_id << ". Skipping." << std::endl;
#endif
                    continue;
                }
                if (subtree_u == subtree_v)
                {
#ifdef DEBUG
                    std::cout << "  Leaves share same child ancestor " << subtree_u->node_id << ". Skipping structural change for this pair." << std::endl;
#endif

                    continue;
                }
#ifdef DEBUG
                std::cout << "  Subtree roots under LCA: " << subtree_u->node_id << " and " << subtree_v->node_id << std::endl;
#endif

                int new_parent_id = next_node_id_++;
                TreeNodePtr new_parent = std::make_shared<TreeNode>(new_parent_id);
                new_parent->is_leaf = false;

                new_parent->intersection_size = new_intersection_size;

                if (new_parent_id >= nodes_.size())
                    nodes_.resize(new_parent_id + 1);
                nodes_[new_parent_id] = new_parent;

                if (new_parent_id >= dsu_parent_.size())
                    dsu_parent_.resize(new_parent_id + 1);
                dsu_parent_[new_parent_id] = new_parent_id;

#ifdef DEBUG
                std::cout << "  Creating new parent " << new_parent_id << " with intersection size " << new_intersection_size << std::endl;
#endif

                removeChild(lca_node, subtree_u);
                removeChild(lca_node, subtree_v);
                lca_node->children.push_back(new_parent);
                new_parent->parent = lca_node;

                new_parent->children.push_back(subtree_u);
                new_parent->children.push_back(subtree_v);
                subtree_u->parent = new_parent;
                subtree_v->parent = new_parent;

                current_pair_structure_changed = true;
#ifdef DEBUG
                std::cout << "  Rewired: " << lca_node->node_id << " -> " << new_parent_id << " -> {" << subtree_u->node_id << ", " << subtree_v_id << "}" << std::endl;
#endif

                bool merged_down = false;

                bool can_merge_u = !subtree_u->is_leaf && new_parent->intersection_size == subtree_u->intersection_size;
                bool can_merge_v = !subtree_v->is_leaf && new_parent->intersection_size == subtree_v->intersection_size;

                if (can_merge_u && can_merge_v)
                {

#ifdef DEBUG
                    std::cout << "  Redundancy: Merging new_parent " << new_parent_id << " with BOTH children " << subtree_u->node_id << " and " << subtree_v->node_id << std::endl;
#endif

                    replaceChild(lca_node, new_parent, subtree_u);
                    subtree_u->parent = lca_node;

                    for (const auto &child_of_v : subtree_v->children)
                    {
                        if (child_of_v)
                        {
                            child_of_v->parent = subtree_u;
                            subtree_u->children.push_back(child_of_v);
                        }
                    }

                    subtree_v->children.clear();

                    nodes_[new_parent_id] = nullptr;
                    nodes_[subtree_v->node_id] = nullptr;
                    merged_down = true;
                }
                else if (can_merge_u)
                {

#ifdef DEBUG
                    std::cout << "  Redundancy: Merging new_parent " << new_parent_id << " with child " << subtree_u->node_id << std::endl;
#endif

                    replaceChild(lca_node, new_parent, subtree_u);
                    subtree_u->parent = lca_node;

                    subtree_u->children.push_back(subtree_v);
                    subtree_v->parent = subtree_u;

                    nodes_[new_parent_id] = nullptr;
                    merged_down = true;
                }
                else if (can_merge_v)
                {

#ifdef DEBUG
                    std::cout << "  Redundancy: Merging new_parent " << new_parent_id << " with child " << subtree_v->node_id << std::endl;
#endif

                    replaceChild(lca_node, new_parent, subtree_v);
                    subtree_v->parent = lca_node;

                    subtree_v->children.push_back(subtree_u);
                    subtree_u->parent = subtree_v;

                    nodes_[new_parent_id] = nullptr;
                    merged_down = true;
                }

                if (!merged_down && !lca_node->is_leaf && new_parent->intersection_size == lca_node->intersection_size)
                {
#ifdef DEBUG
                    std::cout << "  Redundancy: Merging new_parent " << new_parent_id << " upwards into parent " << lca_node->node_id << std::endl;
#endif

                    removeChild(lca_node, new_parent);

                    lca_node->children.push_back(subtree_u);
                    subtree_u->parent = lca_node;
                    lca_node->children.push_back(subtree_v);
                    subtree_v->parent = lca_node;

                    nodes_[new_parent_id] = nullptr;
                }
                else if (!merged_down)
                {

#ifdef DEBUG

#endif
                }

                if (current_pair_structure_changed)
                {
#ifdef DEBUG
                    std::cout << "  Structure changed for pair (" << edge_u_id << ", " << edge_v_id << "). Recomputing LCA locally..." << std::endl;
#endif

                    TreeNodePtr subtree_root_for_recompute = lca_node;

                    if (up_.size() < next_node_id_)
                    {
                        size_t old_size = up_.size();
                        up_.resize(next_node_id_);

                        for (size_t k = old_size; k < next_node_id_; ++k)
                        {

                            if (k < up_.size())
                            {
                                up_[k].assign(MAX_LCA_LOG, nullptr);
                            }
                        }
                    }

                    int start_depth = subtree_root_for_recompute->depth + 1;
                    for (const auto &child : subtree_root_for_recompute->children)
                    {
                        if (child)
                        {
#ifdef DEBUG

#endif
                            precompute_lca(child, start_depth);
                        }
                    }

#ifdef DEBUG
                    std::cout << "  LCA locally recomputed starting from node " << subtree_root_for_recompute->node_id << "'s children." << std::endl;
#endif
                }

                if (max_intersection_k_ >= 1)
                {

                    for (int old_size = 1; old_size <= max_intersection_k_; ++old_size)
                    {
                        int list_idx = old_size - 1;

                        if (edge_u_id >= 0 && edge_u_id < edge_neighbors_by_size_.size() &&
                            list_idx >= 0 && list_idx < edge_neighbors_by_size_[edge_u_id].size())
                        {
                            auto &neighbors_u = edge_neighbors_by_size_[edge_u_id][list_idx];
                            neighbors_u.erase(std::remove(neighbors_u.begin(), neighbors_u.end(), edge_v_id), neighbors_u.end());
                        }
                        if (edge_v_id >= 0 && edge_v_id < edge_neighbors_by_size_.size() &&
                            list_idx >= 0 && list_idx < edge_neighbors_by_size_[edge_v_id].size())
                        {
                            auto &neighbors_v = edge_neighbors_by_size_[edge_v_id][list_idx];
                            neighbors_v.erase(std::remove(neighbors_v.begin(), neighbors_v.end(), edge_u_id), neighbors_v.end());
                        }
                    }

                    if (new_intersection_size >= 1 && new_intersection_size <= max_intersection_k_)
                    {
                        int list_idx = new_intersection_size - 1;

                        if (edge_u_id >= 0 && edge_u_id < edge_neighbors_by_size_.size() &&
                            list_idx >= 0 && list_idx < edge_neighbors_by_size_[edge_u_id].size())
                        {
                            auto &neighbors_u = edge_neighbors_by_size_[edge_u_id][list_idx];

                            if (std::find(neighbors_u.begin(), neighbors_u.end(), edge_v_id) == neighbors_u.end())
                            {
                                neighbors_u.push_back(edge_v_id);
                            }
                        }
                        if (edge_v_id >= 0 && edge_v_id < edge_neighbors_by_size_.size() &&
                            list_idx >= 0 && list_idx < edge_neighbors_by_size_[edge_v_id].size())
                        {
                            auto &neighbors_v = edge_neighbors_by_size_[edge_v_id][list_idx];
                            if (std::find(neighbors_v.begin(), neighbors_v.end(), edge_u_id) == neighbors_v.end())
                            {
                                neighbors_v.push_back(edge_u_id);
                            }
                        }
                    }
                }
            }
        }
    }

    void applyRemoveVertexFromEdges(int vertexId, const std::vector<int> &edgeList)
    {

        hypergraph_.removeVertexFromEdges(vertexId, edgeList);
#ifdef DEBUG
        std::cout << "Hypergraph updated. Updating Tree index incrementally for vertex " << vertexId << " removal..." << std::endl;
#endif

        bool global_structure_changed = false;

        std::set<int> unique_edges_set(edgeList.begin(), edgeList.end());
        std::vector<int> unique_edges(unique_edges_set.begin(), unique_edges_set.end());

        if (max_intersection_k_ >= 1)
        {
#ifdef DEBUG
            std::cout << "  Updating edge neighbor information..." << std::endl;
#endif
            for (size_t i = 0; i < unique_edges.size(); ++i)
            {
                for (size_t j = i + 1; j < unique_edges.size(); ++j)
                {
                    int edge_u_id = unique_edges[i];
                    int edge_v_id = unique_edges[j];

                    if (edge_u_id < 0 || edge_u_id >= hypergraph_.numHyperedges() ||
                        edge_v_id < 0 || edge_v_id >= hypergraph_.numHyperedges())
                    {
#ifdef DEBUG
                        std::cerr << "Warning: Invalid edge ID encountered during neighbor update (" << edge_u_id << ", " << edge_v_id << "). Skipping pair." << std::endl;
#endif
                        continue;
                    }

                    int new_intersection_size = hypergraph_.getHyperedgeIntersection(edge_u_id, edge_v_id).size();

                    auto &neighbors_i = edge_neighbors_by_size_[edge_u_id][new_intersection_size + 1];
                    neighbors_i.erase(std::remove(neighbors_i.begin(), neighbors_i.end(), edge_v_id), neighbors_i.end());
                    if (new_intersection_size != 0)
                        edge_neighbors_by_size_[edge_u_id][new_intersection_size].push_back(edge_v_id);
                    auto &neighbors_j = edge_neighbors_by_size_[edge_v_id][new_intersection_size + 1];
                    neighbors_j.erase(std::remove(neighbors_j.begin(), neighbors_j.end(), edge_u_id), neighbors_j.end());
                    if (new_intersection_size != 0)
                        edge_neighbors_by_size_[edge_v_id][new_intersection_size].push_back(edge_u_id);
                }
            }
#ifdef DEBUG
            std::cout << "  Edge neighbor information updated." << std::endl;
#endif
        }

#ifdef DEBUG
        std::cout << "  Checking LCA connectivity for affected edge pairs..." << std::endl;
#endif
        std::vector<TreeNodePtr> affected_lca_nodes;

        for (size_t i = 0; i < unique_edges.size(); ++i)
        {
            for (size_t j = i + 1; j < unique_edges.size(); ++j)
            {
                int edge_u_id = unique_edges[i];
                int edge_v_id = unique_edges[j];

                if (edge_u_id < 0 || edge_u_id >= hyperedge_to_leaf_.size() || !hyperedge_to_leaf_[edge_u_id] ||
                    edge_v_id < 0 || edge_v_id >= hyperedge_to_leaf_.size() || !hyperedge_to_leaf_[edge_v_id])
                {
                    continue;
                }
                TreeNodePtr leaf_u = hyperedge_to_leaf_[edge_u_id];
                TreeNodePtr leaf_v = hyperedge_to_leaf_[edge_v_id];

                TreeNodePtr lca_node = find_lca(leaf_u, leaf_v);
                if (!lca_node || lca_node->is_leaf)
                {

#ifdef DEBUG
                    std::cout << "    Skipping LCA node " << (lca_node ? lca_node->node_id : -1) << " (is leaf or already processed)." << std::endl;
#endif
                    continue;
                }

                if (std::any_of(affected_lca_nodes.begin(), affected_lca_nodes.end(),
                                [&](const TreeNodePtr &node)
                                { return node && node->node_id == lca_node->node_id; }))
                {
#ifdef DEBUG
                    std::cout << "    Skipping LCA node " << lca_node->node_id << " as it is already in the affected list." << std::endl;
#endif
                    continue;
                }

                int original_lca_size = lca_node->intersection_size;

                int new_intersection_size = hypergraph_.getHyperedgeIntersection(edge_u_id, edge_v_id).size();

                if (original_lca_size == new_intersection_size + 1)
                {
#ifdef DEBUG
                    std::cout << "    LCA node " << lca_node->node_id << " (original size: " << original_lca_size << ") potentially affected by pair (" << edge_u_id << "," << edge_v_id << "). Checking connectivity..." << std::endl;
#endif

                    TreeNodePtr subtree_u = findChildAncestor(leaf_u, lca_node);
                    TreeNodePtr subtree_v = findChildAncestor(leaf_v, lca_node);

                    if (!subtree_u || !subtree_v || subtree_u == subtree_v)
                    {

#ifdef DEBUG

#endif
                        continue;
                    }
                    bool connected = checkSubtreeConnectivity(subtree_u, subtree_v, original_lca_size);

                    if (connected)
                    {
#ifdef DEBUG

#endif
                    }
                    else
                    {
#ifdef DEBUG

#endif

                        affected_lca_nodes.push_back(lca_node);
                    }
                }
            }
        }

        std::sort(affected_lca_nodes.begin(), affected_lca_nodes.end(),
                  [](const TreeNodePtr &a, const TreeNodePtr &b)
                  {
                      return a->depth > b->depth;
                  });

        for (auto root_node : affected_lca_nodes)
        {

            std::vector<TreeNodePtr> subtree_roots;
            for (const auto &child : root_node->children)
            {
                subtree_roots.push_back(child);
            }

            int intersection_size = root_node->intersection_size;

            auto new_parent = std::make_shared<TreeNode>(next_node_id_++);
            new_parent->is_leaf = false;
            new_parent->intersection_size = intersection_size - 1;

            int num_subtrees = subtree_roots.size();
            std::map<int, std::vector<TreeNodePtr>> connected_components;

            if (num_subtrees == 0)
            {
#ifdef DEBUG
                std::cout << "      LCA node " << root_node->node_id << " has no children to partition." << std::endl;
#endif
            }
            else if (num_subtrees == 1)
            {

                connected_components[0] = {subtree_roots[0]};
#ifdef DEBUG
                std::cout << "      LCA node " << root_node->node_id << " has only one child, forming a single component." << std::endl;
#endif
            }
            else
            {

                std::vector<int> dsu_parent(num_subtrees);
                std::iota(dsu_parent.begin(), dsu_parent.end(), 0);

                std::function<int(int)> find = [&](int i)
                {
                    if (dsu_parent[i] == i)
                        return i;
                    return dsu_parent[i] = find(dsu_parent[i]);
                };

                auto unite = [&](int i, int j)
                {
                    int root_i = find(i);
                    int root_j = find(j);
                    if (root_i != root_j)
                    {
                        dsu_parent[root_i] = root_j;
                    }
                };

#ifdef DEBUG

#endif

                for (int i = 0; i < num_subtrees; ++i)
                {
                    for (int j = i + 1; j < num_subtrees; ++j)
                    {

                        if (checkSubtreeConnectivity(subtree_roots[i], subtree_roots[j], intersection_size))
                        {

                            unite(i, j);
#ifdef DEBUG
                            std::cout << "        Subtree " << subtree_roots[i]->node_id << " and subtree " << subtree_roots[j]->node_id << " are connected at size " << intersection_size << "." << std::endl;
#endif
                        }
#ifdef DEBUG
                        else
                        {
                            std::cout << "        Subtree " << subtree_roots[i]->node_id << " and subtree " << subtree_roots[j]->node_id << " are NOT connected at size " << intersection_size << "." << std::endl;
                        }
#endif
                    }
                }

#ifdef DEBUG

                int component_idx = 0;
                for (const auto &[component_id, nodes_in_comp] : connected_components)
                {
                    std::cout << "        Component #" << component_idx++ << " (ID: " << component_id << "): ";
                    for (const auto &node_ptr : nodes_in_comp)
                    {
                        std::cout << node_ptr->node_id << " ";
                    }
                    std::cout << std::endl;
                }
#endif
            }

            for (const auto &[component_id, nodes_in_comp] : connected_components)
            {
                if (nodes_in_comp.size() == 1)
                {

                    if (nodes_in_comp[0]->parent.lock())
                    {
                        removeChild(nodes_in_comp[0]->parent.lock(), nodes_in_comp[0]);
                    }

                    new_parent->children.push_back(nodes_in_comp[0]);
                    nodes_in_comp[0]->parent = new_parent;
#ifdef DEBUG

#endif
                }
                else
                {

                    TreeNodePtr new_subtree_parent = std::make_shared<TreeNode>(next_node_id_++);
                    new_subtree_parent->is_leaf = false;
                    new_subtree_parent->intersection_size = intersection_size;

#ifdef DEBUG

#endif

                    for (const auto &node : nodes_in_comp)
                    {
                        if (node->parent.lock())
                        {
                            removeChild(node->parent.lock(), node);
                        }
                        new_subtree_parent->children.push_back(node);
                        node->parent = new_subtree_parent;
#ifdef DEBUG

#endif
                    }

                    new_parent->children.push_back(new_subtree_parent);
                    new_subtree_parent->parent = new_parent;
#ifdef DEBUG

#endif
                }
            }

            TreeNodePtr original_grandparent = root_node->parent.lock();
            bool merged_up = false;

            if (original_grandparent && !original_grandparent->is_leaf &&
                new_parent->intersection_size == original_grandparent->intersection_size)
            {
#ifdef DEBUG
                std::cout << "      Merging new_parent " << new_parent->node_id << " (size " << new_parent->intersection_size << ") upwards into original parent " << original_grandparent->node_id << " (size " << original_grandparent->intersection_size << ")." << std::endl;
#endif
                merged_up = true;

                for (const auto &child : new_parent->children)
                {
                    if (child)
                    {
                        child->parent = original_grandparent;
                        original_grandparent->children.push_back(child);
                    }
                }

                removeChild(original_grandparent, root_node);

                if (new_parent->node_id >= 0 && new_parent->node_id < nodes_.size())
                {
                    nodes_[new_parent->node_id] = nullptr;
                }
                new_parent.reset();
            }
            else
            {

                if (original_grandparent)
                {

                    new_parent->parent = original_grandparent;
                    original_grandparent->children.push_back(new_parent);

                    removeChild(original_grandparent, root_node);
#ifdef DEBUG
                    std::cout << "      Attaching new_parent " << new_parent->node_id << " under original parent " << original_grandparent->node_id << "." << std::endl;
#endif
                }
                else
                {

                    roots_.erase(std::remove(roots_.begin(), roots_.end(), root_node), roots_.end());

                    roots_.push_back(new_parent);
                    new_parent->parent.reset();
#ifdef DEBUG
                    std::cout << "      Setting new_parent " << new_parent->node_id << " as a new global root." << std::endl;
#endif
                }

                if (new_parent)
                {
                    new_parent->depth = original_grandparent ? (original_grandparent->depth + 1) : 0;
                }
            }

            if (root_node->children.empty())
            {

                if (auto parent = root_node->parent.lock())
                {
                    parent->children.erase(std::remove(parent->children.begin(), parent->children.end(), root_node), parent->children.end());
                }

                if (root_node->node_id >= 0 && root_node->node_id < nodes_.size())
                {
                    nodes_[root_node->node_id] = nullptr;
                }
                root_node.reset();
            }
            else
            {

                throw std::runtime_error("Cannot delete root_node as it still has children.");
            }

            TreeNodePtr recompute_start_node = original_grandparent;
            int start_depth = recompute_start_node->depth;

            if (recompute_start_node)
            {
#ifdef DEBUG
                std::cout << "    Recomputing LCA locally starting from node " << recompute_start_node->node_id << " (depth " << start_depth << ")..." << std::endl;
#endif

                if (up_.size() < next_node_id_)
                {
                    size_t old_size = up_.size();
                    up_.resize(next_node_id_);
                    for (size_t k = old_size; k < next_node_id_; ++k)
                    {
                        if (k < up_.size())
                            up_[k].assign(MAX_LCA_LOG, nullptr);
                    }
#ifdef DEBUG

#endif
                }

                precompute_lca(recompute_start_node, start_depth);

                std::cout << "    LCA locally recomputed for subtree rooted at " << recompute_start_node->node_id << "." << std::endl;
                global_structure_changed = true;
            }
            else
            {
#ifdef DEBUG
                std::cout << "    Skipping local LCA recomputation for this iteration (start node not determined)." << std::endl;
#endif

                global_structure_changed = true;
            }
        }
    }

    bool checkSubtreeConnectivity(TreeNodePtr subtree_root_1, TreeNodePtr subtree_root_2, int intersection_size)
    {

        if (!subtree_root_1 || !subtree_root_2 || subtree_root_1 == subtree_root_2)
        {
            return (subtree_root_1 == subtree_root_2);
        }

        if (intersection_size < 1 || intersection_size > max_intersection_k_)
        {
            return false;
        }

        std::vector<int> leaves_1 = getAllLeafEdges(subtree_root_1);
        std::vector<int> leaves_2 = getAllLeafEdges(subtree_root_2);

        if (leaves_1.empty() || leaves_2.empty())
        {
            return false;
        }

        int list_idx = intersection_size;

        const std::vector<int> *leaves_to_iterate;
        std::unordered_set<int> lookup_set;

        if (leaves_1.size() <= leaves_2.size())
        {
            leaves_to_iterate = &leaves_1;
            lookup_set.insert(leaves_2.begin(), leaves_2.end());
        }
        else
        {
            leaves_to_iterate = &leaves_2;
            lookup_set.insert(leaves_1.begin(), leaves_1.end());
        }

        for (int leaf_edge : *leaves_to_iterate)
        {

            if (leaf_edge >= 0 && leaf_edge < edge_neighbors_by_size_.size() &&
                list_idx >= 0 && list_idx < edge_neighbors_by_size_[leaf_edge].size())
            {

                const auto &neighbors = edge_neighbors_by_size_[leaf_edge][list_idx];

                for (int neighbor_edge : neighbors)
                {
                    if (lookup_set.count(neighbor_edge))
                    {
                        return true;
                    }
                }
            }
        }

        return false;
    }

    double getMemoryUsageMB() const
    {
        size_t memory_bytes = 0;
        memory_bytes += sizeof(*this);

        memory_bytes += sizeof(nodes_);
        memory_bytes += nodes_.capacity() * sizeof(TreeNodePtr);
        for (const auto &node_ptr : nodes_)
        {
            if (node_ptr)
            {
                memory_bytes += sizeof(TreeNode);

                memory_bytes += node_ptr->children.capacity() * sizeof(TreeNodePtr);
            }
        }

        memory_bytes += sizeof(roots_);
        memory_bytes += roots_.capacity() * sizeof(TreeNodePtr);

        memory_bytes += sizeof(hyperedge_to_leaf_);
        memory_bytes += hyperedge_to_leaf_.capacity() * sizeof(TreeNodePtr);

        memory_bytes += sizeof(dsu_parent_);
        memory_bytes += dsu_parent_.capacity() * sizeof(int);

        memory_bytes += sizeof(up_);
        memory_bytes += up_.capacity() * sizeof(std::vector<TreeNodePtr>);
        for (const auto &vec : up_)
        {
            memory_bytes += sizeof(vec);
            memory_bytes += vec.capacity() * sizeof(TreeNodePtr);
        }

        return static_cast<double>(memory_bytes) / (1024.0 * 1024.0);
    }

    const std::vector<int> &getNeighborsBySize(int edge_id, int intersection_size) const
    {

        if (intersection_size < 1 || intersection_size > max_intersection_k_ ||
            edge_id < 0 || edge_id >= static_cast<int>(edge_neighbors_by_size_.size()) ||
            (intersection_size - 1) < 0 || (intersection_size - 1) >= static_cast<int>(edge_neighbors_by_size_[edge_id].size()))
        {

            static const std::vector<int> empty_vec;
            return empty_vec;
        }

        return edge_neighbors_by_size_[edge_id][intersection_size - 1];
    }

    size_t getTotalNodes() const
    {

        return next_node_id_;
    }

private:
    Hypergraph &hypergraph_;
    std::vector<TreeNodePtr> nodes_;
    std::vector<TreeNodePtr> roots_;
    std::vector<TreeNodePtr> hyperedge_to_leaf_;
    int next_node_id_;
    std::vector<int> dsu_parent_;
    static const int MAX_LCA_LOG = 20;
    std::vector<std::vector<TreeNodePtr>> up_;

    int max_intersection_k_ = 10;

    std::vector<std::vector<std::vector<int>>> edge_neighbors_by_size_;

    std::pair<int, std::vector<int>> calculate_intersection(TreeNodePtr n1, TreeNodePtr n2)
    {

        if (!n1 || !n2 || !n1->is_leaf || !n2->is_leaf)
        {
            return {0, {}};
        }

        std::vector<int> intersection_verts = hypergraph_.getHyperedgeIntersection(n1->hyperedge_id, n2->hyperedge_id);
        return {static_cast<int>(intersection_verts.size()), intersection_verts};
    }

    int find_set(int node_id)
    {
        if (node_id < 0 || node_id >= dsu_parent_.size())
        {
            return -1;
        }

        if (dsu_parent_[node_id] == node_id)
        {
            return node_id;
        }

        dsu_parent_[node_id] = find_set(dsu_parent_[node_id]);
        return dsu_parent_[node_id];
    }

    void precompute_lca(TreeNodePtr node, int d)
    {
        if (!node || node->node_id < 0 || node->node_id >= up_.size())
            return;

        node->depth = d;

        TreeNodePtr parent_ptr = node->parent.lock();
        up_[node->node_id][0] = parent_ptr;

        for (int i = 1; i < MAX_LCA_LOG; ++i)
        {

            TreeNodePtr ancestor = up_[node->node_id][i - 1];
            if (ancestor && ancestor->node_id >= 0 && ancestor->node_id < up_.size() && i - 1 < up_[ancestor->node_id].size())
            {

                up_[node->node_id][i] = up_[ancestor->node_id][i - 1];
            }
            else
            {

                up_[node->node_id][i] = nullptr;
            }
        }

        for (const auto &child : node->children)
        {
            precompute_lca(child, d + 1);
        }
    }

    TreeNodePtr find_lca(TreeNodePtr u, TreeNodePtr v)
    {
        if (!u || !v || u->node_id < 0 || v->node_id < 0 || u->node_id >= up_.size() || v->node_id >= up_.size())
            return nullptr;

        if (u->depth < v->depth)
            std::swap(u, v);

        for (int i = MAX_LCA_LOG - 1; i >= 0; --i)
        {

            if (i < up_[u->node_id].size())
            {
                TreeNodePtr ancestor = up_[u->node_id][i];

                if (ancestor && ancestor->depth >= v->depth)
                {
                    u = ancestor;
                }
            }
        }

        if (u->node_id == v->node_id)
            return u;

        for (int i = MAX_LCA_LOG - 1; i >= 0; --i)
        {

            if (i < up_[u->node_id].size() && i < up_[v->node_id].size())
            {
                TreeNodePtr parent_u = up_[u->node_id][i];
                TreeNodePtr parent_v = up_[v->node_id][i];

                if (parent_u && parent_v && parent_u->node_id != parent_v->node_id)
                {
                    u = parent_u;
                    v = parent_v;
                }
            }
        }

        if (!up_[u->node_id].empty())
        {
            return up_[u->node_id][0];
        }
        return nullptr;
    }

    bool saveIndex(const std::string &filename) const
    {
        std::ofstream file(filename);
        if (!file.is_open())
        {
            std::cerr << "Error: Cannot open file for writing index: " << filename << std::endl;
            return false;
        }

        file << "num_nodes " << next_node_id_ << "\n";
        file << "num_hyperedges " << hypergraph_.numHyperedges() << "\n";

        for (int i = 0; i < next_node_id_; ++i)
        {
            const auto &node = nodes_[i];
            if (!node)
                continue;

            file << "node " << node->node_id << " "
                 << node->is_leaf << " "
                 << node->hyperedge_id << " "
                 << node->intersection_size << " ";

            auto parent_ptr = node->parent.lock();
            file << (parent_ptr ? parent_ptr->node_id : -1) << " ";

            file << node->children.size();
            for (const auto &child : node->children)
            {
                if (child)
                    file << " " << child->node_id;
            }
            file << "\n";
        }

        file << "roots " << roots_.size();
        for (const auto &root : roots_)
        {
            if (root)
                file << " " << root->node_id;
        }
        file << "\n";

        file << "hyperedge_to_leaf " << hyperedge_to_leaf_.size();
        for (const auto &leaf_ptr : hyperedge_to_leaf_)
        {
            file << " " << (leaf_ptr ? leaf_ptr->node_id : -1);
        }
        file << "\n";

        return file.good();
    }

    void findRoots()
    {
        roots_.clear();
        std::vector<bool> is_root_added(next_node_id_, false);
        for (int i = 0; i < next_node_id_; ++i)
        {

            if (i < nodes_.size() && nodes_[i] && nodes_[i]->parent.expired())
            {
                int root_id = find_set(i);
                if (root_id != -1 && root_id < is_root_added.size() && !is_root_added[root_id])
                {
                    if (root_id < nodes_.size() && nodes_[root_id])
                    {
                        roots_.push_back(nodes_[root_id]);
                        is_root_added[root_id] = true;
                    }
                }
            }
        }
        if (roots_.empty() && next_node_id_ > 0)
        {
#ifdef DEBUG
            std::cerr << "Warning: No root nodes found after build. Index might be invalid." << std::endl;
#endif
        }
    }

    void clearIndexData()
    {
        nodes_.clear();
        roots_.clear();
        next_node_id_ = 0;

        size_t num_hyperedges = hypergraph_.numHyperedges();
        if (hyperedge_to_leaf_.size() != num_hyperedges)
        {
            hyperedge_to_leaf_.resize(num_hyperedges, nullptr);
        }
        else
        {
            std::fill(hyperedge_to_leaf_.begin(), hyperedge_to_leaf_.end(), nullptr);
        }
        dsu_parent_.clear();
        up_.clear();

        if (max_intersection_k_ >= 1)
        {
            size_t num_hyperedges = hypergraph_.numHyperedges();
            edge_neighbors_by_size_.assign(num_hyperedges, std::vector<std::vector<int>>(max_intersection_k_));

            for (auto &mid_vec : edge_neighbors_by_size_)
            {
                for (auto &inner_vec : mid_vec)
                {
                    inner_vec.clear();
                }
            }
        }
        else
        {
            edge_neighbors_by_size_.clear();
        }
    }

    void mergeNodesSizeOnly(const std::vector<std::tuple<int, int, int>> &merge_candidates)
    {
        for (const auto &candidate : merge_candidates)
        {
            int size = std::get<0>(candidate);
            int node1_id = std::get<1>(candidate);
            int node2_id = std::get<2>(candidate);

            int root1_id = find_set(node1_id);
            int root2_id = find_set(node2_id);

            if (root1_id != root2_id && root1_id != -1 && root2_id != -1)
            {
                TreeNodePtr root1_node = nodes_[root1_id];
                TreeNodePtr root2_node = nodes_[root2_id];
                bool merged = false;

                auto updateDsu = [&](const auto &self, TreeNodePtr node, int new_root_id) -> void
                {
                    if (node->node_id < 0 || node->node_id >= dsu_parent_.size())
                        return;
                    dsu_parent_[node->node_id] = new_root_id;
                    for (auto &child : node->children)
                    {
                        self(self, child, new_root_id);
                    }
                };

                if (!root1_node->is_leaf && root1_node->intersection_size == size)
                {
                    root1_node->children.push_back(root2_node);
                    root2_node->parent = root1_node;
                    updateDsu(updateDsu, root2_node, root1_id);

                    merged = true;
                }
                else if (!root2_node->is_leaf && root2_node->intersection_size == size)
                {
                    root2_node->children.push_back(root1_node);
                    root1_node->parent = root2_node;
                    updateDsu(updateDsu, root1_node, root2_id);

                    merged = true;
                }

                if (!merged)
                {
                    int parent_node_id = next_node_id_++;
                    TreeNodePtr parent = std::make_shared<TreeNode>(parent_node_id);
                    parent->is_leaf = false;
                    parent->intersection_size = size;

                    parent->children.push_back(root1_node);
                    parent->children.push_back(root2_node);
                    root1_node->parent = parent;
                    root2_node->parent = parent;

                    if (parent_node_id >= nodes_.size())
                        nodes_.resize(parent_node_id + 1);
                    nodes_[parent_node_id] = parent;
                    if (parent_node_id >= dsu_parent_.size())
                        dsu_parent_.resize(parent_node_id + 1);

                    dsu_parent_[parent_node_id] = parent_node_id;
                    dsu_parent_[root1_id] = parent_node_id;
                    dsu_parent_[root2_id] = parent_node_id;
                }
            }
        }
    }

    void mergeNodes(const std::vector<std::tuple<int, std::vector<int>, int, int>> &merge_candidates)
    {
        for (const auto &candidate : merge_candidates)
        {
            int size = std::get<0>(candidate);
            const auto &vertices = std::get<1>(candidate);
            int node1_id = std::get<2>(candidate);
            int node2_id = std::get<3>(candidate);

            int root1_id = find_set(node1_id);
            int root2_id = find_set(node2_id);

            if (root1_id != root2_id && root1_id != -1 && root2_id != -1)
            {
                TreeNodePtr root1_node = nodes_[root1_id];
                TreeNodePtr root2_node = nodes_[root2_id];
                bool merged = false;

                auto updateDsu = [&](const auto &self, TreeNodePtr node, int new_root_id) -> void
                {
                    if (node->node_id < 0 || node->node_id >= dsu_parent_.size())
                        return;
                    dsu_parent_[node->node_id] = new_root_id;
                    for (auto &child : node->children)
                    {
                        if (child)
                            self(self, child, new_root_id);
                    }
                };

                if (!root1_node->is_leaf && root1_node->intersection_size == size && root1_node->intersection_vertices == vertices)
                {
                    root1_node->children.push_back(root2_node);
                    root2_node->parent = root1_node;
                    updateDsu(updateDsu, root2_node, root1_id);

                    merged = true;
                }
                else if (!root2_node->is_leaf && root2_node->intersection_size == size && root2_node->intersection_vertices == vertices)
                {
                    root2_node->children.push_back(root1_node);
                    root1_node->parent = root2_node;
                    updateDsu(updateDsu, root1_node, root2_id);

                    merged = true;
                }

                if (!merged)
                {
                    int parent_node_id = next_node_id_++;
                    TreeNodePtr parent = std::make_shared<TreeNode>(parent_node_id);
                    parent->is_leaf = false;
                    parent->intersection_size = size;
                    parent->intersection_vertices = vertices;
                    parent->children.push_back(root1_node);
                    parent->children.push_back(root2_node);
                    root1_node->parent = parent;
                    root2_node->parent = parent;

                    if (parent_node_id >= nodes_.size())
                        nodes_.resize(parent_node_id + 1);
                    nodes_[parent_node_id] = parent;
                    if (parent_node_id >= dsu_parent_.size())
                        dsu_parent_.resize(parent_node_id + 1);

                    dsu_parent_[parent_node_id] = parent_node_id;
                    dsu_parent_[root1_id] = parent_node_id;
                    dsu_parent_[root2_id] = parent_node_id;
                }
            }
        }
    }

    void sortMergeCandidatesSizeOnly(std::vector<std::tuple<int, int, int>> &merge_candidates)
    {
        std::sort(merge_candidates.begin(), merge_candidates.end(),
                  [](const auto &a, const auto &b)
                  {
                      return std::get<0>(a) > std::get<0>(b);
                  });
    }

    void calculateAndPrepareMergeCandidates(
        const std::vector<TreeNodePtr> &current_leaf_nodes,
        std::vector<std::tuple<int, std::vector<int>, int, int>> &merge_candidates)
    {
        std::mutex merge_candidates_mutex;
        size_t num_leaves = current_leaf_nodes.size();
        unsigned int num_threads = std::thread::hardware_concurrency();
        if (num_threads == 0)
            num_threads = 1;
        std::vector<std::thread> threads(num_threads);

        merge_candidates.reserve(num_leaves * 1000);

        auto calculate_leaf_intersections =
            [&](size_t start_idx, size_t end_idx)
        {
            std::vector<std::tuple<int, std::vector<int>, int, int>> local_candidates;
            local_candidates.reserve((end_idx - start_idx) * num_leaves / 2);

            for (size_t i = start_idx; i < end_idx; ++i)
            {
                for (size_t j = i + 1; j < num_leaves; ++j)
                {
                    auto intersection_result = calculate_intersection(current_leaf_nodes[i], current_leaf_nodes[j]);
                    int intersection_size = intersection_result.first;
                    std::vector<int> intersection_verts = std::move(intersection_result.second);

                    if (intersection_size > 0)
                    {
                        local_candidates.emplace_back(intersection_size, std::move(intersection_verts), current_leaf_nodes[i]->node_id, current_leaf_nodes[j]->node_id);
                    }
                }
            }
            std::lock_guard<std::mutex> lock(merge_candidates_mutex);
            merge_candidates.insert(merge_candidates.end(),
                                    std::make_move_iterator(local_candidates.begin()),
                                    std::make_move_iterator(local_candidates.end()));
        };

        size_t chunk_size = (num_leaves + num_threads - 1) / num_threads;
        size_t current_start = 0;
        for (unsigned int t = 0; t < num_threads; ++t)
        {
            size_t current_end = std::min(current_start + chunk_size, num_leaves);
            if (current_start >= current_end)
                break;
            threads[t] = std::thread(calculate_leaf_intersections, current_start, current_end);
            current_start = current_end;
        }

        for (unsigned int t = 0; t < threads.size(); ++t)
        {
            if (threads[t].joinable())
            {
                threads[t].join();
            }
        }
    }

    void calculateAndPrepareMergeCandidatesSizeOnly(
        const std::vector<TreeNodePtr> &current_leaf_nodes,
        std::vector<std::tuple<int, int, int>> &merge_candidates)
    {
        std::mutex merge_candidates_mutex;
        size_t num_leaves = current_leaf_nodes.size();
        unsigned int num_threads = std::thread::hardware_concurrency();
        if (num_threads == 0)
            num_threads = 1;
        std::vector<std::thread> threads(num_threads);

        merge_candidates.reserve(num_leaves * 100);

        auto calculate_leaf_intersections_size_only =
            [&](size_t start_idx, size_t end_idx)
        {
            std::vector<std::tuple<int, int, int>> local_candidates;
            local_candidates.reserve(10000);

            for (size_t i = start_idx; i < end_idx; ++i)
            {
                for (size_t j = i + 1; j < num_leaves; ++j)
                {
                    auto intersection_result = calculate_intersection(current_leaf_nodes[i], current_leaf_nodes[j]);
                    int intersection_size = intersection_result.first;

                    if (intersection_size > 0)
                    {
                        local_candidates.emplace_back(intersection_size, current_leaf_nodes[i]->node_id, current_leaf_nodes[j]->node_id);
                    }
                }
            }
            std::lock_guard<std::mutex> lock(merge_candidates_mutex);
            merge_candidates.insert(merge_candidates.end(), local_candidates.begin(), local_candidates.end());
        };

        size_t chunk_size = (num_leaves + num_threads - 1) / num_threads;
        size_t current_start = 0;
        for (unsigned int t = 0; t < num_threads; ++t)
        {
            size_t current_end = std::min(current_start + chunk_size, num_leaves);
            if (current_start >= current_end)
                break;
            threads[t] = std::thread(calculate_leaf_intersections_size_only, current_start, current_end);
            current_start = current_end;
        }

        for (unsigned int t = 0; t < threads.size(); ++t)
        {
            if (threads[t].joinable())
            {
                threads[t].join();
            }
        }
    }

    void sortMergeCandidates(std::vector<std::tuple<int, std::vector<int>, int, int>> &merge_candidates)
    {
        std::sort(merge_candidates.begin(), merge_candidates.end(),
                  [](const auto &a, const auto &b)
                  {
                      if (std::get<0>(a) != std::get<0>(b))
                      {
                          return std::get<0>(a) > std::get<0>(b);
                      }

                      return std::get<1>(a) < std::get<1>(b);
                  });
    }

    void initializeDSU()
    {
        dsu_parent_.resize(next_node_id_);
        for (int i = 0; i < next_node_id_; ++i)
        {
            dsu_parent_[i] = i;
        }
    }

    void recomputeAllLCA()
    {
        if (next_node_id_ == 0)
            return;

        if (up_.size() < next_node_id_)
        {
            up_.resize(next_node_id_);
        }

        for (int i = 0; i < next_node_id_; ++i)
        {
            if (i < up_.size())
            {
                up_[i].assign(MAX_LCA_LOG, nullptr);
            }

            if (i < nodes_.size() && nodes_[i])
            {
                nodes_[i]->depth = 0;
            }
        }

        for (const auto &root : roots_)
        {
            if (root)
            {
                precompute_lca(root, 0);
            }
        }
    }

    bool loadIndex(const std::string &filename)
    {
        std::ifstream file(filename);
        if (!file.is_open())
        {
            return false;
        }

        std::string line_type;
        int num_nodes = 0;
        size_t num_hyperedges_expected = 0;

        file >> line_type >> num_nodes;
        if (file.fail() || line_type != "num_nodes" || num_nodes <= 0)
        {
            std::cerr << "Error: Invalid cache file format (num_nodes)." << std::endl;
            return false;
        }
        file >> line_type >> num_hyperedges_expected;
        if (file.fail() || line_type != "num_hyperedges" || num_hyperedges_expected != hypergraph_.numHyperedges())
        {
            std::cerr << "Error: Cache file hyperedge count mismatch with current hypergraph." << std::endl;
            return false;
        }

        clearIndexData();
        nodes_.resize(num_nodes, nullptr);
        hyperedge_to_leaf_.resize(num_hyperedges_expected, nullptr);
        next_node_id_ = num_nodes;

        std::vector<int> parent_ids(num_nodes, -1);
        std::vector<std::vector<int>> children_ids(num_nodes);

        std::string line;
        std::getline(file, line);

        while (std::getline(file, line))
        {
            std::istringstream iss(line);
            iss >> line_type;

            if (line_type == "node")
            {
                int id, h_id, i_size, p_id, child_count;
                bool is_leaf;

                iss >> id >> is_leaf >> h_id >> i_size >> p_id >> child_count;
                if (iss.fail() || id < 0 || id >= num_nodes)
                {
                    std::cerr << "Error: Invalid node data in cache." << std::endl;
                    return false;
                }

                nodes_[id] = std::make_shared<TreeNode>(id);
                nodes_[id]->is_leaf = is_leaf;
                nodes_[id]->hyperedge_id = h_id;
                nodes_[id]->intersection_size = i_size;

                auto parent_ptr = nodes_[p_id];
                if (p_id != -1 && parent_ptr)
                {
                    nodes_[id]->parent = parent_ptr;
                    parent_ptr->children.push_back(nodes_[id]);
                }

                children_ids[id].resize(child_count);
                for (int i = 0; i < child_count; ++i)
                {
                    if (!(iss >> children_ids[id][i]))
                    {
                        std::cerr << "Error: Reading child IDs failed in cache." << std::endl;
                        return false;
                    }
                }
            }
            else if (line_type == "roots")
            {
                int count, root_id;
                iss >> count;
                if (iss.fail())
                {
                    std::cerr << "Error: Reading root count failed." << std::endl;
                    return false;
                }
                roots_.reserve(count);
                for (int i = 0; i < count; ++i)
                {
                    iss >> root_id;
                    if (iss.fail() || root_id < 0 || root_id >= num_nodes || !nodes_[root_id])
                    {
                        std::cerr << "Error: Invalid root ID in cache." << std::endl;
                        return false;
                    }
                    roots_.push_back(nodes_[root_id]);
                }
            }
            else if (line_type == "hyperedge_to_leaf")
            {
                size_t count;
                int leaf_id;
                iss >> count;
                if (iss.fail() || count != hyperedge_to_leaf_.size())
                {
                    std::cerr << "Error: hyperedge_to_leaf count mismatch in cache." << std::endl;
                    return false;
                }
                for (size_t i = 0; i < count; ++i)
                {
                    iss >> leaf_id;
                    if (iss.fail())
                    {
                        std::cerr << "Error: Reading leaf ID failed." << std::endl;
                        return false;
                    }
                    if (leaf_id >= 0 && leaf_id < num_nodes && nodes_[leaf_id])
                    {
                        hyperedge_to_leaf_[i] = nodes_[leaf_id];
                    }
                    else if (leaf_id != -1)
                    {
                        std::cerr << "Error: Invalid leaf ID in hyperedge_to_leaf map." << std::endl;
                        return false;
                    }
                }
            }
        }

        if (file.bad())
        {
            std::cerr << "Error: File read error occurred during cache loading." << std::endl;
            return false;
        }

        for (int i = 0; i < num_nodes; ++i)
        {
            if (!nodes_[i])
                continue;

            if (parent_ids[i] != -1)
            {
                if (parent_ids[i] >= 0 && parent_ids[i] < num_nodes && nodes_[parent_ids[i]])
                {
                    nodes_[i]->parent = nodes_[parent_ids[i]];
                }
                else
                {
                    std::cerr << "Error: Invalid parent link ID (" << parent_ids[i] << ") for node " << i << " in cache." << std::endl;
                    return false;
                }
            }

            nodes_[i]->children.reserve(children_ids[i].size());
            for (int child_id : children_ids[i])
            {
                if (child_id >= 0 && child_id < num_nodes && nodes_[child_id])
                {
                    nodes_[i]->children.push_back(nodes_[child_id]);
                }
                else
                {
                    std::cerr << "Error: Invalid child link ID (" << child_id << ") for node " << i << " in cache." << std::endl;
                    return false;
                }
            }
        }

        if (roots_.empty() && num_nodes > 0)
        {
            std::cerr << "Error: No roots loaded from cache, but nodes exist." << std::endl;

            findRoots();
            if (roots_.empty())
                return false;
        }

        recomputeAllLCA();

        return true;
    }

    void precomputeAllLCA()
    {
        if (next_node_id_ == 0)
            return;

        if (up_.size() < next_node_id_)
        {
            up_.resize(next_node_id_);
        }

        for (int i = 0; i < next_node_id_; ++i)
        {
            if (i < up_.size())
            {
                up_[i].assign(MAX_LCA_LOG, nullptr);
            }

            if (i < nodes_.size() && nodes_[i])
            {
                nodes_[i]->depth = 0;
            }
        }

        for (const auto &root : roots_)
        {
            if (root)
            {
                precompute_lca(root, 0);
            }
        }
    }

    void populateEdgeNeighbors(const std::vector<std::tuple<int, std::vector<int>, int, int>> &merge_candidates)
    {
        if (max_intersection_k_ < 1)
            return;

        for (auto &mid_vec : edge_neighbors_by_size_)
        {
            for (auto &inner_vec : mid_vec)
            {
                inner_vec.clear();
            }
        }

        for (const auto &candidate : merge_candidates)
        {
            int size = std::get<0>(candidate);
            int node1_id = std::get<2>(candidate);
            int node2_id = std::get<3>(candidate);

            if (size >= 1 && size <= max_intersection_k_ &&
                node1_id >= 0 && node1_id < nodes_.size() && nodes_[node1_id]->is_leaf &&
                node2_id >= 0 && node2_id < nodes_.size() && nodes_[node2_id]->is_leaf)
            {
                int edge1_id = nodes_[node1_id]->hyperedge_id;
                int edge2_id = nodes_[node2_id]->hyperedge_id;

                if (edge1_id >= 0 && edge1_id < edge_neighbors_by_size_.size() &&
                    edge2_id >= 0 && edge2_id < edge_neighbors_by_size_.size() &&
                    (size - 1) >= 0 && (size - 1) < edge_neighbors_by_size_[edge1_id].size() &&
                    (size - 1) < edge_neighbors_by_size_[edge2_id].size())
                {
                    edge_neighbors_by_size_[edge1_id][size - 1].push_back(edge2_id);
                    edge_neighbors_by_size_[edge2_id][size - 1].push_back(edge1_id);
                }
            }
        }
    }

    void populateEdgeNeighborsSizeOnly(const std::vector<std::tuple<int, int, int>> &merge_candidates_size_only)
    {
        if (max_intersection_k_ < 1)
            return;

        for (auto &mid_vec : edge_neighbors_by_size_)
        {
            for (auto &inner_vec : mid_vec)
            {
                inner_vec.clear();
            }
        }

        for (const auto &candidate : merge_candidates_size_only)
        {
            int size = std::get<0>(candidate);
            int node1_id = std::get<1>(candidate);
            int node2_id = std::get<2>(candidate);

            if (size >= 1 && size <= max_intersection_k_ &&
                node1_id >= 0 && node1_id < nodes_.size() && nodes_[node1_id]->is_leaf &&
                node2_id >= 0 && node2_id < nodes_.size() && nodes_[node2_id]->is_leaf)
            {
                int edge1_id = nodes_[node1_id]->hyperedge_id;
                int edge2_id = nodes_[node2_id]->hyperedge_id;

                if (edge1_id >= 0 && edge1_id < edge_neighbors_by_size_.size() &&
                    edge2_id >= 0 && edge2_id < edge_neighbors_by_size_.size() &&
                    (size) > 0 && (size) < edge_neighbors_by_size_[edge1_id].size() &&
                    (size) < edge_neighbors_by_size_[edge2_id].size())
                {
                    edge_neighbors_by_size_[edge1_id][size].push_back(edge2_id);
                    edge_neighbors_by_size_[edge2_id][size].push_back(edge1_id);
                }
            }
        }
    }

    void recalculateAndPopulateNeighbors()
    {
        if (max_intersection_k_ < 2)
            return;

        std::vector<TreeNodePtr> current_leaf_nodes;
        current_leaf_nodes.reserve(hyperedge_to_leaf_.size());
        for (const auto &leaf_ptr : hyperedge_to_leaf_)
        {
            if (leaf_ptr)
            {
                current_leaf_nodes.push_back(leaf_ptr);
            }
        }
        if (current_leaf_nodes.empty())
            return;

        std::vector<std::tuple<int, std::vector<int>, int, int>> merge_candidates;
        calculateAndPrepareMergeCandidates(current_leaf_nodes, merge_candidates);

        populateEdgeNeighbors(merge_candidates);
    }

    void recalculateAndPopulateNeighborsSizeOnly()
    {
        if (max_intersection_k_ < 2)
            return;

        std::vector<TreeNodePtr> current_leaf_nodes;
        current_leaf_nodes.reserve(hyperedge_to_leaf_.size());
        for (const auto &leaf_ptr : hyperedge_to_leaf_)
        {
            if (leaf_ptr)
            {
                current_leaf_nodes.push_back(leaf_ptr);
            }
        }
        if (current_leaf_nodes.empty())
            return;

        std::vector<std::tuple<int, int, int>> merge_candidates_size_only;
        calculateAndPrepareMergeCandidatesSizeOnly(current_leaf_nodes, merge_candidates_size_only);

        populateEdgeNeighborsSizeOnly(merge_candidates_size_only);
    }

    void buildIndexCache(std::string cache_path = "")
    {

        bool loaded_from_cache = false;
        if (!cache_path.empty())
        {
            cache_path = cache_path + "hypergraph_tree_index";
            try
            {
                loaded_from_cache = loadIndex(cache_path);
                if (loaded_from_cache)
                {
                    std::cout << "HypergraphTreeIndex loaded from cache: " << cache_path << std::endl;
#ifdef DEBUG
                    std::cout << "Recalculating neighbor information after loading from cache..." << std::endl;
#endif
                    recalculateAndPopulateNeighbors();
#ifdef DEBUG
                    std::cout << "Neighbor information recalculated." << std::endl;
#endif

                    return;
                }
                else
                {

                    std::cout << "Info: HypergraphTreeIndex cache not found or invalid at '" << cache_path << "'. Building index." << std::endl;
                }
            }
            catch (const std::exception &e)
            {
                std::cerr << "Warning: Failed to load HypergraphTreeIndex cache '" << cache_path << "'. Building index. Error: " << e.what() << std::endl;
            }
        }

        std::cout << "Building HypergraphTreeIndex..." << std::endl;
        int num_hyperedges = hypergraph_.numHyperedges();
        if (num_hyperedges == 0)
        {
            std::cout << "Warning: Hypergraph has no hyperedges. Index is empty." << std::endl;
            return;
        }

        clearIndexData();

        std::vector<TreeNodePtr> current_leaf_nodes;
        current_leaf_nodes.reserve(num_hyperedges);
        for (int i = 0; i < num_hyperedges; ++i)
        {
            if (hypergraph_.getHyperedge(i).size() == 0)
                continue;

            TreeNodePtr leaf = std::make_shared<TreeNode>(next_node_id_++);
            leaf->is_leaf = true;
            leaf->hyperedge_id = i;
            nodes_.push_back(leaf);

            if (i < hyperedge_to_leaf_.size())
            {
                hyperedge_to_leaf_[i] = leaf;
            }
            else
            {
                throw std::out_of_range("Hyperedge ID out of range for hyperedge_to_leaf_ vector");
            }
            current_leaf_nodes.push_back(leaf);
        }

        if (current_leaf_nodes.empty())
        {
            std::cout << "Warning: No valid hyperedges found. Index is empty." << std::endl;
            return;
        }

        std::vector<std::tuple<int, std::vector<int>, int, int>> merge_candidates;
        calculateAndPrepareMergeCandidates(current_leaf_nodes, merge_candidates);

        sortMergeCandidates(merge_candidates);

        populateEdgeNeighbors(merge_candidates);

        dsu_parent_.resize(next_node_id_);
        for (int i = 0; i < next_node_id_; ++i)
        {
            dsu_parent_[i] = i;
        }

        for (const auto &candidate : merge_candidates)
        {
            int size = std::get<0>(candidate);
            const auto &vertices = std::get<1>(candidate);
            int node1_id = std::get<2>(candidate);
            int node2_id = std::get<3>(candidate);

            int root1_id = find_set(node1_id);
            int root2_id = find_set(node2_id);

            if (root1_id != root2_id && root1_id != -1 && root2_id != -1)
            {
                TreeNodePtr root1_node = nodes_[root1_id];
                TreeNodePtr root2_node = nodes_[root2_id];

                bool merged = false;

                auto updateDsu = [&](const auto &self, TreeNodePtr node, int new_root_id) -> void
                {
                    dsu_parent_[node->node_id] = new_root_id;
                    for (auto &child : node->children)
                    {
                        self(self, child, new_root_id);
                    }
                };

                if (!root1_node->is_leaf && root1_node->intersection_size == size)
                {

                    root1_node->children.push_back(root2_node);
                    root2_node->parent = root1_node;
                    updateDsu(updateDsu, root2_node, root1_id);
                    dsu_parent_[root2_id] = root1_id;
                    merged = true;
                }
                else if (!root2_node->is_leaf && root2_node->intersection_size == size)
                {

                    root2_node->children.push_back(root1_node);
                    root1_node->parent = root2_node;
                    updateDsu(updateDsu, root1_node, root2_id);
                    dsu_parent_[root1_id] = root2_id;
                    merged = true;
                }

                if (!merged)
                {
                    int parent_node_id = next_node_id_++;
                    TreeNodePtr parent = std::make_shared<TreeNode>(parent_node_id);
                    parent->is_leaf = false;
                    parent->intersection_size = size;

                    parent->children.push_back(root1_node);
                    parent->children.push_back(root2_node);
                    root1_node->parent = parent;
                    root2_node->parent = parent;

                    if (parent_node_id >= nodes_.size())
                        nodes_.resize(parent_node_id + 1);
                    nodes_[parent_node_id] = parent;

                    if (parent_node_id >= dsu_parent_.size())
                        dsu_parent_.resize(parent_node_id + 1);

                    dsu_parent_[parent_node_id] = parent_node_id;
                    dsu_parent_[root1_id] = parent_node_id;
                    dsu_parent_[root2_id] = parent_node_id;
                }
            }
        }

        roots_.clear();
        std::vector<bool> is_root_added(next_node_id_, false);
        for (int i = 0; i < next_node_id_; ++i)
        {
            if (i < nodes_.size() && nodes_[i] && nodes_[i]->parent.expired())
            {
                int root_id = find_set(i);
                if (root_id != -1 && root_id < is_root_added.size() && !is_root_added[root_id])
                {
                    if (root_id < nodes_.size() && nodes_[root_id])
                    {
                        roots_.push_back(nodes_[root_id]);
                        is_root_added[root_id] = true;
                    }
                }
            }
        }

        if (up_.size() < next_node_id_)
        {
            up_.resize(next_node_id_);
        }
        for (int i = 0; i < next_node_id_; ++i)
        {
            up_[i].assign(MAX_LCA_LOG, nullptr);
        }
        for (const auto &root : roots_)
        {
            precompute_lca(root, 0);
        }
    }

    TreeNodePtr findChildAncestor(TreeNodePtr node, TreeNodePtr ancestor)
    {
        if (!node || !ancestor || node == ancestor)
        {
            return nullptr;
        }
        TreeNodePtr current = node;
        while (current)
        {
            TreeNodePtr parent = current->parent.lock();
            if (parent == ancestor)
            {
                return current;
            }
            if (parent == nullptr || parent == node || parent == current)
            {
                break;
            }
            current = parent;
        }

        return nullptr;
    }

    void replaceChild(TreeNodePtr parent, TreeNodePtr old_child, TreeNodePtr new_child)
    {
        if (!parent || !old_child || !new_child)
            return;
        auto &children = parent->children;
        bool found = false;
        for (auto &child : children)
        {

            if (child && child->node_id == old_child->node_id)
            {
                child = new_child;
                found = true;
                break;
            }
        }
        if (!found)
        {
#ifdef DEBUG
            std::cerr << "Warning: Tried to replace child " << old_child->node_id << " under parent " << parent->node_id << ", but it was not found." << std::endl;
#endif
        }
    }

    void removeChild(TreeNodePtr parent, TreeNodePtr child_to_remove)
    {
        if (!parent || !child_to_remove)
            return;
        auto &children = parent->children;

        children.erase(std::remove_if(children.begin(), children.end(),
                                      [&](const TreeNodePtr &child)
                                      {
                                          return child && child->node_id == child_to_remove->node_id;
                                      }),
                       children.end());
    }

    void getAllLeafEdgesRecursive(TreeNodePtr node, std::vector<int> &leaf_edges)
    {
        if (!node)
            return;
        if (node->is_leaf)
        {
            if (node->hyperedge_id != -1)
            {
                leaf_edges.push_back(node->hyperedge_id);
            }
        }
        else
        {
            for (const auto &child : node->children)
            {
                getAllLeafEdgesRecursive(child, leaf_edges);
            }
        }
    }

    std::vector<int> getAllLeafEdges(TreeNodePtr subtree_root)
    {
        std::vector<int> leaf_edges;
        getAllLeafEdgesRecursive(subtree_root, leaf_edges);

        return leaf_edges;
    }

    bool checkConnectivity(const std::vector<int> &leaves1, const std::vector<int> &leaves2, int required_size)
    {
        if (required_size < 1 || required_size > max_intersection_k_)
        {
            return false;
        }

        int list_idx = required_size - 1;

        std::unordered_set<int> leaves2_set(leaves2.begin(), leaves2.end());
        if (leaves2_set.empty() || leaves1.empty())
            return false;

        for (int edge1 : leaves1)
        {

            if (edge1 < 0 || edge1 >= edge_neighbors_by_size_.size() ||
                list_idx < 0 || list_idx >= edge_neighbors_by_size_[edge1].size())
            {
                continue;
            }

            const auto &neighbors = edge_neighbors_by_size_[edge1][list_idx];
            for (int neighbor_edge : neighbors)
            {
                if (leaves2_set.count(neighbor_edge))
                {
                    return true;
                }
            }
        }
        return false;
    }
};

#endif
