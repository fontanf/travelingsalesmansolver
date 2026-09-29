#include "travelingsalesmansolver/candidates/alpha_nearness.hpp"

#include "travelingsalesmansolver/distances/ball_tree.hpp"

#include <algorithm>
#include <limits>
#include <numeric>
#include <queue>
#include <type_traits>

using namespace travelingsalesmansolver;

namespace
{

template <typename T, typename = void>
struct HasCoordinates: std::false_type { };

template <typename T>
struct HasCoordinates<T, std::void_t<decltype(std::declval<const T&>().coordinates(VertexId(0)))>>: std::true_type { };

constexpr Distance infinity = std::numeric_limits<Distance>::max();

/**
 * 2D k-d tree, to find the nearest neighbors of a point in each of the 8
 * octants around it.
 */
class KdTree
{

public:

    KdTree(const std::vector<Coordinates2D>& points):
        points_(points),
        ids_(points.size())
    {
        std::iota(ids_.begin(), ids_.end(), 0);
        nodes_.push_back({0, (VertexId)points.size()});
        for (NodeId node_id = 0; node_id < (NodeId)nodes_.size(); ++node_id) {
            Node node = nodes_[node_id];
            node.x_min = node.y_min = std::numeric_limits<double>::infinity();
            node.x_max = node.y_max = -std::numeric_limits<double>::infinity();
            for (VertexId pos = node.begin; pos < node.end; ++pos) {
                const Coordinates2D& point = points_[ids_[pos]];
                node.x_min = std::min(node.x_min, point.x);
                node.x_max = std::max(node.x_max, point.x);
                node.y_min = std::min(node.y_min, point.y);
                node.y_max = std::max(node.y_max, point.y);
            }
            if (node.end - node.begin > maximum_leaf_size
                    && (node.x_min < node.x_max || node.y_min < node.y_max)) {
                bool split_x = (node.x_max - node.x_min >= node.y_max - node.y_min);
                VertexId middle = (node.begin + node.end) / 2;
                std::nth_element(
                        ids_.begin() + node.begin,
                        ids_.begin() + middle,
                        ids_.begin() + node.end,
                        [this, split_x](VertexId vertex_id_1, VertexId vertex_id_2)
                        {
                            return (split_x)?
                                points_[vertex_id_1].x < points_[vertex_id_2].x:
                                points_[vertex_id_1].y < points_[vertex_id_2].y;
                        });
                node.child_1 = nodes_.size();
                nodes_.push_back({node.begin, middle});
                node.child_2 = nodes_.size();
                nodes_.push_back({middle, node.end});
            }
            nodes_[node_id] = node;
        }
    }

    /**
     * Get the 'number_of_neighbors' nearest points to a point in one of the
     * octants around it (octant o covers the angles in [45 o, 45 (o + 1))
     * degrees; points at the same coordinates are in none).
     */
    void octant_neighbors(
            VertexId vertex_id,
            int octant,
            VertexId number_of_neighbors,
            std::vector<VertexId>& neighbors)
    {
        static constexpr double directions[9][2] = {
            {1, 0}, {1, 1}, {0, 1}, {-1, 1}, {-1, 0}, {-1, -1}, {0, -1}, {1, -1}, {1, 0}};
        const double* ray_1 = directions[octant];
        const double* ray_2 = directions[octant + 1];
        const Coordinates2D& point = points_[vertex_id];

        // Best points found, by increasing squared distance.
        best_.clear();
        auto worst = [this, number_of_neighbors]()
        {
            return ((VertexId)best_.size() < number_of_neighbors)?
                std::numeric_limits<double>::infinity():
                best_.back().first;
        };

        stack_.clear();
        stack_.push_back(0);
        while (!stack_.empty()) {
            const Node& node = nodes_[stack_.back()];
            stack_.pop_back();

            // Prune the nodes farther than the worst point found.
            double dx = std::max({node.x_min - point.x, 0.0, point.x - node.x_max});
            double dy = std::max({node.y_min - point.y, 0.0, point.y - node.y_max});
            if (dx * dx + dy * dy >= worst())
                continue;
            // Prune the nodes outside one of the two half-planes delimiting
            // the octant.
            double cross_1 = -std::numeric_limits<double>::infinity();
            double cross_2 = -std::numeric_limits<double>::infinity();
            for (double x: {node.x_min, node.x_max}) {
                for (double y: {node.y_min, node.y_max}) {
                    double xd = x - point.x;
                    double yd = y - point.y;
                    cross_1 = std::max(cross_1, ray_1[0] * yd - ray_1[1] * xd);
                    cross_2 = std::max(cross_2, xd * ray_2[1] - yd * ray_2[0]);
                }
            }
            if (cross_1 < 0 || cross_2 <= 0)
                continue;

            if (node.child_1 != -1) {
                stack_.push_back(node.child_1);
                stack_.push_back(node.child_2);
                continue;
            }
            for (VertexId pos = node.begin; pos < node.end; ++pos) {
                VertexId other_vertex_id = ids_[pos];
                double xd = points_[other_vertex_id].x - point.x;
                double yd = points_[other_vertex_id].y - point.y;
                if (ray_1[0] * yd - ray_1[1] * xd < 0
                        || xd * ray_2[1] - yd * ray_2[0] <= 0)
                    continue;
                double distance = xd * xd + yd * yd;
                if (distance >= worst())
                    continue;
                if ((VertexId)best_.size() == number_of_neighbors)
                    best_.pop_back();
                auto it = std::upper_bound(
                        best_.begin(),
                        best_.end(),
                        std::pair<double, VertexId>(distance, other_vertex_id));
                best_.insert(it, {distance, other_vertex_id});
            }
        }
        for (const auto& p: best_)
            neighbors.push_back(p.second);
    }

private:

    using NodeId = int64_t;

    static constexpr VertexId maximum_leaf_size = 8;

    struct Node
    {
        VertexId begin;
        VertexId end;
        NodeId child_1 = -1;
        NodeId child_2 = -1;
        double x_min = 0;
        double x_max = 0;
        double y_min = 0;
        double y_max = 0;
    };

    const std::vector<Coordinates2D>& points_;

    std::vector<VertexId> ids_;

    std::vector<Node> nodes_;

    std::vector<NodeId> stack_;

    std::vector<std::pair<double, VertexId>> best_;
};

/** Sparse graph, in compressed sparse row format. */
struct Graph
{
    std::vector<VertexId> offsets;

    std::vector<VertexId> neighbors;

    /** Distance of each edge. */
    std::vector<Distance> distances;
};

template <typename T>
Graph make_graph(
        const T& distances,
        VertexId number_of_vertices,
        const std::vector<std::pair<VertexId, VertexId>>& edges);

VertexId find(std::vector<VertexId>& roots, VertexId vertex_id)
{
    while (roots[vertex_id] != vertex_id) {
        roots[vertex_id] = roots[roots[vertex_id]];
        vertex_id = roots[vertex_id];
    }
    return vertex_id;
}

template <typename T>
Graph build_graph(
        const T& distances,
        const Distances& distances_wrapper,
        const AlphaNearnessParameters& parameters)
{
    VertexId number_of_vertices = distances_wrapper.number_of_vertices();
    std::vector<std::pair<VertexId, VertexId>> edges;
    auto add_edge = [&edges](VertexId vertex_id_1, VertexId vertex_id_2)
    {
        if (vertex_id_1 < vertex_id_2)
            edges.push_back({vertex_id_1, vertex_id_2});
        else if (vertex_id_2 < vertex_id_1)
            edges.push_back({vertex_id_2, vertex_id_1});
    };

    // Nearest neighbors.
    if (parameters.number_of_nearest_neighbors > 0) {
        BallTree ball_tree(distances_wrapper);
        for (VertexId vertex_id = 0; vertex_id < number_of_vertices; ++vertex_id) {
            for (VertexId neighbor_id: ball_tree.nearest_neighbors(
                        vertex_id,
                        parameters.number_of_nearest_neighbors)) {
                add_edge(vertex_id, neighbor_id);
            }
        }
    }

    // Octant neighbors.
    if constexpr (HasCoordinates<T>::value) {
        if (parameters.number_of_octant_neighbors > 0) {
            std::vector<Coordinates2D> points(number_of_vertices);
            for (VertexId vertex_id = 0; vertex_id < number_of_vertices; ++vertex_id)
                points[vertex_id] = distances.coordinates(vertex_id);
            KdTree kd_tree(points);
            std::vector<VertexId> neighbors;
            for (VertexId vertex_id = 0; vertex_id < number_of_vertices; ++vertex_id) {
                neighbors.clear();
                for (int octant = 0; octant < 8; ++octant) {
                    kd_tree.octant_neighbors(
                            vertex_id,
                            octant,
                            parameters.number_of_octant_neighbors,
                            neighbors);
                }
                for (VertexId neighbor_id: neighbors)
                    add_edge(vertex_id, neighbor_id);
            }
        }
    }
    std::sort(edges.begin(), edges.end());
    edges.erase(std::unique(edges.begin(), edges.end()), edges.end());

    // Connect the components, if any, by their shortest edges to the rest of
    // the graph (Boruvka rounds, by brute force).
    std::vector<VertexId> roots(number_of_vertices);
    std::iota(roots.begin(), roots.end(), 0);
    VertexId number_of_components = number_of_vertices;
    for (const auto& edge: edges) {
        VertexId root_1 = find(roots, edge.first);
        VertexId root_2 = find(roots, edge.second);
        if (root_1 != root_2) {
            roots[root_1] = root_2;
            number_of_components--;
        }
    }
    while (number_of_components > 1) {
        std::vector<VertexId> components(number_of_vertices);
        for (VertexId vertex_id = 0; vertex_id < number_of_vertices; ++vertex_id)
            components[vertex_id] = find(roots, vertex_id);
        std::vector<Distance> best_distances(number_of_vertices, infinity);
        std::vector<std::pair<VertexId, VertexId>> best_edges(number_of_vertices);
        for (VertexId vertex_id_1 = 0; vertex_id_1 < number_of_vertices; ++vertex_id_1) {
            VertexId component = components[vertex_id_1];
            for (VertexId vertex_id_2 = 0; vertex_id_2 < number_of_vertices; ++vertex_id_2) {
                if (components[vertex_id_2] == component)
                    continue;
                Distance distance = distances.distance(vertex_id_1, vertex_id_2);
                if (best_distances[component] > distance) {
                    best_distances[component] = distance;
                    best_edges[component] = {vertex_id_1, vertex_id_2};
                }
            }
        }
        for (VertexId vertex_id = 0; vertex_id < number_of_vertices; ++vertex_id) {
            if (components[vertex_id] != vertex_id)
                continue;
            auto edge = best_edges[vertex_id];
            add_edge(edge.first, edge.second);
            VertexId root_1 = find(roots, edge.first);
            VertexId root_2 = find(roots, edge.second);
            if (root_1 != root_2) {
                roots[root_1] = root_2;
                number_of_components--;
            }
        }
        std::sort(edges.begin(), edges.end());
        edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
    }
    return make_graph(distances, number_of_vertices, edges);
}

/** Build a graph from its edges (sorted, without duplicates). */
template <typename T>
Graph make_graph(
        const T& distances,
        VertexId number_of_vertices,
        const std::vector<std::pair<VertexId, VertexId>>& edges)
{
    Graph graph;
    graph.offsets.assign(number_of_vertices + 1, 0);
    for (const auto& edge: edges) {
        graph.offsets[edge.first + 1]++;
        graph.offsets[edge.second + 1]++;
    }
    for (VertexId vertex_id = 0; vertex_id < number_of_vertices; ++vertex_id)
        graph.offsets[vertex_id + 1] += graph.offsets[vertex_id];
    graph.neighbors.resize(2 * edges.size());
    graph.distances.resize(2 * edges.size());
    std::vector<VertexId> positions(graph.offsets.begin(), graph.offsets.end() - 1);
    for (const auto& edge: edges) {
        Distance distance = distances.distance(edge.first, edge.second);
        graph.neighbors[positions[edge.first]] = edge.second;
        graph.distances[positions[edge.first]++] = distance;
        graph.neighbors[positions[edge.second]] = edge.first;
        graph.distances[positions[edge.second]++] = distance;
    }
    return graph;
}

/**
 * Minimum 1-tree: a minimum spanning tree, plus an edge between a leaf (the
 * special vertex) and its nearest vertex besides its tree neighbor, the leaf
 * maximizing the cost of that edge.
 */
struct OneTree
{
    /** Parent of each vertex in the tree ('-1' for the root). */
    std::vector<VertexId> parents;

    /** Cost of the edge between each vertex and its parent. */
    std::vector<Distance> parent_costs;

    /** Vertices, parents first. */
    std::vector<VertexId> order;

    /** Special vertex ('-1' if none has an edge besides its tree edge). */
    VertexId special_vertex_id = -1;

    /** Other end of the extra edge of the special vertex. */
    VertexId special_neighbor_id = -1;

    /** Cost of the extra edge. */
    Distance special_cost = 0;

    /** Degrees of the vertices. */
    std::vector<VertexId> degrees;

    /** Total cost. */
    Distance cost = 0;

    /** Is an edge in the 1-tree? */
    bool contains(VertexId vertex_id_1, VertexId vertex_id_2) const
    {
        return parents[vertex_id_1] == vertex_id_2
            || parents[vertex_id_2] == vertex_id_1
            || (vertex_id_1 == special_vertex_id && vertex_id_2 == special_neighbor_id)
            || (vertex_id_2 == special_vertex_id && vertex_id_1 == special_neighbor_id);
    }

    /** Get the tree neighbor of a leaf. */
    VertexId leaf_neighbor(VertexId vertex_id) const
    {
        return (parents[vertex_id] != -1)? parents[vertex_id]: order[1];
    }
};

template <typename T>
class AlphaNearness
{

public:

    AlphaNearness(
            const T& distances,
            const Distances& distances_wrapper,
            const AlphaNearnessParameters& parameters):
        distances_(distances),
        parameters_(parameters),
        number_of_vertices_(distances_wrapper.number_of_vertices()),
        precision_(parameters.precision)
    {
        complete_ = (parameters.number_of_nearest_neighbors < 0
                || parameters.number_of_nearest_neighbors >= number_of_vertices_ - 1);
        dense_ = complete_;
        if (!complete_)
            graph_ = build_graph(distances, distances_wrapper, parameters);
        tree_.parents.resize(number_of_vertices_);
        tree_.parent_costs.resize(number_of_vertices_);
        tree_.degrees.resize(number_of_vertices_);
        keys_.resize(number_of_vertices_);
        heap_positions_.resize(number_of_vertices_);
    }

    AlphaNearnessOutput run()
    {
        AlphaNearnessOutput output;
        std::vector<Distance> pi(number_of_vertices_, 0);
        if (complete_
                && parameters_.number_of_ascent_candidates >= 0
                && parameters_.number_of_ascent_candidates < number_of_vertices_ - 1) {
            // The ascent is done on the graph of the alpha-nearest edges of
            // the first 1-tree.
            compute_one_tree(pi);
            CandidateLists candidates;
            std::vector<std::vector<Distance>> alphas;
            compute_candidates(pi, parameters_.number_of_ascent_candidates, candidates, alphas);
            std::vector<std::pair<VertexId, VertexId>> edges;
            for (VertexId vertex_id = 0; vertex_id < number_of_vertices_; ++vertex_id) {
                for (VertexId other_vertex_id: candidates[vertex_id]) {
                    edges.push_back({
                            std::min(vertex_id, other_vertex_id),
                            std::max(vertex_id, other_vertex_id)});
                }
            }
            std::sort(edges.begin(), edges.end());
            edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
            graph_ = make_graph(distances_, number_of_vertices_, edges);
            dense_ = false;
        }
        ascent(pi, output);
        dense_ = complete_;
        compute_one_tree(output.pi);
        if (complete_) {
            // Bound of the minimum 1-tree of the complete graph.
            Distance pi_sum = 0;
            for (Distance p: output.pi)
                pi_sum += p;
            output.lower_bound = (double)(tree_.cost - 2 * pi_sum) / precision_;
        }
        compute_candidates(
                output.pi,
                parameters_.number_of_candidates,
                output.candidates,
                output.alphas);
        return output;
    }

private:

    inline Distance cost(
            const std::vector<Distance>& pi,
            VertexId vertex_id_1,
            VertexId vertex_id_2,
            Distance distance) const
    {
        return precision_ * distance + pi[vertex_id_1] + pi[vertex_id_2];
    }

    inline Distance cost(
            const std::vector<Distance>& pi,
            VertexId vertex_id_1,
            VertexId vertex_id_2) const
    {
        return cost(pi, vertex_id_1, vertex_id_2, distances_.distance(vertex_id_1, vertex_id_2));
    }

    /** Prim's algorithm, O(n^2). */
    void minimum_spanning_tree_dense(const std::vector<Distance>& pi)
    {
        std::vector<VertexId>& remaining = remaining_;
        remaining.resize(number_of_vertices_ - 1);
        std::iota(remaining.begin(), remaining.end(), 1);
        for (VertexId vertex_id = 0; vertex_id < number_of_vertices_; ++vertex_id)
            keys_[vertex_id] = infinity;
        tree_.parents[0] = -1;
        tree_.parent_costs[0] = 0;
        tree_.order.push_back(0);
        VertexId vertex_id = 0;
        while (!remaining.empty()) {
            VertexId best_pos = -1;
            Distance best_key = infinity;
            for (VertexId pos = 0; pos < (VertexId)remaining.size(); ++pos) {
                VertexId other_vertex_id = remaining[pos];
                Distance c = cost(pi, vertex_id, other_vertex_id);
                if (keys_[other_vertex_id] > c) {
                    keys_[other_vertex_id] = c;
                    tree_.parents[other_vertex_id] = vertex_id;
                }
                if (best_key > keys_[other_vertex_id]) {
                    best_key = keys_[other_vertex_id];
                    best_pos = pos;
                }
            }
            vertex_id = remaining[best_pos];
            remaining[best_pos] = remaining.back();
            remaining.pop_back();
            tree_.parent_costs[vertex_id] = best_key;
            tree_.order.push_back(vertex_id);
        }
    }

    /** Move a vertex up the heap. */
    inline void heap_up(VertexId heap_pos)
    {
        VertexId vertex_id = heap_[heap_pos];
        Distance key = keys_[vertex_id];
        while (heap_pos > 0) {
            VertexId parent_pos = (heap_pos - 1) / 2;
            VertexId parent_vertex_id = heap_[parent_pos];
            if (keys_[parent_vertex_id] <= key)
                break;
            heap_[heap_pos] = parent_vertex_id;
            heap_positions_[parent_vertex_id] = heap_pos;
            heap_pos = parent_pos;
        }
        heap_[heap_pos] = vertex_id;
        heap_positions_[vertex_id] = heap_pos;
    }

    /** Remove the vertex of smallest key from the heap. */
    inline VertexId heap_pop()
    {
        VertexId top_vertex_id = heap_[0];
        VertexId vertex_id = heap_.back();
        heap_.pop_back();
        VertexId heap_size = heap_.size();
        if (heap_size > 0) {
            Distance key = keys_[vertex_id];
            VertexId heap_pos = 0;
            for (;;) {
                VertexId child_pos = 2 * heap_pos + 1;
                if (child_pos >= heap_size)
                    break;
                if (child_pos + 1 < heap_size
                        && keys_[heap_[child_pos + 1]] < keys_[heap_[child_pos]]) {
                    child_pos++;
                }
                VertexId child_vertex_id = heap_[child_pos];
                if (keys_[child_vertex_id] >= key)
                    break;
                heap_[heap_pos] = child_vertex_id;
                heap_positions_[child_vertex_id] = heap_pos;
                heap_pos = child_pos;
            }
            heap_[heap_pos] = vertex_id;
            heap_positions_[vertex_id] = heap_pos;
        }
        return top_vertex_id;
    }

    /** Prim's algorithm with an indexed binary heap, O(m log n). */
    void minimum_spanning_tree_sparse(const std::vector<Distance>& pi)
    {
        // heap_positions_: -1 if never reached, -2 if in the tree.
        std::fill(heap_positions_.begin(), heap_positions_.end(), -1);
        heap_.clear();
        keys_[0] = 0;
        tree_.parents[0] = -1;
        heap_.push_back(0);
        heap_positions_[0] = 0;
        while (!heap_.empty()) {
            VertexId vertex_id = heap_pop();
            heap_positions_[vertex_id] = -2;
            tree_.parent_costs[vertex_id] = keys_[vertex_id];
            tree_.order.push_back(vertex_id);
            Distance pi_vertex = pi[vertex_id];
            for (VertexId edge_pos = graph_.offsets[vertex_id];
                    edge_pos < graph_.offsets[vertex_id + 1];
                    ++edge_pos) {
                VertexId other_vertex_id = graph_.neighbors[edge_pos];
                VertexId heap_pos = heap_positions_[other_vertex_id];
                if (heap_pos == -2)
                    continue;
                Distance c = precision_ * graph_.distances[edge_pos] + pi_vertex + pi[other_vertex_id];
                if (heap_pos == -1) {
                    keys_[other_vertex_id] = c;
                    tree_.parents[other_vertex_id] = vertex_id;
                    heap_.push_back(other_vertex_id);
                    heap_up(heap_.size() - 1);
                } else if (keys_[other_vertex_id] > c) {
                    keys_[other_vertex_id] = c;
                    tree_.parents[other_vertex_id] = vertex_id;
                    heap_up(heap_pos);
                }
            }
        }
    }

    /** Compute a minimum 1-tree for penalties 'pi' in 'tree_'. */
    void compute_one_tree(const std::vector<Distance>& pi)
    {
        tree_.order.clear();
        if (dense_) {
            minimum_spanning_tree_dense(pi);
        } else {
            minimum_spanning_tree_sparse(pi);
        }

        std::fill(tree_.degrees.begin(), tree_.degrees.end(), 0);
        tree_.cost = 0;
        for (VertexId vertex_id = 0; vertex_id < number_of_vertices_; ++vertex_id) {
            if (tree_.parents[vertex_id] == -1)
                continue;
            tree_.degrees[vertex_id]++;
            tree_.degrees[tree_.parents[vertex_id]]++;
            tree_.cost += tree_.parent_costs[vertex_id];
        }

        // Special vertex.
        tree_.special_vertex_id = -1;
        tree_.special_neighbor_id = -1;
        tree_.special_cost = -infinity;
        for (VertexId vertex_id = 0; vertex_id < number_of_vertices_; ++vertex_id) {
            if (tree_.degrees[vertex_id] != 1)
                continue;
            VertexId tree_neighbor_id = tree_.leaf_neighbor(vertex_id);
            VertexId best_neighbor_id = -1;
            Distance best_cost = infinity;
            if (dense_) {
                for (VertexId other_vertex_id = 0;
                        other_vertex_id < number_of_vertices_;
                        ++other_vertex_id) {
                    if (other_vertex_id == vertex_id || other_vertex_id == tree_neighbor_id)
                        continue;
                    Distance c = cost(pi, vertex_id, other_vertex_id);
                    if (best_cost > c) {
                        best_cost = c;
                        best_neighbor_id = other_vertex_id;
                    }
                }
            } else {
                for (VertexId edge_pos = graph_.offsets[vertex_id];
                        edge_pos < graph_.offsets[vertex_id + 1];
                        ++edge_pos) {
                    VertexId other_vertex_id = graph_.neighbors[edge_pos];
                    if (other_vertex_id == tree_neighbor_id)
                        continue;
                    Distance c = cost(pi, vertex_id, other_vertex_id, graph_.distances[edge_pos]);
                    if (best_cost > c) {
                        best_cost = c;
                        best_neighbor_id = other_vertex_id;
                    }
                }
            }
            if (best_neighbor_id != -1 && tree_.special_cost < best_cost) {
                tree_.special_vertex_id = vertex_id;
                tree_.special_neighbor_id = best_neighbor_id;
                tree_.special_cost = best_cost;
            }
        }
        if (tree_.special_vertex_id != -1) {
            tree_.degrees[tree_.special_vertex_id]++;
            tree_.degrees[tree_.special_neighbor_id]++;
            tree_.cost += tree_.special_cost;
        }
    }

    /**
     * Subgradient ascent: the penalties of the vertices are moved in the
     * direction of their degree minus 2 in the 1-tree, with step sizes
     * following Helsgaun (2000):
     * - the step size is constant during each period; at the end of a
     *   period, both are halved;
     * - during the initial phase, at the beginning of the first period, the
     *   step size is doubled after each improvement of the bound; the phase
     *   ends at the first iteration without improvement in the second half
     *   of the period, which then restarts with a step size reduced by a
     *   quarter;
     * - a period is doubled (up to the initial period) when its last
     *   iteration improved the bound.
     * An iteration improves the bound if it increases it, or leaves it
     * unchanged with a 1-tree closer to a tour (smaller norm of the degrees
     * minus 2).
     */
    void ascent(
            std::vector<Distance>& pi,
            AlphaNearnessOutput& output)
    {
        std::vector<VertexId> v(number_of_vertices_);
        std::vector<VertexId> v_previous(number_of_vertices_);
        Distance pi_sum = 0;

        // Compute a 1-tree, and the degrees minus 2; return the norm of the
        // latter.
        auto iteration = [&]()
        {
            compute_one_tree(pi);
            output.number_of_iterations++;
            Distance norm = 0;
            for (VertexId vertex_id = 0; vertex_id < number_of_vertices_; ++vertex_id) {
                v[vertex_id] = tree_.degrees[vertex_id] - 2;
                norm += v[vertex_id] * v[vertex_id];
            }
            return norm;
        };
        auto bound = [&]() { return tree_.cost - 2 * pi_sum; };

        Distance norm = iteration();
        Distance best_bound = bound();
        Distance best_norm = norm;
        output.pi = pi;
        v_previous = v;

        int64_t initial_period = parameters_.initial_period;
        if (initial_period < 0)
            initial_period = std::max((int64_t)100, (int64_t)number_of_vertices_ / 2);
        int64_t period = initial_period;
        Distance step = parameters_.initial_step_size * precision_;
        bool initial_phase = true;
        auto limit_reached = [&]()
        {
            return parameters_.maximum_number_of_iterations >= 0
                && output.number_of_iterations >= parameters_.maximum_number_of_iterations;
        };
        // Penalties are bounded to avoid overflows.
        const Distance maximum_pi = std::numeric_limits<Distance>::max() / (8 * number_of_vertices_);
        while (norm != 0 && step > 0 && period > 0 && !limit_reached()) {
            for (int64_t p = 1;
                    norm != 0 && step > 0 && p <= period && !limit_reached();
                    ++p) {
                for (VertexId vertex_id = 0; vertex_id < number_of_vertices_; ++vertex_id) {
                    Distance delta = step * (7 * v[vertex_id] + 3 * v_previous[vertex_id]) / 10;
                    Distance pi_new = std::max(-maximum_pi, std::min(maximum_pi, pi[vertex_id] + delta));
                    pi_sum += pi_new - pi[vertex_id];
                    pi[vertex_id] = pi_new;
                }
                v_previous = v;
                norm = iteration();
                Distance b = bound();
                if (best_bound < b || (best_bound == b && norm < best_norm)) {
                    best_bound = b;
                    best_norm = norm;
                    output.pi = pi;
                    if (initial_phase)
                        step *= 2;
                    if (p == period)
                        period = std::min(2 * period, initial_period);
                } else if (initial_phase && p > period / 2) {
                    initial_phase = false;
                    p = 0;
                    step = 3 * step / 4;
                }
            }
            initial_phase = false;
            step /= 2;
            period /= 2;
        }
        output.tour_found = (best_norm == 0);
        output.lower_bound = (double)best_bound / precision_;
    }

    /** Compute the alpha values and the candidates from 'tree_'. */
    void compute_candidates(
            const std::vector<Distance>& pi,
            VertexId number_of_candidates,
            CandidateLists& candidates,
            std::vector<std::vector<Distance>>& alphas)
    {
        number_of_candidates = std::min(number_of_candidates, number_of_vertices_ - 1);
        candidates.assign(number_of_vertices_, {});
        alphas.assign(number_of_vertices_, {});
        if (number_of_candidates <= 0)
            return;

        // (alpha, cost, vertex).
        std::vector<std::tuple<Distance, Distance, VertexId>> values;
        auto select = [&](VertexId vertex_id)
        {
            VertexId k = std::min(number_of_candidates, (VertexId)values.size());
            std::partial_sort(values.begin(), values.begin() + k, values.end());
            for (VertexId pos = 0; pos < k; ++pos) {
                candidates[vertex_id].push_back(std::get<2>(values[pos]));
                alphas[vertex_id].push_back(std::get<0>(values[pos]));
            }
        };
        // Alpha of an edge incident to the special vertex.
        auto special_alpha = [this](VertexId vertex_id_1, VertexId vertex_id_2, Distance c)
        {
            return (tree_.contains(vertex_id_1, vertex_id_2))? 0: c - tree_.special_cost;
        };
        VertexId special_vertex_id = tree_.special_vertex_id;

        if (dense_) {
            // For each vertex i, beta(i, j) for all j: the largest cost on the
            // tree path between them, computed from the path from i to the
            // root, then in the order of the tree for the other vertices.
            std::vector<Distance> beta(number_of_vertices_);
            std::vector<VertexId> marks(number_of_vertices_, -1);
            for (VertexId vertex_id = 0; vertex_id < number_of_vertices_; ++vertex_id) {
                values.clear();
                if (vertex_id == special_vertex_id) {
                    for (VertexId other_vertex_id = 0;
                            other_vertex_id < number_of_vertices_;
                            ++other_vertex_id) {
                        if (other_vertex_id == vertex_id)
                            continue;
                        Distance c = cost(pi, vertex_id, other_vertex_id);
                        values.push_back({special_alpha(vertex_id, other_vertex_id, c), c, other_vertex_id});
                    }
                    select(vertex_id);
                    continue;
                }
                beta[vertex_id] = -infinity;
                marks[vertex_id] = vertex_id;
                for (VertexId current_vertex_id = vertex_id;
                        tree_.parents[current_vertex_id] != -1;
                        current_vertex_id = tree_.parents[current_vertex_id]) {
                    VertexId parent_id = tree_.parents[current_vertex_id];
                    beta[parent_id] = std::max(
                            beta[current_vertex_id],
                            tree_.parent_costs[current_vertex_id]);
                    marks[parent_id] = vertex_id;
                }
                for (VertexId other_vertex_id: tree_.order) {
                    if (marks[other_vertex_id] != vertex_id) {
                        beta[other_vertex_id] = std::max(
                                beta[tree_.parents[other_vertex_id]],
                                tree_.parent_costs[other_vertex_id]);
                    }
                    if (other_vertex_id == vertex_id)
                        continue;
                    Distance c = cost(pi, vertex_id, other_vertex_id);
                    Distance alpha = (other_vertex_id == special_vertex_id)?
                        special_alpha(vertex_id, other_vertex_id, c):
                        c - beta[other_vertex_id];
                    values.push_back({alpha, c, other_vertex_id});
                }
                select(vertex_id);
            }
            return;
        }

        // Sparse graph: beta by binary lifting.
        int number_of_levels = 1;
        while (((VertexId)1 << number_of_levels) < number_of_vertices_)
            number_of_levels++;
        std::vector<VertexId> depths(number_of_vertices_, 0);
        std::vector<VertexId> ancestors(number_of_levels * number_of_vertices_);
        std::vector<Distance> maximums(number_of_levels * number_of_vertices_);
        for (VertexId vertex_id: tree_.order) {
            VertexId parent_id = tree_.parents[vertex_id];
            if (parent_id == -1) {
                ancestors[vertex_id] = vertex_id;
                maximums[vertex_id] = -infinity;
            } else {
                depths[vertex_id] = depths[parent_id] + 1;
                ancestors[vertex_id] = parent_id;
                maximums[vertex_id] = tree_.parent_costs[vertex_id];
            }
        }
        for (int level = 1; level < number_of_levels; ++level) {
            VertexId* a = ancestors.data() + level * number_of_vertices_;
            const VertexId* a_previous = a - number_of_vertices_;
            Distance* m = maximums.data() + level * number_of_vertices_;
            const Distance* m_previous = m - number_of_vertices_;
            for (VertexId vertex_id = 0; vertex_id < number_of_vertices_; ++vertex_id) {
                VertexId middle = a_previous[vertex_id];
                a[vertex_id] = a_previous[middle];
                m[vertex_id] = std::max(m_previous[vertex_id], m_previous[middle]);
            }
        }
        auto beta = [&](VertexId vertex_id_1, VertexId vertex_id_2)
        {
            Distance result = -infinity;
            if (depths[vertex_id_1] < depths[vertex_id_2])
                std::swap(vertex_id_1, vertex_id_2);
            VertexId difference = depths[vertex_id_1] - depths[vertex_id_2];
            for (int level = 0; difference > 0; ++level, difference >>= 1) {
                if (difference & 1) {
                    result = std::max(result, maximums[level * number_of_vertices_ + vertex_id_1]);
                    vertex_id_1 = ancestors[level * number_of_vertices_ + vertex_id_1];
                }
            }
            if (vertex_id_1 == vertex_id_2)
                return result;
            for (int level = number_of_levels - 1; level >= 0; --level) {
                VertexId a1 = ancestors[level * number_of_vertices_ + vertex_id_1];
                VertexId a2 = ancestors[level * number_of_vertices_ + vertex_id_2];
                if (a1 != a2) {
                    result = std::max({
                            result,
                            maximums[level * number_of_vertices_ + vertex_id_1],
                            maximums[level * number_of_vertices_ + vertex_id_2]});
                    vertex_id_1 = a1;
                    vertex_id_2 = a2;
                }
            }
            return std::max({result, maximums[vertex_id_1], maximums[vertex_id_2]});
        };
        for (VertexId vertex_id = 0; vertex_id < number_of_vertices_; ++vertex_id) {
            values.clear();
            for (VertexId edge_pos = graph_.offsets[vertex_id];
                    edge_pos < graph_.offsets[vertex_id + 1];
                    ++edge_pos) {
                VertexId other_vertex_id = graph_.neighbors[edge_pos];
                Distance c = cost(pi, vertex_id, other_vertex_id, graph_.distances[edge_pos]);
                Distance alpha = 0;
                if (vertex_id == special_vertex_id || other_vertex_id == special_vertex_id) {
                    alpha = special_alpha(vertex_id, other_vertex_id, c);
                } else if (!tree_.contains(vertex_id, other_vertex_id)) {
                    alpha = c - beta(vertex_id, other_vertex_id);
                }
                values.push_back({alpha, c, other_vertex_id});
            }
            select(vertex_id);
        }
    }

    const T& distances_;

    const AlphaNearnessParameters& parameters_;

    VertexId number_of_vertices_;

    Distance precision_;

    /** Is the graph complete? */
    bool complete_;

    /** Are the 1-trees computed on the complete graph (or on 'graph_')? */
    bool dense_;

    Graph graph_;

    OneTree tree_;

    std::vector<Distance> keys_;

    /** Position of each vertex in 'heap_'. */
    std::vector<VertexId> heap_positions_;

    std::vector<VertexId> remaining_;

    /** Binary heap of the vertices, by key. */
    std::vector<VertexId> heap_;
};

template <typename T>
AlphaNearnessOutput compute_alpha_nearness_candidates(
        const T& distances,
        const Distances& distances_wrapper,
        const AlphaNearnessParameters& parameters)
{
    VertexId number_of_vertices = distances_wrapper.number_of_vertices();
    if (number_of_vertices <= 2) {
        AlphaNearnessOutput output;
        output.candidates.resize(number_of_vertices);
        output.alphas.resize(number_of_vertices);
        output.pi.assign(number_of_vertices, 0);
        if (parameters.number_of_candidates > 0) {
            for (VertexId vertex_id_1 = 0; vertex_id_1 < number_of_vertices; ++vertex_id_1) {
                for (VertexId vertex_id_2 = 0; vertex_id_2 < number_of_vertices; ++vertex_id_2) {
                    if (vertex_id_2 == vertex_id_1)
                        continue;
                    output.candidates[vertex_id_1].push_back(vertex_id_2);
                    output.alphas[vertex_id_1].push_back(0);
                }
            }
        }
        return output;
    }
    AlphaNearness<T> alpha_nearness(distances, distances_wrapper, parameters);
    return alpha_nearness.run();
}

}

AlphaNearnessOutput travelingsalesmansolver::alpha_nearness_candidates(
        const Distances& distances,
        const AlphaNearnessParameters& parameters)
{
    // Coordinate-based distances first, for the octant neighbors.
    return FUNCTION_WITH_DISTANCES_R(
            compute_alpha_nearness_candidates,
            distances,
            distances,
            parameters);
}
