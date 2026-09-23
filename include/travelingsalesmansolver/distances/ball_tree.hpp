#pragma once

#include "travelingsalesmansolver/distances/distances.hpp"

#include "optimizationtools/containers/indexed_binary_heap.hpp"

namespace travelingsalesmansolver
{

class BallTree
{

public:

    /** Constructor. */
    BallTree(const Distances& distances);

    /**
     * Get the k nearest neighbors of a point (excluding the point itself),
     * sorted by increasing distance.
     *
     * The search relies on the triangle inequality; if the distances don't
     * satisfy it, the returned neighbors are only approximate.
     *
     * Fewer than k vertices are returned if the instance has fewer than k + 1
     * vertices.
     */
    std::vector<VertexId> nearest_neighbors(
            VertexId vertex_id,
            VertexId number_of_neighbors);

private:

    using NodeId = int64_t;

    struct Node
    {
        VertexId center_vertex_id = -1;

        Distance median_distance = 0;

        Distance radius = 0;

        NodeId lesser_child_id = -1;

        NodeId greater_child_id = -1;

        std::vector<VertexId> vertex_ids;
    };

    struct StackElement
    {
        NodeId node_id;
        std::vector<VertexId> vertex_ids;
    };

    /** Build the tree. */
    void build_tree();

    /** Build the tree. */
    template <typename T>
    void build_tree(
            const T& distances);

    /** Get the k nearest neighbors of a point. */
    template <typename T>
    std::vector<VertexId> nearest_neighbors(
            const T& distances,
            VertexId vertex_id,
            VertexId number_of_neighbors);

    /** Distances. */
    const Distances& distances_;

    /** Tree. */
    std::vector<Node> tree_;

    optimizationtools::IndexedBinaryHeap<Distance> neighbors_tmp_;
};

}
