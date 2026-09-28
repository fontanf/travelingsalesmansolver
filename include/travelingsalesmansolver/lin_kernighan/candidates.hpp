#pragma once

#include "travelingsalesmansolver/distances/ball_tree.hpp"

#include <vector>

namespace travelingsalesmansolver
{

/**
 * Candidate lists: for each vertex, the vertices its edges are tried with,
 * sorted by increasing distance.
 */
using CandidateLists = std::vector<std::vector<VertexId>>;

/**
 * Compute, for each vertex, its 'number_of_candidates' nearest neighbors
 * (ball tree search; approximate when the distances don't satisfy the
 * triangle inequality).
 */
inline CandidateLists nearest_neighbor_candidates(
        const Distances& distances,
        VertexId number_of_candidates)
{
    VertexId number_of_vertices = distances.number_of_vertices();
    if (number_of_candidates > number_of_vertices - 1)
        number_of_candidates = number_of_vertices - 1;
    CandidateLists candidates(number_of_vertices);
    if (number_of_candidates <= 0)
        return candidates;
    BallTree ball_tree(distances);
    for (VertexId vertex_id = 0; vertex_id < number_of_vertices; ++vertex_id)
        candidates[vertex_id] = ball_tree.nearest_neighbors(vertex_id, number_of_candidates);
    return candidates;
}

}
