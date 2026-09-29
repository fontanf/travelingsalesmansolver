#include "travelingsalesmansolver/candidates/nearest_neighbor.hpp"

#include "travelingsalesmansolver/distances/ball_tree.hpp"

using namespace travelingsalesmansolver;

CandidateLists travelingsalesmansolver::nearest_neighbor_candidates(
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
