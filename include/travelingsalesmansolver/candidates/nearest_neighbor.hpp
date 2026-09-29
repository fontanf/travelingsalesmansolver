#pragma once

#include "travelingsalesmansolver/candidates/candidate_lists.hpp"
#include "travelingsalesmansolver/distances/distances.hpp"

namespace travelingsalesmansolver
{

/**
 * Compute, for each vertex, its 'number_of_candidates' nearest neighbors,
 * nearest first (ball tree search; approximate when the distances don't
 * satisfy the triangle inequality).
 */
CandidateLists nearest_neighbor_candidates(
        const Distances& distances,
        VertexId number_of_candidates);

}
