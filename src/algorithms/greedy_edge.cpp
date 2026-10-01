#include "travelingsalesmansolver/algorithms/greedy_edge.hpp"

using namespace travelingsalesmansolver;

const Output travelingsalesmansolver::greedy_edge(
        const Instance& instance,
        const GreedyEdgeParameters& parameters,
        const CandidateLists* candidates)
{
    return FUNCTION_WITH_DISTANCES(
            greedy_edge,
            instance.distances(),
            instance,
            parameters,
            candidates);
}
