#include "travelingsalesmansolver/algorithms/greedy.hpp"

using namespace travelingsalesmansolver;

const Output travelingsalesmansolver::greedy(
        const Instance& instance,
        const GreedyParameters& parameters,
        const CandidateLists* candidates)
{
    return FUNCTION_WITH_DISTANCES(
            greedy,
            instance.distances(),
            instance,
            parameters,
            candidates);
}
