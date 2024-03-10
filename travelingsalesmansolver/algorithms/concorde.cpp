#include "travelingsalesmansolver/algorithms/concorde.hpp"

using namespace travelingsalesmansolver;

const Output travelingsalesmansolver::concorde(
        const Instance& instance,
        const Parameters& parameters)
{
    return FUNCTION_WITH_DISTANCES(
            concorde,
            instance.distances(),
            instance,
            parameters);
}
