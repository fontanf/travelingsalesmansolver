#include "travelingsalesmansolver/algorithms/eax.hpp"

using namespace travelingsalesmansolver;

const Output travelingsalesmansolver::eax(
        const Instance& instance,
        const EaxParameters& parameters)
{
    return FUNCTION_WITH_DISTANCES(
            eax,
            instance.distances(),
            instance,
            parameters);
}
