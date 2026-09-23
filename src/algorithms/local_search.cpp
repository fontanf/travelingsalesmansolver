// The EAX genetic algorithm implemented in 'local_search.hpp' is adapted from
// Shujia Liu's C++ implementation (https://github.com/Sugia/GA-for-TSP),
// licensed under the Apache License, Version 2.0; see
// 'licenses/eax-ga/LICENSE' and 'licenses/eax-ga/NOTICE.md' for the license
// text and the full list of changes made to the original source.

#include "travelingsalesmansolver/algorithms/local_search.hpp"

using namespace travelingsalesmansolver;

const Output travelingsalesmansolver::local_search(
        const Instance& instance,
        const LocalSearchParameters& parameters)
{
    return FUNCTION_WITH_DISTANCES(
            local_search,
            instance.distances(),
            instance,
            parameters);
}
