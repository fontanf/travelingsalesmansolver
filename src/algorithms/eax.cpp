// The EAX genetic algorithm implemented in 'eax.hpp' is adapted from
// Shujia Liu's C++ implementation (https://github.com/Sugia/GA-for-TSP),
// licensed under the Apache License, Version 2.0; see
// 'licenses/eax-ga/LICENSE' and 'licenses/eax-ga/NOTICE.md' for the license
// text and the full list of changes made to the original source.

#include "travelingsalesmansolver/algorithms/eax.hpp"

using namespace travelingsalesmansolver;

const Output travelingsalesmansolver::eax(
        const Instance& instance,
        std::mt19937_64& generator,
        const EaxParameters& parameters)
{
    return FUNCTION_WITH_DISTANCES(
            eax,
            instance.distances(),
            instance,
            generator,
            parameters);
}
