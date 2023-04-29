#pragma once

#include "travelingsalesmansolver/solution.hpp"

#include "travelingsalesmansolver/algorithms/concorde.hpp"
#include "travelingsalesmansolver/algorithms/lkh.hpp"

namespace travelingsalesmansolver
{

Output run(
        std::string algorithm,
        const Instance& instance,
        const Solution& initial_solution,
        std::mt19937_64& generator,
        optimizationtools::Info info);

}

