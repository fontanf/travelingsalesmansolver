#pragma once

#include "travelingsalesmansolver/solution.hpp"

namespace travelingsalesmansolver
{

struct LkhOptionalParameters
{
    /** CANDIDATE_SET_TYPE */
    std::string candidate_set_type;

    /** INITIAL_PERIOD */
    std::string initial_period;

    /** RUNS */
    std::string runs;

    /** MAX_TRIALS */
    std::string max_trials;

    /** Info structure. */
    optimizationtools::Info info = optimizationtools::Info();
};

Output lkh(
        const Instance& instance,
        LkhOptionalParameters parameters = {});

}

