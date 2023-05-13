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

    /** SEED */
    std::string seed;

    /** CANDIDATE_FILE content. */
    std::string candidate_file_content;

    /** Info structure. */
    optimizationtools::Info info = optimizationtools::Info();
};

struct LkhOutput: Output
{
    LkhOutput(
            const Instance& instance,
            optimizationtools::Info& info):
        Output(instance, info) { }

    //void print_statistics(
    //        optimizationtools::Info& info) const override;

    std::string candidate_file_content;
};

LkhOutput lkh(
        const Instance& instance,
        LkhOptionalParameters parameters = {});

}

