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

    /** MAX_CANDIDATES. */
    std::string max_candidates;

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

const LkhOutput lkh(
        const Instance& instance,
        LkhOptionalParameters parameters = {});

struct LkhCandidateEdge
{
    /** Other end of the edge. */
    VertexId vertex_id;

    /** Alpha value. */
    double alpha;
};

struct LkhCandidate
{
    /**
     * Parent of the vertex in the minimum spanning tree.
     *
     * -1 if the vertex has no parent.
     */
    VertexId parent_id = -1;

    /** Candidate edges. */
    std::vector<LkhCandidateEdge> edges;
};

std::vector<LkhCandidate> read_candidates(std::string candidate_file_content);

}

