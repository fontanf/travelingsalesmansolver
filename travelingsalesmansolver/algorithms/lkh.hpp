#pragma once

#include "travelingsalesmansolver/solution.hpp"

#include <iomanip>

namespace travelingsalesmansolver
{

struct LkhParameters: Parameters
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


    virtual int format_width() const override { return 37; }

    virtual void format(std::ostream& os) const override
    {
        Parameters::format(os);
        int width = format_width();
        os
            << std::setw(width) << std::left << "Candidate set type: " << candidate_set_type << std::endl
            << std::setw(width) << std::left << "Initial period:  " << initial_period << std::endl
            << std::setw(width) << std::left << "Runs:  " << runs << std::endl
            << std::setw(width) << std::left << "Seed:  " << seed << std::endl
            << std::setw(width) << std::left << "Has candidate flle content:  " << (!candidate_file_content.empty()) << std::endl
            << std::setw(width) << std::left << "Max candidates:  " << max_candidates << std::endl
            ;
    }

    virtual nlohmann::json to_json() const override
    {
        nlohmann::json json = Parameters::to_json();
        json.merge_patch({
                {"CandidateSetType", candidate_set_type},
                {"InitialPeriod", initial_period},
                {"seed", seed},
                {"HasCandidateFileContent", (!candidate_file_content.empty())},
                {"MaxCandidates", max_candidates},
                });
        return json;
    }
};

struct LkhOutput: Output
{
    LkhOutput(const Instance& instance):
        Output(instance) { }


    std::string candidate_file_content;
};

const LkhOutput lkh(
        const Instance& instance,
        const LkhParameters& parameters = {});

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

std::vector<LkhCandidate> read_candidates(
        const std::string& candidate_file_content);

}

