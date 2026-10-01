#pragma once

#include "travelingsalesmansolver/algorithm_formatter.hpp"
#include "travelingsalesmansolver/candidates/candidate_lists.hpp"

#include <random>

namespace travelingsalesmansolver
{

struct RandomWalkParameters: Parameters
{
    /** Number of candidate edges per vertex, if no candidates are given. */
    VertexId number_of_candidates = 5;

    /**
     * Candidates, if none are given: "alpha-nearness"
     * ('alpha_nearness_candidates') or "nearest-neighbor"
     * ('nearest_neighbor_candidates').
     */
    std::string candidates = "alpha-nearness";


    virtual nlohmann::json to_json() const override
    {
        nlohmann::json json = Parameters::to_json();
        json.merge_patch({
                {"NumberOfCandidates", number_of_candidates},
                {"Candidates", candidates},
                });
        return json;
    }

    virtual int format_width() const override { return 37; }

    virtual void format(std::ostream& os) const override
    {
        Parameters::format(os);
        int width = format_width();
        os
            << std::setw(width) << std::left << "Number of candidates: " << number_of_candidates << std::endl
            << std::setw(width) << std::left << "Candidates: " << candidates << std::endl
            ;
    }
};

/**
 * Random walk on the candidate graph (as LKH's 'WALK' initial tour): from a
 * random vertex, go to a random unvisited candidate of the current vertex,
 * or to a random unvisited vertex if there is none.
 *
 * 'candidates', if provided, are the candidate edges of each vertex;
 * otherwise, they are computed as given by 'parameters.candidates'.
 */
const Output random_walk(
        const Instance& instance,
        std::mt19937_64& generator,
        const RandomWalkParameters& parameters = {},
        const CandidateLists* candidates = nullptr);

}
