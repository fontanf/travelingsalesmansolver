#pragma once

#include "travelingsalesmansolver/algorithm_formatter.hpp"
#include "travelingsalesmansolver/solution_builder.hpp"
#include "travelingsalesmansolver/candidates/nearest_neighbor.hpp"
#include "travelingsalesmansolver/lin_kernighan/greedy.hpp"

namespace travelingsalesmansolver
{

struct GreedyParameters: Parameters
{
    /** Number of candidate edges per vertex (nearest neighbors). */
    VertexId number_of_candidates = 10;


    virtual nlohmann::json to_json() const override
    {
        nlohmann::json json = Parameters::to_json();
        json.merge_patch({
                {"NumberOfCandidates", number_of_candidates},
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
            ;
    }
};

/**
 * Greedy edge (multiple fragment) tour: the candidate edges (nearest
 * neighbors) are added in order of increasing length, skipping those which
 * would give a vertex a third edge or close a cycle; the fragments are then
 * joined (see 'greedy_edge_tour').
 *
 * 'candidates', if provided, are the candidate edges (the
 * 'parameters.number_of_candidates' nearest neighbors otherwise).
 */
const Output greedy(
        const Instance& instance,
        const GreedyParameters& parameters = {},
        const CandidateLists* candidates = nullptr);

template <typename Distances>
const Output greedy(
        const Distances& distances,
        const Instance& instance,
        const GreedyParameters& parameters = {},
        const CandidateLists* candidates = nullptr)
{
    Output output(instance);
    AlgorithmFormatter algorithm_formatter(parameters, output);
    algorithm_formatter.start("Greedy");
    algorithm_formatter.print_header();

    CandidateLists nearest_neighbors;
    if (candidates == nullptr) {
        nearest_neighbors = nearest_neighbor_candidates(
                instance.distances(),
                parameters.number_of_candidates);
        candidates = &nearest_neighbors;
    }
    std::vector<VertexId> tour = greedy_edge_tour(distances, *candidates);
    SolutionBuilder solution_builder(instance);
    for (VertexId vertex_id: tour)
        solution_builder.add_vertex(vertex_id);
    algorithm_formatter.update_solution(solution_builder.build(), "");

    algorithm_formatter.end();
    return output;
}

}
