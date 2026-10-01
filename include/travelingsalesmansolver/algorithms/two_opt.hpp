#pragma once

#include "travelingsalesmansolver/algorithm_formatter.hpp"
#include "travelingsalesmansolver/solution_builder.hpp"
#include "travelingsalesmansolver/candidates/nearest_neighbor.hpp"
#include "travelingsalesmansolver/lin_kernighan/two_level_list.hpp"
#include "travelingsalesmansolver/algorithms/random_permutation.hpp"

#include <numeric>
#include <random>

namespace travelingsalesmansolver
{

struct TwoOptParameters: Parameters
{
    /**
     * Number of candidate edges per vertex (nearest neighbors), if no
     * candidates are given (as the 2-opt of the original EAX; from a random
     * tour, 10 gave tours about twice as long as the optimum on d1291).
     */
    VertexId number_of_candidates = 50;


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
 * 2-opt local search: apply improving 2-opt moves until none remains. The
 * moves remove a tour edge (t1, t2) and add an edge (t1, t4) to a candidate
 * t4 of t1 (the first improving one is applied); the vertices are processed
 * with don't-look bits (a vertex is reconsidered when a move removes or adds
 * an edge at a vertex having it among its candidates).
 *
 * The generator gives the initial tour (if none is given) and the order in
 * which the vertices are processed.
 *
 * 'initial_solution', if provided, is the initial tour (a
 * 'random_permutation' otherwise). It must be a feasible solution of the same
 * 'instance'.
 *
 * 'candidates', if provided, are the candidate edges of each vertex (the
 * 'parameters.number_of_candidates' nearest neighbors otherwise).
 */
const Output two_opt(
        const Instance& instance,
        std::mt19937_64& generator,
        const TwoOptParameters& parameters = {},
        const Solution* initial_solution = nullptr,
        const CandidateLists* candidates = nullptr);

template <typename Distances>
const Output two_opt(
        const Distances& distances,
        const Instance& instance,
        std::mt19937_64& generator,
        const TwoOptParameters& parameters = {},
        const Solution* initial_solution = nullptr,
        const CandidateLists* candidates = nullptr)
{
    Output output(instance);
    AlgorithmFormatter algorithm_formatter(parameters, output);
    algorithm_formatter.start("2-opt");
    algorithm_formatter.print_header();

    VertexId number_of_vertices = instance.number_of_vertices();

    // Initial tour.
    std::vector<VertexId> initial_tour;
    if (initial_solution != nullptr) {
        for (VertexId pos = 0; pos < initial_solution->number_of_vertices(); ++pos)
            initial_tour.push_back(initial_solution->vertex_id(pos));
    } else {
        RandomPermutationParameters random_permutation_parameters;
        random_permutation_parameters.verbosity_level = 0;
        Output random_permutation_output = random_permutation(instance, generator, random_permutation_parameters);
        for (VertexId pos = 0; pos < number_of_vertices; ++pos)
            initial_tour.push_back(random_permutation_output.solution.vertex_id(pos));
    }

    CandidateLists default_candidates;
    if (candidates == nullptr) {
        default_candidates = nearest_neighbor_candidates(
                instance.distances(),
                parameters.number_of_candidates);
        candidates = &default_candidates;
    }

    if (number_of_vertices >= 4) {
        // Vertices having each vertex among their candidates.
        std::vector<std::vector<VertexId>> inverse_candidates(number_of_vertices);
        for (VertexId vertex_id = 0; vertex_id < number_of_vertices; ++vertex_id)
            for (VertexId candidate_id: (*candidates)[vertex_id])
                inverse_candidates[candidate_id].push_back(vertex_id);

        TwoLevelList tour(initial_tour);
        // Queue of the active vertices, initially all of them, in random
        // order.
        std::vector<VertexId> queue(number_of_vertices);
        std::iota(queue.begin(), queue.end(), 0);
        std::shuffle(queue.begin(), queue.end(), generator);
        std::vector<char> in_queue(number_of_vertices, true);
        VertexId queue_start = 0;
        auto activate = [&](VertexId vertex_id)
        {
            for (VertexId other_vertex_id: inverse_candidates[vertex_id]) {
                if (!in_queue[other_vertex_id]) {
                    in_queue[other_vertex_id] = true;
                    queue.push_back(other_vertex_id);
                }
            }
            if (!in_queue[vertex_id]) {
                in_queue[vertex_id] = true;
                queue.push_back(vertex_id);
            }
        };
        while (queue_start < (VertexId)queue.size()) {
            VertexId t1 = queue[queue_start];
            queue_start++;
            in_queue[t1] = false;
            bool improved = false;
            for (int side = 0; side < 2 && !improved; ++side) {
                // Remove (t1, t2) and (t4, t3), add (t1, t4) and (t2, t3);
                // t3 is on the same side of t4 as t2 of t1.
                VertexId t2 = (side == 0)? tour.next(t1): tour.previous(t1);
                Distance gain_1 = distances.distance(t1, t2);
                for (VertexId t4: (*candidates)[t1]) {
                    Distance gain_2 = gain_1 - distances.distance(t1, t4);
                    if (gain_2 <= 0)
                        continue;
                    if (t4 == t2)
                        continue;
                    VertexId t3 = (side == 0)? tour.next(t4): tour.previous(t4);
                    if (t3 == t1)
                        continue;
                    Distance gain = gain_2
                        + distances.distance(t4, t3)
                        - distances.distance(t2, t3);
                    if (gain > 0) {
                        tour.two_opt_move(t1, t2, t3, t4);
                        for (VertexId vertex_id: {t1, t2, t3, t4})
                            activate(vertex_id);
                        improved = true;
                        break;
                    }
                }
            }
            if (queue_start > 4096 && 2 * queue_start > (VertexId)queue.size()) {
                queue.erase(queue.begin(), queue.begin() + queue_start);
                queue_start = 0;
            }
        }
        initial_tour = tour.tour();
    }

    SolutionBuilder solution_builder(instance);
    for (VertexId vertex_id: initial_tour)
        solution_builder.add_vertex(vertex_id);
    algorithm_formatter.update_solution(solution_builder.build(), "");

    algorithm_formatter.end();
    return output;
}

}
