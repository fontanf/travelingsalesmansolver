#pragma once

#include "travelingsalesmansolver/algorithm_formatter.hpp"
#include "travelingsalesmansolver/solution_builder.hpp"
#include "travelingsalesmansolver/lin_kernighan/engine.hpp"
#include "travelingsalesmansolver/algorithms/greedy.hpp"

namespace travelingsalesmansolver
{

struct LinKernighanParameters: Parameters
{
    /** Number of candidate edges per vertex (nearest neighbors). */
    VertexId number_of_candidates = 10;

    /**
     * Size of the steps of a chain: moves of up to that many edges (3 to 5;
     * LKH's default is 5).
     */
    int move_type = 5;

    /** Maximum number of steps in a chain. */
    int maximum_depth = 50;

    /** See 'lin_kernighan_engine::EngineParameters::perturbation'. */
    std::string perturbation = "walks";

    /** See 'lin_kernighan_engine::EngineParameters::restricted_search'. */
    bool restricted_search = true;

    /** See 'lin_kernighan_engine::EngineParameters::non_sequential_moves'. */
    bool non_sequential_moves = true;

    /** Maximum number of trials or kicks ('-1': no limit, until the time limit). */
    int64_t maximum_number_of_trials = -1;

    /** Maximum length of each of the two segments swapped by a segment swap. */
    VertexId kick_segment_length = 50;

    /** Seed. */
    int seed = 0;


    virtual nlohmann::json to_json() const override
    {
        nlohmann::json json = Parameters::to_json();
        json.merge_patch({
                {"NumberOfCandidates", number_of_candidates},
                {"MoveType", move_type},
                {"MaximumDepth", maximum_depth},
                {"Perturbation", perturbation},
                {"RestrictedSearch", restricted_search},
                {"NonSequentialMoves", non_sequential_moves},
                {"MaximumNumberOfTrials", maximum_number_of_trials},
                {"KickSegmentLength", kick_segment_length},
                {"Seed", seed},
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
            << std::setw(width) << std::left << "Move type: " << move_type << std::endl
            << std::setw(width) << std::left << "Maximum depth: " << maximum_depth << std::endl
            << std::setw(width) << std::left << "Perturbation: " << perturbation << std::endl
            << std::setw(width) << std::left << "Restricted search: " << restricted_search << std::endl
            << std::setw(width) << std::left << "Non-sequential moves: " << non_sequential_moves << std::endl
            << std::setw(width) << std::left << "Maximum number of trials: " << maximum_number_of_trials << std::endl
            << std::setw(width) << std::left << "Kick segment length: " << kick_segment_length << std::endl
            << std::setw(width) << std::left << "Seed: " << seed << std::endl
            ;
    }
};

/**
 * LKH-style Lin-Kernighan for the TSP (see 'lin_kernighan/engine.hpp'):
 * nearest neighbor candidates and the engine, minimizing the length of the
 * tour.
 *
 * 'initial_solution', if provided, is the initial tour (a 'greedy' tour is
 * built otherwise). It must be a feasible solution of the same 'instance'.
 */
const Output lin_kernighan(
        const Instance& instance,
        const LinKernighanParameters& parameters = {},
        const Solution* initial_solution = nullptr);

template <typename Distances>
const Output lin_kernighan(
        const Distances& distances,
        const Instance& instance,
        const LinKernighanParameters& parameters = {},
        const Solution* initial_solution = nullptr);

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/** The TSP, for the Lin-Kernighan engine: the objective is the tour length. */
template <typename Distances>
struct LinKernighanTspProblem
{
    static constexpr bool evaluates_moves = false;
    static constexpr bool tracks_direction = false;
    static constexpr bool uses_weights = false;
    static constexpr bool has_dynamic_candidates = false;

    using Objective = Distance;

    const Distances& distances;

    const Instance& instance;

    AlgorithmFormatter& algorithm_formatter;

    inline VertexId number_of_vertices() const { return instance.number_of_vertices(); }

    inline Distance cost(VertexId vertex_id_1, VertexId vertex_id_2) const
    {
        return distances.distance(vertex_id_1, vertex_id_2);
    }

    inline bool candidate_edge(VertexId, VertexId) const { return true; }

    Objective evaluate_tour(const std::vector<VertexId>& tour) const
    {
        return lin_kernighan_engine::tour_cost(*this, tour);
    }

    std::vector<VertexId> trial_start_tour(
            lin_kernighan_engine::Engine<LinKernighanTspProblem>& engine,
            const std::vector<VertexId>& best_tour)
    {
        return lin_kernighan_engine::random_walk_tour(engine, best_tour);
    }

    void report(
            const std::vector<VertexId>& tour,
            const std::string& comment)
    {
        SolutionBuilder solution_builder(instance);
        for (VertexId vertex_id: tour)
            solution_builder.add_vertex(vertex_id);
        algorithm_formatter.update_solution(solution_builder.build(), comment);
    }
};

template <typename Distances>
const Output lin_kernighan(
        const Distances& distances,
        const Instance& instance,
        const LinKernighanParameters& parameters,
        const Solution* initial_solution)
{
    Output output(instance);
    AlgorithmFormatter algorithm_formatter(parameters, output);
    algorithm_formatter.start("Lin-Kernighan");
    algorithm_formatter.print_header();

    std::vector<VertexId> initial_tour;
    if (initial_solution != nullptr) {
        for (VertexId pos = 0; pos < initial_solution->number_of_vertices(); ++pos)
            initial_tour.push_back(initial_solution->vertex_id(pos));
    } else {
        GreedyParameters greedy_parameters;
        greedy_parameters.verbosity_level = 0;
        greedy_parameters.number_of_candidates = parameters.number_of_candidates;
        Output greedy_output = greedy(distances, instance, greedy_parameters);
        for (VertexId pos = 0; pos < greedy_output.solution.number_of_vertices(); ++pos)
            initial_tour.push_back(greedy_output.solution.vertex_id(pos));
    }
    CandidateLists candidates = nearest_neighbor_candidates(
            instance.distances(),
            parameters.number_of_candidates);

    lin_kernighan_engine::EngineParameters engine_parameters;
    engine_parameters.move_type = parameters.move_type;
    engine_parameters.maximum_depth = parameters.maximum_depth;
    engine_parameters.perturbation = parameters.perturbation;
    engine_parameters.restricted_search = parameters.restricted_search;
    engine_parameters.non_sequential_moves = parameters.non_sequential_moves;
    engine_parameters.maximum_number_of_trials = parameters.maximum_number_of_trials;
    engine_parameters.kick_segment_length = parameters.kick_segment_length;
    engine_parameters.seed = parameters.seed;
    engine_parameters.needs_to_end = [&parameters]() { return parameters.timer.needs_to_end(); };

    LinKernighanTspProblem<Distances> problem{distances, instance, algorithm_formatter};
    lin_kernighan_engine::run(problem, engine_parameters, std::move(candidates), initial_tour);

    algorithm_formatter.end();
    return output;
}

}
