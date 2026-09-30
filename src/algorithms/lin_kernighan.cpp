#include "travelingsalesmansolver/algorithms/lin_kernighan.hpp"

using namespace travelingsalesmansolver;

const Output travelingsalesmansolver::lin_kernighan(
        const Instance& instance,
        const LinKernighanParameters& parameters,
        const Solution* initial_solution,
        const CandidateLists* candidates,
        const VertexPenalties* penalties)
{
    if (initial_solution != nullptr
            && (&initial_solution->instance() != &instance
                || !initial_solution->feasible())) {
        throw std::invalid_argument(
                "travelingsalesmansolver::lin_kernighan: "
                "the initial solution must be a feasible solution of the instance.");
    }
    if (parameters.move_type < 3
            || parameters.move_type > SequentialMovePatterns::maximum_number_of_removed_edges) {
        throw std::invalid_argument(
                "travelingsalesmansolver::lin_kernighan: "
                "the move type must be between 3 and "
                + std::to_string(SequentialMovePatterns::maximum_number_of_removed_edges)
                + "; move type: " + std::to_string(parameters.move_type) + ".");
    }
    if (candidates != nullptr
            && (VertexId)candidates->size() != instance.number_of_vertices()) {
        throw std::invalid_argument(
                "travelingsalesmansolver::lin_kernighan: "
                "wrong number of candidate lists; "
                "candidates->size(): " + std::to_string(candidates->size()) + "; "
                "instance.number_of_vertices(): " + std::to_string(instance.number_of_vertices()) + ".");
    }

    if (parameters.candidates != "alpha-nearness"
            && parameters.candidates != "nearest-neighbor") {
        throw std::invalid_argument(
                "travelingsalesmansolver::lin_kernighan: "
                "unknown candidates \"" + parameters.candidates + "\".");
    }

    // Default candidates, computed before the distance matrix (the
    // alpha-nearness graph uses the coordinates).
    std::unique_ptr<CandidateLists> default_candidates;
    std::unique_ptr<VertexPenalties> default_penalties;
    if (candidates == nullptr) {
        if (parameters.candidates == "alpha-nearness") {
            AlphaNearnessParameters alpha_nearness_parameters;
            alpha_nearness_parameters.number_of_candidates = parameters.number_of_candidates;
            if (parameters.ascent_initial_period >= 0) {
                alpha_nearness_parameters.initial_period = parameters.ascent_initial_period;
            } else if (instance.number_of_vertices() >= 10000) {
                alpha_nearness_parameters.initial_period = 100;
            }
            AlphaNearnessOutput alpha_nearness_output = alpha_nearness_candidates(
                    instance.distances(),
                    alpha_nearness_parameters);
            default_candidates = std::make_unique<CandidateLists>(
                    std::move(alpha_nearness_output.candidates));
            if (parameters.penalized_costs && penalties == nullptr) {
                default_penalties = std::make_unique<VertexPenalties>();
                default_penalties->pi = std::move(alpha_nearness_output.pi);
                default_penalties->precision = alpha_nearness_parameters.precision;
                penalties = default_penalties.get();
            }
        } else {
            default_candidates = std::make_unique<CandidateLists>(
                    nearest_neighbor_candidates(instance.distances(), parameters.number_of_candidates));
        }
        candidates = default_candidates.get();
    }

    // Like LKH, precompute the distance matrix when the instance is small
    // enough, to avoid recomputing distances in the inner loop.
    if (instance.number_of_vertices() <= 16000) {
        instance.distances().compute_distances_explicit();
    } else if (instance.number_of_vertices() <= 20000) {
        instance.distances().compute_distances_explicit_triangle();
    }

    return FUNCTION_WITH_DISTANCES(
            lin_kernighan,
            instance.distances(),
            instance,
            parameters,
            initial_solution,
            candidates,
            penalties);
}
