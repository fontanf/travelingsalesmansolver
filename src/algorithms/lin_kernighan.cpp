#include "travelingsalesmansolver/algorithms/lin_kernighan.hpp"

using namespace travelingsalesmansolver;

const Output travelingsalesmansolver::lin_kernighan(
        const Instance& instance,
        const LinKernighanParameters& parameters,
        const Solution* initial_solution,
        const CandidateLists* candidates)
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
            candidates);
}
