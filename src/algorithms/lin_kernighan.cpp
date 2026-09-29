#include "travelingsalesmansolver/algorithms/lin_kernighan.hpp"

using namespace travelingsalesmansolver;

const Output travelingsalesmansolver::lin_kernighan(
        const Instance& instance,
        const LinKernighanParameters& parameters)
{
    if (parameters.move_type < 2
            || parameters.move_type > SequentialMovePatterns::maximum_number_of_removed_edges) {
        throw std::invalid_argument(
                "travelingsalesmansolver::lin_kernighan: "
                "the move type must be between 2 and "
                + std::to_string(SequentialMovePatterns::maximum_number_of_removed_edges)
                + "; move type: " + std::to_string(parameters.move_type) + ".");
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
            parameters);
}
