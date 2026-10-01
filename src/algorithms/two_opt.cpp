#include "travelingsalesmansolver/algorithms/two_opt.hpp"

using namespace travelingsalesmansolver;

const Output travelingsalesmansolver::two_opt(
        const Instance& instance,
        std::mt19937_64& generator,
        const TwoOptParameters& parameters,
        const Solution* initial_solution,
        const CandidateLists* candidates)
{
    if (initial_solution != nullptr
            && (&initial_solution->instance() != &instance
                || !initial_solution->feasible())) {
        throw std::invalid_argument(
                "travelingsalesmansolver::two_opt: "
                "the initial solution must be a feasible solution of the instance.");
    }
    if (candidates != nullptr
            && (VertexId)candidates->size() != instance.number_of_vertices()) {
        throw std::invalid_argument(
                "travelingsalesmansolver::two_opt: "
                "wrong number of candidate lists; "
                "candidates->size(): " + std::to_string(candidates->size()) + "; "
                "instance.number_of_vertices(): " + std::to_string(instance.number_of_vertices()) + ".");
    }
    return FUNCTION_WITH_DISTANCES(
            two_opt,
            instance.distances(),
            instance,
            generator,
            parameters,
            initial_solution,
            candidates);
}
