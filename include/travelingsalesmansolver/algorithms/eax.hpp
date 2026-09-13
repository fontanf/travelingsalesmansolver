#pragma once

#include "travelingsalesmansolver/solution.hpp"
#include "travelingsalesmansolver/algorithm_formatter.hpp"

#include "travelingsalesmansolver/algorithms/eax/environment.hpp"

namespace travelingsalesmansolver
{

struct EaxParameters: Parameters
{
    /** Population size. */
    int population_size = 100;

    /** Number of children generated per generation. */
    int number_of_children = 30;

    /** Seed of the random number generator. */
    int seed = 0;
};

const Output eax(
        const Instance& instance,
        const EaxParameters& parameters = {});

template <typename Distances>
const Output eax(
        const Distances& distances,
        const Instance& instance,
        const EaxParameters& parameters = {});

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <typename Distances>
const Output eax(
        const Distances& distances,
        const Instance& instance,
        const EaxParameters& parameters)
{
    Output output(instance);
    AlgorithmFormatter algorithm_formatter(parameters, output);
    algorithm_formatter.start("EAX");
    algorithm_formatter.print_header();

    VertexId number_of_vertices = instance.number_of_vertices();

    // The vendored implementation relies on process-wide global state
    // ('tRand', 'tSort'), so this function is not re-entrant / thread-safe.
    // 'InitSort()' is never called anywhere in the original upstream
    // repository either (including its own 'main.cpp'), leaving 'tSort'
    // null until the GA reaches its "Block2" eset stage (only reached
    // after enough stalled generations, so small/quick-converging
    // instances never hit it) and dereferences it in 'TCross'.
    eax_ga::InitURandom(parameters.seed);
    eax_ga::InitSort();

    // 'distances' is looked up directly by 'TEvaluator<Distances>' (see
    // 'evaluator.hpp'), dispatched once here via the 'Distances' template
    // parameter -- no distance is ever copied into a separate matrix, the
    // same way every other algorithm in this library consumes distances.
    eax_ga::TEnvironment<Distances> environment(distances, number_of_vertices);
    environment.Npop = parameters.population_size;
    environment.Nch = parameters.number_of_children;
    environment.fTimeLimit = parameters.timer.remaining_time();
    environment.define();
    environment.doIt();

    std::vector<int> tour = environment.fEvaluator->getTour(environment.tBest);

    // 'tour[0]' is always city 0 (the traversal in 'getTour' starts there),
    // and the 'Solution' constructor already starts with vertex 0, so it
    // must be skipped here to avoid visiting it twice.
    Solution solution(instance);
    for (std::size_t i = 1; i < tour.size(); ++i)
        solution.add_vertex(distances, tour[i]);

    algorithm_formatter.update_solution(solution, "final solution");

    algorithm_formatter.end();
    return output;
}

}
