#include "travelingsalesmansolver/algorithms/random_permutation.hpp"

#include "travelingsalesmansolver/solution_builder.hpp"

#include <algorithm>
#include <numeric>
#include <random>

using namespace travelingsalesmansolver;

const Output travelingsalesmansolver::random_permutation(
        const Instance& instance,
        std::mt19937_64& generator,
        const RandomPermutationParameters& parameters)
{
    Output output(instance);
    AlgorithmFormatter algorithm_formatter(parameters, output);
    algorithm_formatter.start("Random permutation");
    algorithm_formatter.print_header();

    std::vector<VertexId> tour(instance.number_of_vertices());
    std::iota(tour.begin(), tour.end(), 0);
    std::shuffle(tour.begin(), tour.end(), generator);

    SolutionBuilder solution_builder(instance);
    for (VertexId vertex_id: tour)
        solution_builder.add_vertex(vertex_id);
    algorithm_formatter.update_solution(solution_builder.build(), "");

    algorithm_formatter.end();
    return output;
}
