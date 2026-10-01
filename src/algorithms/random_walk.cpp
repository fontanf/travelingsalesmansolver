#include "travelingsalesmansolver/algorithms/random_walk.hpp"

#include "travelingsalesmansolver/solution_builder.hpp"
#include "travelingsalesmansolver/candidates/nearest_neighbor.hpp"
#include "travelingsalesmansolver/candidates/alpha_nearness.hpp"

#include <numeric>
#include <random>

using namespace travelingsalesmansolver;

const Output travelingsalesmansolver::random_walk(
        const Instance& instance,
        std::mt19937_64& generator,
        const RandomWalkParameters& parameters,
        const CandidateLists* candidates)
{
    VertexId number_of_vertices = instance.number_of_vertices();
    if (candidates != nullptr && (VertexId)candidates->size() != number_of_vertices) {
        throw std::invalid_argument(
                "travelingsalesmansolver::random_walk: "
                "wrong number of candidate lists; "
                "candidates->size(): " + std::to_string(candidates->size()) + "; "
                "instance.number_of_vertices(): " + std::to_string(number_of_vertices) + ".");
    }

    Output output(instance);
    AlgorithmFormatter algorithm_formatter(parameters, output);
    algorithm_formatter.start("Random walk");
    algorithm_formatter.print_header();

    CandidateLists default_candidates;
    if (candidates == nullptr) {
        if (parameters.candidates == "alpha-nearness") {
            AlphaNearnessParameters alpha_nearness_parameters;
            alpha_nearness_parameters.number_of_candidates = parameters.number_of_candidates;
            if (number_of_vertices >= 10000)
                alpha_nearness_parameters.initial_period = 100;
            default_candidates = alpha_nearness_candidates(
                    instance.distances(),
                    alpha_nearness_parameters).candidates;
        } else if (parameters.candidates == "nearest-neighbor") {
            default_candidates = nearest_neighbor_candidates(
                    instance.distances(),
                    parameters.number_of_candidates);
        } else {
            throw std::invalid_argument(
                    "travelingsalesmansolver::random_walk: "
                    "unknown candidates \"" + parameters.candidates + "\".");
        }
        candidates = &default_candidates;
    }

    // Unvisited vertices, and the position of each vertex among them ('-1'
    // once visited).
    std::vector<VertexId> unvisited(number_of_vertices);
    std::iota(unvisited.begin(), unvisited.end(), 0);
    std::vector<VertexId> positions(number_of_vertices);
    std::iota(positions.begin(), positions.end(), 0);
    SolutionBuilder solution_builder(instance);
    std::vector<VertexId> options;
    VertexId vertex_id = (number_of_vertices > 0)?
        std::uniform_int_distribution<VertexId>(0, number_of_vertices - 1)(generator):
        -1;
    while (vertex_id != -1) {
        solution_builder.add_vertex(vertex_id);
        VertexId pos = positions[vertex_id];
        VertexId last_vertex_id = unvisited.back();
        unvisited[pos] = last_vertex_id;
        positions[last_vertex_id] = pos;
        unvisited.pop_back();
        positions[vertex_id] = -1;
        if (unvisited.empty())
            break;
        options.clear();
        for (VertexId candidate_id: (*candidates)[vertex_id])
            if (positions[candidate_id] != -1)
                options.push_back(candidate_id);
        if (!options.empty()) {
            vertex_id = options[std::uniform_int_distribution<size_t>(0, options.size() - 1)(generator)];
        } else {
            vertex_id = unvisited[std::uniform_int_distribution<size_t>(0, unvisited.size() - 1)(generator)];
        }
    }
    algorithm_formatter.update_solution(solution_builder.build(), "");

    algorithm_formatter.end();
    return output;
}
