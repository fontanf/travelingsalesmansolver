#pragma once

#include "travelingsalesmansolver/algorithm_formatter.hpp"
#include "travelingsalesmansolver/solution_builder.hpp"
#include "travelingsalesmansolver/candidates/nearest_neighbor.hpp"

#include <algorithm>
#include <array>
#include <numeric>
#include <tuple>
#include <vector>

namespace travelingsalesmansolver
{

/**
 * Greedy-edge tour (Bentley's multiple fragment heuristic): candidate edges
 * are added by increasing distance, unless they would give a vertex a third
 * edge or close a cycle too early. The remaining fragments (paths, or
 * isolated vertices) are then chained, each time to the nearest free end of
 * another fragment.
 *
 * Returns the tour as a sequence of vertices.
 */
template <typename Distances>
std::vector<VertexId> greedy_edge_tour(
        const Distances& distances,
        const CandidateLists& candidates)
{
    VertexId number_of_vertices = candidates.size();
    if (number_of_vertices <= 3) {
        std::vector<VertexId> tour(number_of_vertices);
        std::iota(tour.begin(), tour.end(), 0);
        return tour;
    }

    // Candidate edges, each once, by increasing distance.
    std::vector<std::tuple<Distance, VertexId, VertexId>> edges;
    for (VertexId vertex_id_1 = 0; vertex_id_1 < number_of_vertices; ++vertex_id_1) {
        for (VertexId vertex_id_2: candidates[vertex_id_1]) {
            if (vertex_id_1 < vertex_id_2) {
                edges.emplace_back(distances.distance(vertex_id_1, vertex_id_2), vertex_id_1, vertex_id_2);
            } else {
                edges.emplace_back(distances.distance(vertex_id_2, vertex_id_1), vertex_id_2, vertex_id_1);
            }
        }
    }
    std::sort(edges.begin(), edges.end());
    edges.erase(std::unique(edges.begin(), edges.end()), edges.end());

    // Fragments: the neighbors of each vertex ('-1' if none), and a
    // union-find structure to detect cycles.
    std::vector<std::array<VertexId, 2>> neighbors(number_of_vertices, {-1, -1});
    std::vector<VertexId> parents(number_of_vertices);
    std::iota(parents.begin(), parents.end(), 0);
    auto find = [&parents](VertexId vertex_id)
    {
        while (parents[vertex_id] != vertex_id) {
            parents[vertex_id] = parents[parents[vertex_id]];
            vertex_id = parents[vertex_id];
        }
        return vertex_id;
    };
    auto degree = [&neighbors](VertexId vertex_id)
    {
        return (neighbors[vertex_id][0] != -1) + (neighbors[vertex_id][1] != -1);
    };
    auto link = [&neighbors](VertexId vertex_id_1, VertexId vertex_id_2)
    {
        neighbors[vertex_id_1][(neighbors[vertex_id_1][0] == -1)? 0: 1] = vertex_id_2;
        neighbors[vertex_id_2][(neighbors[vertex_id_2][0] == -1)? 0: 1] = vertex_id_1;
    };
    VertexId number_of_edges = 0;
    for (const auto& edge: edges) {
        VertexId vertex_id_1 = std::get<1>(edge);
        VertexId vertex_id_2 = std::get<2>(edge);
        if (degree(vertex_id_1) == 2 || degree(vertex_id_2) == 2)
            continue;
        VertexId root_1 = find(vertex_id_1);
        VertexId root_2 = find(vertex_id_2);
        if (root_1 == root_2)
            continue;
        parents[root_1] = root_2;
        link(vertex_id_1, vertex_id_2);
        number_of_edges++;
        if (number_of_edges == number_of_vertices - 1)
            break;
    }

    // Chain the fragments. The ends of each fragment: vertices of degree
    // 0 (both ends) or 1.
    auto other_end = [&neighbors](VertexId end_vertex_id)
    {
        VertexId previous_vertex_id = -1;
        VertexId vertex_id = end_vertex_id;
        for (;;) {
            VertexId next_vertex_id = (neighbors[vertex_id][0] != previous_vertex_id)?
                neighbors[vertex_id][0]:
                neighbors[vertex_id][1];
            if (next_vertex_id == -1)
                return vertex_id;
            previous_vertex_id = vertex_id;
            vertex_id = next_vertex_id;
        }
    };
    std::vector<VertexId> free_ends;
    for (VertexId vertex_id = 0; vertex_id < number_of_vertices; ++vertex_id)
        if (degree(vertex_id) < 2)
            free_ends.push_back(vertex_id);
    if (!free_ends.empty()) {
        std::vector<bool> used(number_of_vertices, false);
        VertexId first_end = free_ends.front();
        VertexId current_end = other_end(first_end);
        used[first_end] = used[current_end] = true;
        for (;;) {
            VertexId best_end = -1;
            Distance best_distance = 0;
            for (VertexId end_vertex_id: free_ends) {
                if (used[end_vertex_id])
                    continue;
                Distance distance = distances.distance(current_end, end_vertex_id);
                if (best_end == -1 || distance < best_distance) {
                    best_end = end_vertex_id;
                    best_distance = distance;
                }
            }
            if (best_end == -1)
                break;
            VertexId next_end = other_end(best_end);
            used[best_end] = used[next_end] = true;
            link(current_end, best_end);
            current_end = next_end;
        }
        link(current_end, first_end);
    }

    // Read the tour.
    std::vector<VertexId> tour;
    tour.reserve(number_of_vertices);
    VertexId previous_vertex_id = neighbors[0][1];
    VertexId vertex_id = 0;
    for (VertexId pos = 0; pos < number_of_vertices; ++pos) {
        tour.push_back(vertex_id);
        VertexId next_vertex_id = (neighbors[vertex_id][0] != previous_vertex_id)?
            neighbors[vertex_id][0]:
            neighbors[vertex_id][1];
        previous_vertex_id = vertex_id;
        vertex_id = next_vertex_id;
    }
    return tour;
}


struct GreedyEdgeParameters: Parameters
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
const Output greedy_edge(
        const Instance& instance,
        const GreedyEdgeParameters& parameters = {},
        const CandidateLists* candidates = nullptr);

template <typename Distances>
const Output greedy_edge(
        const Distances& distances,
        const Instance& instance,
        const GreedyEdgeParameters& parameters = {},
        const CandidateLists* candidates = nullptr)
{
    Output output(instance);
    AlgorithmFormatter algorithm_formatter(parameters, output);
    algorithm_formatter.start("Greedy edge");
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
