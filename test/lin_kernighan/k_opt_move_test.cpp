#include "travelingsalesmansolver/lin_kernighan/k_opt_move.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <numeric>
#include <random>
#include <set>

using namespace travelingsalesmansolver;

namespace
{

using Edge = std::pair<VertexId, VertexId>;

Edge make_edge(VertexId vertex_id_1, VertexId vertex_id_2)
{
    return {(std::min)(vertex_id_1, vertex_id_2), (std::max)(vertex_id_1, vertex_id_2)};
}

std::multiset<Edge> tour_edges(const TwoLevelList& tour)
{
    std::multiset<Edge> edges;
    for (VertexId vertex_id = 0; vertex_id < tour.number_of_vertices(); ++vertex_id)
        edges.insert(make_edge(vertex_id, tour.next(vertex_id)));
    return edges;
}

/** Edges after the move, or an empty set if a removed edge isn't in the tour. */
std::multiset<Edge> edges_after_move(
        const TwoLevelList& tour,
        const std::vector<VertexId>& t)
{
    std::multiset<Edge> edges = tour_edges(tour);
    int k = t.size() / 2;
    for (int i = 0; i < k; ++i) {
        auto it = edges.find(make_edge(t[2 * i], t[2 * i + 1]));
        if (it == edges.end())
            return {};
        edges.erase(it);
    }
    for (int i = 0; i < k; ++i)
        edges.insert(make_edge(t[2 * i + 1], t[(2 * i + 2) % (2 * k)]));
    return edges;
}

/** Is the edge set a single Hamiltonian cycle? */
bool is_hamiltonian_cycle(
        const std::multiset<Edge>& edges,
        VertexId number_of_vertices)
{
    if ((VertexId)edges.size() != number_of_vertices)
        return false;
    std::vector<std::vector<VertexId>> neighbors(number_of_vertices);
    for (const Edge& edge: edges) {
        if (edge.first == edge.second)
            return false;
        neighbors[edge.first].push_back(edge.second);
        neighbors[edge.second].push_back(edge.first);
    }
    for (const auto& vertex_neighbors: neighbors)
        if (vertex_neighbors.size() != 2 || vertex_neighbors[0] == vertex_neighbors[1])
            return false;
    VertexId previous_vertex_id = -1;
    VertexId vertex_id = 0;
    for (VertexId pos = 0; pos < number_of_vertices; ++pos) {
        VertexId next_vertex_id = (neighbors[vertex_id][0] != previous_vertex_id)?
            neighbors[vertex_id][0]:
            neighbors[vertex_id][1];
        previous_vertex_id = vertex_id;
        vertex_id = next_vertex_id;
        if (vertex_id == 0 && pos < number_of_vertices - 1)
            return false;
    }
    return vertex_id == 0;
}

}

TEST(SequentialMove, RandomMoves)
{
    std::mt19937_64 generator(0);
    SequentialMove move;
    int number_of_valid_moves = 0;
    for (VertexId n: {8, 12, 30, 200}) {
        std::uniform_int_distribution<VertexId> distribution(0, n - 1);
        for (int k = 2; k <= 5; ++k) {
            for (int trial = 0; trial < 3000; ++trial) {
                std::vector<VertexId> initial_tour(n);
                std::iota(initial_tour.begin(), initial_tour.end(), 0);
                std::shuffle(initial_tour.begin(), initial_tour.end(), generator);
                TwoLevelList tour(initial_tour, 3);
                // Random removed tour edges; on small tours, close to each
                // other so that valid moves aren't too rare.
                std::vector<VertexId> t(2 * k);
                std::set<Edge> removed_edges;
                bool distinct = true;
                for (int i = 0; i < k; ++i) {
                    t[2 * i] = distribution(generator);
                    t[2 * i + 1] = (generator() % 2 == 0)?
                        tour.next(t[2 * i]):
                        tour.previous(t[2 * i]);
                    if (!removed_edges.insert(make_edge(t[2 * i], t[2 * i + 1])).second)
                        distinct = false;
                }
                if (!distinct)
                    continue;
                std::multiset<Edge> expected_edges = edges_after_move(tour, t);
                bool expected_valid = is_hamiltonian_cycle(expected_edges, n);
                bool valid = move.analyze(tour, t.data(), k);
                ASSERT_EQ(valid, expected_valid) << "n " << n << " k " << k << " trial " << trial;
                if (!valid)
                    continue;
                number_of_valid_moves++;
                move.apply(tour, [&tour](VertexId vertex_id_1, VertexId vertex_id_2)
                {
                    tour.reverse(vertex_id_1, vertex_id_2);
                });
                ASSERT_EQ(tour_edges(tour), expected_edges) << "n " << n << " k " << k << " trial " << trial;
            }
        }
    }
    EXPECT_GT(number_of_valid_moves, 1000);
}
