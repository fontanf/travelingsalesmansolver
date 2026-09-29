#include "travelingsalesmansolver/candidates/alpha_nearness.hpp"
#include "travelingsalesmansolver/distances/distances_builder.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <numeric>
#include <random>
#include <set>

using namespace travelingsalesmansolver;

namespace
{

Distances random_euclidean_distances(
        VertexId number_of_vertices,
        uint64_t seed)
{
    std::mt19937_64 generator(seed);
    std::uniform_real_distribution<double> coordinate(0, 1000);
    DistancesEuc2DBuilder euc_2d_builder;
    euc_2d_builder.set_number_of_vertices(number_of_vertices);
    for (VertexId vertex_id = 0; vertex_id < number_of_vertices; ++vertex_id)
        euc_2d_builder.set_coordinates(vertex_id, coordinate(generator), coordinate(generator));
    DistancesBuilder distances_builder;
    distances_builder.set_number_of_vertices(number_of_vertices);
    distances_builder.set_distances_euc_2d(euc_2d_builder.build());
    return distances_builder.build();
}

/** Cost of a minimum spanning tree of the complete graph (Prim). */
Distance minimum_spanning_tree_cost(const std::vector<std::vector<Distance>>& costs)
{
    VertexId n = costs.size();
    std::vector<Distance> keys(n, std::numeric_limits<Distance>::max());
    std::vector<bool> in_tree(n, false);
    keys[0] = 0;
    Distance cost = 0;
    for (VertexId iteration = 0; iteration < n; ++iteration) {
        VertexId best = -1;
        for (VertexId vertex_id = 0; vertex_id < n; ++vertex_id)
            if (!in_tree[vertex_id] && (best == -1 || keys[vertex_id] < keys[best]))
                best = vertex_id;
        in_tree[best] = true;
        cost += keys[best];
        for (VertexId vertex_id = 0; vertex_id < n; ++vertex_id)
            if (!in_tree[vertex_id])
                keys[vertex_id] = std::min(keys[vertex_id], costs[best][vertex_id]);
    }
    return cost;
}

void check_candidates(
        const AlphaNearnessOutput& output,
        VertexId number_of_vertices,
        VertexId number_of_candidates)
{
    ASSERT_EQ((VertexId)output.candidates.size(), number_of_vertices);
    for (VertexId vertex_id = 0; vertex_id < number_of_vertices; ++vertex_id) {
        const auto& candidates = output.candidates[vertex_id];
        EXPECT_EQ((VertexId)candidates.size(), number_of_candidates);
        std::set<VertexId> distinct(candidates.begin(), candidates.end());
        EXPECT_EQ(distinct.size(), candidates.size());
        EXPECT_EQ(distinct.count(vertex_id), 0);
        EXPECT_TRUE(std::is_sorted(output.alphas[vertex_id].begin(), output.alphas[vertex_id].end()));
        EXPECT_GE(output.alphas[vertex_id].front(), 0);
    }
}

}

TEST(AlphaNearness, AlphaValuesMatchTheirDefinition)
{
    // The alpha value of an edge not incident to the special vertex is the
    // increase of the cost of a minimum spanning tree forced to contain it.
    VertexId n = 60;
    Distances distances = random_euclidean_distances(n, 0);
    AlphaNearnessParameters parameters;
    parameters.number_of_nearest_neighbors = -1;
    parameters.number_of_candidates = 8;
    AlphaNearnessOutput output = alpha_nearness_candidates(distances, parameters);
    check_candidates(output, n, 8);

    const DistancesEuc2D& euc_2d = *distances.distances_euc_2d();
    std::vector<std::vector<Distance>> costs(n, std::vector<Distance>(n));
    for (VertexId i = 0; i < n; ++i)
        for (VertexId j = 0; j < n; ++j)
            costs[i][j] = parameters.precision * euc_2d.distance(i, j) + output.pi[i] + output.pi[j];
    Distance mst_cost = minimum_spanning_tree_cost(costs);

    // Pairs whose alpha differs from the definition: they must all be
    // incident to a single (special) vertex.
    std::vector<std::pair<VertexId, VertexId>> mismatches;
    for (VertexId i = 0; i < n; ++i) {
        for (VertexId pos = 0; pos < (VertexId)output.candidates[i].size(); ++pos) {
            VertexId j = output.candidates[i][pos];
            // Force (i, j): contract it.
            std::vector<std::vector<Distance>> forced = costs;
            forced[i][j] = forced[j][i] = std::numeric_limits<Distance>::min() / 4;
            Distance forced_cost = minimum_spanning_tree_cost(forced)
                - forced[i][j] + costs[i][j];
            if (forced_cost - mst_cost != output.alphas[i][pos])
                mismatches.push_back({i, j});
        }
    }
    if (!mismatches.empty()) {
        VertexId special = -1;
        for (VertexId candidate: {mismatches[0].first, mismatches[0].second}) {
            bool all = true;
            for (const auto& mismatch: mismatches)
                if (mismatch.first != candidate && mismatch.second != candidate)
                    all = false;
            if (all)
                special = candidate;
        }
        EXPECT_NE(special, -1);
    }
}

TEST(AlphaNearness, LowerBound)
{
    for (uint64_t seed = 0; seed < 5; ++seed) {
        VertexId n = 9;
        Distances distances = random_euclidean_distances(n, seed);
        const DistancesEuc2D& euc_2d = *distances.distances_euc_2d();
        std::vector<VertexId> tour(n);
        std::iota(tour.begin(), tour.end(), 0);
        Distance optimal = std::numeric_limits<Distance>::max();
        do {
            Distance length = 0;
            for (VertexId pos = 0; pos < n; ++pos)
                length += euc_2d.distance(tour[pos], tour[(pos + 1) % n]);
            optimal = std::min(optimal, length);
        } while (std::next_permutation(tour.begin() + 1, tour.end()));

        AlphaNearnessParameters parameters;
        parameters.number_of_nearest_neighbors = -1;
        AlphaNearnessOutput output = alpha_nearness_candidates(distances, parameters);
        check_candidates(output, n, 5);
        EXPECT_LE(output.lower_bound, optimal + 1e-9);
        EXPECT_GE(output.lower_bound, 0.9 * optimal);
        if (output.tour_found)
            EXPECT_NEAR(output.lower_bound, optimal, 1e-9);
    }
}

TEST(AlphaNearness, SparseGraph)
{
    VertexId n = 500;
    Distances distances = random_euclidean_distances(n, 1);
    AlphaNearnessParameters parameters_complete;
    parameters_complete.number_of_nearest_neighbors = -1;
    AlphaNearnessOutput output_complete = alpha_nearness_candidates(distances, parameters_complete);
    check_candidates(output_complete, n, 5);
    AlphaNearnessOutput output_sparse = alpha_nearness_candidates(distances);
    check_candidates(output_sparse, n, 5);
    EXPECT_NEAR(output_sparse.lower_bound, output_complete.lower_bound, 1e-3 * output_complete.lower_bound);

    // Most candidates are the same.
    VertexId common = 0;
    for (VertexId vertex_id = 0; vertex_id < n; ++vertex_id) {
        std::set<VertexId> complete(
                output_complete.candidates[vertex_id].begin(),
                output_complete.candidates[vertex_id].end());
        for (VertexId other_vertex_id: output_sparse.candidates[vertex_id])
            common += complete.count(other_vertex_id);
    }
    EXPECT_GE(common, 0.9 * 5 * n);
}

TEST(AlphaNearness, DisconnectedNearestNeighborGraph)
{
    // Two clusters, far apart: the graph of the 2 nearest neighbors isn't
    // connected.
    VertexId n = 20;
    DistancesExplicitBuilder explicit_builder;
    explicit_builder.set_number_of_vertices(n);
    for (VertexId i = 0; i < n; ++i) {
        for (VertexId j = 0; j < n; ++j) {
            Distance distance = (i == j)? 0: (i < n / 2) == (j < n / 2)? 1 + (i + j) % 3: 1000;
            explicit_builder.set_distance(i, j, distance);
        }
    }
    DistancesBuilder distances_builder;
    distances_builder.set_number_of_vertices(n);
    distances_builder.set_distances_explicit(explicit_builder.build());
    Distances distances = distances_builder.build();

    AlphaNearnessParameters parameters;
    parameters.number_of_nearest_neighbors = 2;
    parameters.number_of_candidates = 3;
    AlphaNearnessOutput output = alpha_nearness_candidates(distances, parameters);
    ASSERT_EQ((VertexId)output.candidates.size(), n);
    for (VertexId vertex_id = 0; vertex_id < n; ++vertex_id)
        EXPECT_FALSE(output.candidates[vertex_id].empty());
    // A tour crosses between the clusters twice.
    EXPECT_GE(output.lower_bound, 2000);
}
