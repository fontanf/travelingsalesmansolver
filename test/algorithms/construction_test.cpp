#include "travelingsalesmansolver/algorithms/random_permutation.hpp"
#include "travelingsalesmansolver/algorithms/random_walk.hpp"
#include "travelingsalesmansolver/algorithms/two_opt.hpp"
#include "travelingsalesmansolver/algorithms/eax.hpp"
#include "travelingsalesmansolver/distances/distances_builder.hpp"

#include <gtest/gtest.h>

#include <random>

using namespace travelingsalesmansolver;

namespace
{

Instance random_euclidean_instance(
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
    return Instance(std::make_shared<const Distances>(distances_builder.build()));
}

}

TEST(RandomPermutation, Feasible)
{
    Instance instance = random_euclidean_instance(100, 0);
    RandomPermutationParameters parameters;
    parameters.verbosity_level = 0;
    std::mt19937_64 generator(0);
    for (int i = 0; i < 5; ++i) {
        Output output = random_permutation(instance, generator, parameters);
        EXPECT_TRUE(output.solution.feasible());
    }
}

TEST(RandomWalk, Feasible)
{
    Instance instance = random_euclidean_instance(200, 1);
    for (std::string candidates: {"alpha-nearness", "nearest-neighbor"}) {
        RandomWalkParameters parameters;
        parameters.verbosity_level = 0;
        parameters.candidates = candidates;
        std::mt19937_64 generator(0);
        Output output = random_walk(instance, generator, parameters);
        EXPECT_TRUE(output.solution.feasible());
    }
}

TEST(TwoOpt, LocalOptimum)
{
    // No 2-opt move with an added edge between candidates (both added edges
    // among the candidates of an end) improves the tour.
    Instance instance = random_euclidean_instance(300, 2);
    const DistancesEuc2D& distances = *instance.distances().distances_euc_2d();
    TwoOptParameters parameters;
    parameters.verbosity_level = 0;
    std::mt19937_64 generator(0);
    Output output = two_opt(instance, generator, parameters);
    ASSERT_TRUE(output.solution.feasible());
    CandidateLists candidates = nearest_neighbor_candidates(instance.distances(), parameters.number_of_candidates);
    VertexId n = instance.number_of_vertices();
    std::vector<VertexId> positions(n);
    for (VertexId pos = 0; pos < n; ++pos)
        positions[output.solution.vertex_id(pos)] = pos;
    auto next = [&](VertexId v) { return output.solution.vertex_id((positions[v] + 1) % n); };
    auto previous = [&](VertexId v) { return output.solution.vertex_id((positions[v] + n - 1) % n); };
    for (VertexId t1 = 0; t1 < n; ++t1) {
        for (int side = 0; side < 2; ++side) {
            VertexId t2 = (side == 0)? next(t1): previous(t1);
            for (VertexId t4: candidates[t1]) {
                VertexId t3 = (side == 0)? next(t4): previous(t4);
                if (t4 == t2 || t3 == t1)
                    continue;
                Distance gain = distances.distance(t1, t2) + distances.distance(t3, t4)
                    - distances.distance(t1, t4) - distances.distance(t2, t3);
                EXPECT_LE(gain, 0) << t1 << " " << t2 << " " << t3 << " " << t4;
            }
        }
    }
}

TEST(Eax, InitialPopulations)
{
    Instance instance = random_euclidean_instance(150, 3);
    for (std::string initial_tours: {"random-permutation", "random-walk"}) {
        for (std::string initial_local_search: {"2-opt", "lin-kernighan", "none"}) {
            EaxParameters parameters;
            parameters.verbosity_level = 0;
            parameters.population_size = 10;
            parameters.number_of_children = 10;
            parameters.initial_tours = initial_tours;
            parameters.initial_local_search = initial_local_search;
            parameters.timer.set_time_limit(5);
            std::mt19937_64 generator(0);
            Output output = eax(instance, generator, parameters);
            EXPECT_TRUE(output.solution.feasible());
        }
    }
}
