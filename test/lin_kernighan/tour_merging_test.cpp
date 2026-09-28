#include "travelingsalesmansolver/lin_kernighan/tour_merging.hpp"
#include "travelingsalesmansolver/distances/distances_builder.hpp"

#include <gtest/gtest.h>

#include <numeric>
#include <random>

using namespace travelingsalesmansolver;

namespace
{

Distance tour_length(const DistancesEuc2D& distances, const std::vector<VertexId>& tour)
{
    Distance length = 0;
    for (size_t pos = 0; pos < tour.size(); ++pos)
        length += distances.distance(tour[pos], tour[(pos + 1) % tour.size()]);
    return length;
}

/** Improve a tour with 2-opt moves until none improves it. */
void two_opt(const DistancesEuc2D& distances, std::vector<VertexId>& tour)
{
    VertexId n = tour.size();
    bool improved = true;
    while (improved) {
        improved = false;
        for (VertexId i = 0; i < n - 1; ++i) {
            for (VertexId j = i + 2; j < n; ++j) {
                VertexId a = tour[i], b = tour[i + 1], c = tour[j], d = tour[(j + 1) % n];
                if (d == a)
                    continue;
                if (distances.distance(a, c) + distances.distance(b, d)
                        < distances.distance(a, b) + distances.distance(c, d)) {
                    std::reverse(tour.begin() + i + 1, tour.begin() + j + 1);
                    improved = true;
                }
            }
        }
    }
}

}

TEST(TourMerging, PartitionCrossover)
{
    std::mt19937_64 generator(0);
    std::uniform_real_distribution<double> coordinate(0, 1000);
    int number_of_improvements = 0;
    for (int instance = 0; instance < 20; ++instance) {
        VertexId n = 200;
        DistancesEuc2DBuilder builder;
        builder.set_number_of_vertices(n);
        for (VertexId vertex_id = 0; vertex_id < n; ++vertex_id)
            builder.set_coordinates(vertex_id, coordinate(generator), coordinate(generator));
        DistancesEuc2D distances = builder.build();
        for (bool optimized: {false, true}) {
            std::vector<VertexId> tour_1(n), tour_2(n);
            std::iota(tour_1.begin(), tour_1.end(), 0);
            std::iota(tour_2.begin(), tour_2.end(), 0);
            std::shuffle(tour_1.begin(), tour_1.end(), generator);
            std::shuffle(tour_2.begin(), tour_2.end(), generator);
            if (optimized) {
                // Two similar local optima, as in a trial loop: 'tour_2' is
                // 'tour_1' perturbed by a few reversals, then optimized.
                two_opt(distances, tour_1);
                tour_2 = tour_1;
                std::uniform_int_distribution<VertexId> position(0, n - 1);
                for (int k = 0; k < 5; ++k) {
                    VertexId i = position(generator), j = position(generator);
                    std::reverse(tour_2.begin() + (std::min)(i, j), tour_2.begin() + (std::max)(i, j) + 1);
                }
                two_opt(distances, tour_2);
            }
            std::vector<VertexId> child = partition_crossover(distances, tour_1, tour_2);
            // A permutation, no longer than 'tour_1'.
            std::vector<VertexId> sorted = child;
            std::sort(sorted.begin(), sorted.end());
            std::vector<VertexId> identity(n);
            std::iota(identity.begin(), identity.end(), 0);
            ASSERT_EQ(sorted, identity);
            Distance length_1 = tour_length(distances, tour_1);
            Distance child_length = tour_length(distances, child);
            EXPECT_LE(child_length, length_1);
            if (child_length < length_1 && child_length < tour_length(distances, tour_2))
                number_of_improvements++;
        }
    }
    // On 2-opt local optima, merging finds tours better than both parents.
    EXPECT_GT(number_of_improvements, 0);
}

TEST(TourMerging, IterativePartialTranscription)
{
    std::mt19937_64 generator(1);
    std::uniform_real_distribution<double> coordinate(0, 1000);
    int number_of_improvements = 0;
    for (int instance = 0; instance < 20; ++instance) {
        VertexId n = 200;
        DistancesEuc2DBuilder builder;
        builder.set_number_of_vertices(n);
        for (VertexId vertex_id = 0; vertex_id < n; ++vertex_id)
            builder.set_coordinates(vertex_id, coordinate(generator), coordinate(generator));
        DistancesEuc2D distances = builder.build();
        // Two independent 2-opt local optima.
        std::vector<VertexId> tour_1(n), tour_2(n);
        std::iota(tour_1.begin(), tour_1.end(), 0);
        std::iota(tour_2.begin(), tour_2.end(), 0);
        std::shuffle(tour_1.begin(), tour_1.end(), generator);
        std::shuffle(tour_2.begin(), tour_2.end(), generator);
        two_opt(distances, tour_1);
        two_opt(distances, tour_2);
        std::vector<VertexId> child = iterative_partial_transcription(distances, tour_1, tour_2, generator);
        std::vector<VertexId> sorted = child;
        std::sort(sorted.begin(), sorted.end());
        std::vector<VertexId> identity(n);
        std::iota(identity.begin(), identity.end(), 0);
        ASSERT_EQ(sorted, identity);
        Distance length_1 = tour_length(distances, tour_1);
        Distance length_2 = tour_length(distances, tour_2);
        Distance child_length = tour_length(distances, child);
        EXPECT_LE(child_length, (std::min)(length_1, length_2));
        if (child_length < (std::min)(length_1, length_2))
            number_of_improvements++;
    }
    EXPECT_GT(number_of_improvements, 5);
}

