#include "travelingsalesmansolver/distances/distances_builder.hpp"

#include <gtest/gtest.h>

using namespace travelingsalesmansolver;

namespace
{

Distances euclidean_distances()
{
    DistancesEuc2DBuilder euc_2d_builder;
    euc_2d_builder.set_number_of_vertices(4);
    euc_2d_builder.set_coordinates(0, 0, 0);
    euc_2d_builder.set_coordinates(1, 3, 0);
    euc_2d_builder.set_coordinates(2, 3, 4);
    euc_2d_builder.set_coordinates(3, 0, 4);
    DistancesBuilder distances_builder;
    distances_builder.set_number_of_vertices(4);
    distances_builder.set_distances_euc_2d(euc_2d_builder.build());
    return distances_builder.build();
}

}

TEST(Distances, ComputeDistancesExplicitTwice)
{
    // The second call must keep the matrix computed by the first one, which
    // might be referenced.
    Distances distances = euclidean_distances();
    distances.compute_distances_explicit();
    const DistancesExplicit* distances_explicit = distances.distances_explicit();
    ASSERT_NE(distances_explicit, nullptr);
    EXPECT_EQ(distances_explicit->distance(0, 2), 5);

    distances.compute_distances_explicit();
    EXPECT_EQ(distances.distances_explicit(), distances_explicit);
    distances.compute_distances_explicit_triangle();
    EXPECT_EQ(distances.distances_explicit(), distances_explicit);
    EXPECT_EQ(distances.distances_explicit_triangle(), nullptr);
}

TEST(Distances, ComputeDistancesExplicitTriangleTwice)
{
    Distances distances = euclidean_distances();
    distances.compute_distances_explicit_triangle();
    const DistancesExplicitTriangle* distances_explicit_triangle = distances.distances_explicit_triangle();
    ASSERT_NE(distances_explicit_triangle, nullptr);
    EXPECT_EQ(distances_explicit_triangle->distance(0, 2), 5);

    distances.compute_distances_explicit_triangle();
    EXPECT_EQ(distances.distances_explicit_triangle(), distances_explicit_triangle);
}
