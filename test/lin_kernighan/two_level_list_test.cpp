#include "travelingsalesmansolver/lin_kernighan/two_level_list.hpp"

#include <gtest/gtest.h>

#include <numeric>
#include <random>

using namespace travelingsalesmansolver;

namespace
{

/** Reference implementation: the tour as a plain array. */
class ArrayTour
{

public:

    explicit ArrayTour(const std::vector<VertexId>& tour):
        tour_(tour),
        positions_(tour.size())
    {
        update_positions();
    }

    VertexId next(VertexId vertex_id) const
    {
        return tour_[(positions_[vertex_id] + 1) % tour_.size()];
    }

    VertexId previous(VertexId vertex_id) const
    {
        return tour_[(positions_[vertex_id] + tour_.size() - 1) % tour_.size()];
    }

    /** Reverse the path from 'vertex_id_1' to 'vertex_id_2' (following 'next'). */
    void reverse(VertexId vertex_id_1, VertexId vertex_id_2)
    {
        VertexId n = tour_.size();
        VertexId pos_1 = positions_[vertex_id_1];
        VertexId length = (positions_[vertex_id_2] - pos_1 + n) % n + 1;
        for (VertexId i = 0; i < length / 2; ++i)
            std::swap(tour_[(pos_1 + i) % n], tour_[(pos_1 + length - 1 - i) % n]);
        update_positions();
    }

private:

    void update_positions()
    {
        for (VertexId pos = 0; pos < (VertexId)tour_.size(); ++pos)
            positions_[tour_[pos]] = pos;
    }

    std::vector<VertexId> tour_;

    std::vector<VertexId> positions_;

};

std::vector<VertexId> random_tour(VertexId number_of_vertices, std::mt19937_64& generator)
{
    std::vector<VertexId> tour(number_of_vertices);
    std::iota(tour.begin(), tour.end(), 0);
    std::shuffle(tour.begin(), tour.end(), generator);
    return tour;
}

/**
 * Check that both tours are the same cycle (possibly in opposite
 * directions), and that 'next', 'previous', 'between' and 'tour' are
 * consistent with each other.
 */
void check(const TwoLevelList& list, const ArrayTour& reference, std::mt19937_64& generator)
{
    VertexId n = list.number_of_vertices();
    if (n <= 2)
        return;
    // Same cycle, in one direction for all vertices.
    bool same_direction = (list.next(0) == reference.next(0));
    for (VertexId vertex_id = 0; vertex_id < n; ++vertex_id) {
        ASSERT_EQ(list.previous(list.next(vertex_id)), vertex_id);
        if (same_direction) {
            ASSERT_EQ(list.next(vertex_id), reference.next(vertex_id));
        } else {
            ASSERT_EQ(list.next(vertex_id), reference.previous(vertex_id));
        }
    }
    // 'between' against positions in the list's own tour.
    std::vector<VertexId> tour = list.tour(0);
    ASSERT_EQ((VertexId)tour.size(), n);
    std::vector<VertexId> positions(n);
    for (VertexId pos = 0; pos < n; ++pos)
        positions[tour[pos]] = pos;
    std::uniform_int_distribution<VertexId> distribution(0, n - 1);
    for (int k = 0; k < 20; ++k) {
        VertexId a = distribution(generator);
        VertexId b = distribution(generator);
        VertexId c = distribution(generator);
        VertexId offset_b = (positions[b] - positions[a] + n) % n;
        VertexId offset_c = (positions[c] - positions[a] + n) % n;
        ASSERT_EQ(list.between(a, b, c), offset_b <= offset_c)
            << "a " << a << " b " << b << " c " << c;
    }
}

}

TEST(TwoLevelList, InitialTour)
{
    std::mt19937_64 generator(0);
    for (VertexId n: {1, 2, 3, 10, 100}) {
        std::vector<VertexId> tour = random_tour(n, generator);
        TwoLevelList list(tour);
        EXPECT_EQ(list.tour(tour[0]), tour);
        for (VertexId pos = 0; pos < n; ++pos) {
            EXPECT_EQ(list.next(tour[pos]), tour[(pos + 1) % n]);
            EXPECT_EQ(list.previous(tour[pos]), tour[(pos + n - 1) % n]);
        }
    }
}

TEST(TwoLevelList, InvalidTour)
{
    EXPECT_THROW(TwoLevelList({0, 1, 1}), std::invalid_argument);
    EXPECT_THROW(TwoLevelList({0, 1, 3}), std::invalid_argument);
}

TEST(TwoLevelList, RandomReversals)
{
    std::mt19937_64 generator(1);
    for (VertexId n: {3, 4, 5, 10, 50, 200, 1000}) {
        for (VertexId segment_size: {(VertexId)-1, (VertexId)1, (VertexId)2, (VertexId)3, n}) {
            std::vector<VertexId> tour = random_tour(n, generator);
            TwoLevelList list(tour, segment_size);
            ArrayTour reference(tour);
            std::uniform_int_distribution<VertexId> distribution(0, n - 1);
            for (int k = 0; k < 500; ++k) {
                VertexId vertex_id_1 = distribution(generator);
                VertexId vertex_id_2 = distribution(generator);
                // The reference reverses the same path, in the list's own
                // direction.
                if (list.next(0) == reference.next(0)) {
                    reference.reverse(vertex_id_1, vertex_id_2);
                } else {
                    reference.reverse(vertex_id_2, vertex_id_1);
                }
                list.reverse(vertex_id_1, vertex_id_2);
                check(list, reference, generator);
                if (::testing::Test::HasFatalFailure()) {
                    FAIL() << "n " << n << " segment size " << segment_size << " step " << k;
                }
                if (k % 100 == 99)
                    list.rebuild();
            }
        }
    }
}

TEST(TwoLevelList, RandomTwoOptMoves)
{
    std::mt19937_64 generator(2);
    for (VertexId n: {4, 5, 10, 100, 1000}) {
        std::vector<VertexId> tour = random_tour(n, generator);
        TwoLevelList list(tour);
        std::uniform_int_distribution<VertexId> distribution(0, n - 1);
        for (int k = 0; k < 1000; ++k) {
            VertexId t1 = distribution(generator);
            VertexId t3 = distribution(generator);
            bool forward = (distribution(generator) % 2 == 0);
            VertexId t2 = (forward)? list.next(t1): list.previous(t1);
            VertexId t4 = (forward)? list.previous(t3): list.next(t3);
            if (t3 == t1 || t3 == t2 || t4 == t2)
                continue;
            list.two_opt_move(t1, t2, t3, t4);
            // The edges (t2, t3) and (t4, t1) are in the tour, the removed
            // ones aren't.
            auto adjacent = [&list](VertexId u, VertexId v)
            {
                return list.next(u) == v || list.previous(u) == v;
            };
            ASSERT_TRUE(adjacent(t2, t3));
            ASSERT_TRUE(adjacent(t4, t1));
            ASSERT_FALSE(adjacent(t1, t2));
            ASSERT_FALSE(adjacent(t3, t4));
        }
        // Still a single cycle through every vertex.
        std::vector<VertexId> final_tour = list.tour(0);
        std::vector<bool> seen(n, false);
        for (VertexId vertex_id: final_tour)
            seen[vertex_id] = true;
        EXPECT_EQ(std::count(seen.begin(), seen.end(), true), n);
    }
}
