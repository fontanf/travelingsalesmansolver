#include "travelingsalesmansolver/lin_kernighan/sequential_move_patterns.hpp"
#include "travelingsalesmansolver/lin_kernighan/k_opt_move.hpp"

#include <gtest/gtest.h>

#include <numeric>

using namespace travelingsalesmansolver;

TEST(SequentialMovePatterns, ClosableMatchesAnalysis)
{
    // For each pattern, build a move realizing it, with no two ends at the
    // same vertex: consecutive ends are adjacent vertices if they are the
    // ends of a removed edge, else two vertices apart.
    const SequentialMovePatterns& patterns = SequentialMovePatterns::get();
    for (int pattern_id = 0; pattern_id < patterns.number_of_patterns(); ++pattern_id) {
        const SequentialMovePatterns::Pattern& pattern = patterns.pattern(pattern_id);
        int number_of_ends = pattern.ends.size();
        int m = pattern.number_of_removed_edges;
        std::vector<VertexId> t(number_of_ends);
        VertexId vertex_id = 0;
        for (int pos = 0; pos < number_of_ends; ++pos) {
            int end = pattern.ends[pos];
            if (pos > 0)
                vertex_id += ((pattern.ends[pos - 1] ^ 1) == end)? 1: 3;
            t[end] = vertex_id;
        }
        std::vector<VertexId> vertices(vertex_id + 3);
        std::iota(vertices.begin(), vertices.end(), 0);
        for (bool reversed: {false, true}) {
            // In a reversed tour, the move goes backward.
            std::vector<VertexId> order = vertices;
            if (reversed)
                std::reverse(order.begin(), order.end());
            TwoLevelList tour(order);
            SequentialMove move;
            EXPECT_EQ(move.analyze(tour, t.data(), m), pattern.closable)
                << "pattern " << pattern_id << " reversed " << reversed;
        }
    }
}

TEST(SequentialMovePatterns, NumberOfClosablePatterns)
{
    // Number of cases of valid sequential 2-, 3-, 4- and 5-opt moves (148
    // for 5-opt moves, as reported by Helsgaun).
    const SequentialMovePatterns& patterns = SequentialMovePatterns::get();
    std::vector<int> counts(SequentialMovePatterns::maximum_number_of_removed_edges + 1, 0);
    for (int pattern_id = 0; pattern_id < patterns.number_of_patterns(); ++pattern_id) {
        const SequentialMovePatterns::Pattern& pattern = patterns.pattern(pattern_id);
        if (pattern.closable)
            counts[pattern.number_of_removed_edges]++;
    }
    EXPECT_EQ(counts[2], 1);
    EXPECT_EQ(counts[3], 4);
    EXPECT_EQ(counts[4], 20);
    EXPECT_EQ(counts[5], 148);
}
