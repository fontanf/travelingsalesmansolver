#pragma once

#include <array>
#include <map>
#include <vector>

namespace travelingsalesmansolver
{

/**
 * Patterns of sequential moves: all the cases of the relative positions,
 * along the tour, of the ends of their first removed edges. They let a
 * search enumerate sequential moves without analyzing each of them.
 *
 * A move t[0], ..., t[2m - 1] removes the tour edges (t[2i], t[2i + 1]) and
 * adds the edges (t[2i + 1], t[2i + 2]). Its pattern is the sequence of its
 * 2m ends (end e being t[e]) in tour order, starting at end 0, in the
 * direction where t[1] follows t[0] ("forward"); two ends at the same vertex
 * are ordered as the vertex's two removed edges (the one going backward
 * first). The pattern always starts with ends 0 and 1, and the two ends of a
 * removed edge are consecutive.
 *
 * The gap j of a pattern is the part of the tour between its ends of
 * positions j and j + 1 (modulo 2m); it's a segment (a path of the tour) if
 * these two ends aren't the ends of a removed edge.
 *
 * Adding a removed edge (t_new, t_next), t_new being in segment j, gives the
 * pattern with ends 2m and 2m + 1 inserted in gap j, in the order of t_new
 * and t_next along the tour.
 *
 * With m removed edges, and the m - 1 added edges between them, the tour
 * becomes a path from t[0] to t[2m - 1], and possibly cycles. The move
 * closed by adding (t[2m - 1], t[0]) is valid if the result is a tour.
 * Since each further removed edge can join at most one cycle to the path, a
 * move with c cycles can't be closed validly with fewer than m + c removed
 * edges.
 */
class SequentialMovePatterns
{

public:

    /** Maximum number of removed edges of a pattern. */
    static constexpr int maximum_number_of_removed_edges = 5;

    struct Pattern
    {
        /** Number of removed edges. */
        int number_of_removed_edges = 0;

        /** Ends, in tour order from end 0. */
        std::vector<int> ends;

        /**
         * For each gap j and side s, the pattern obtained by removing an edge
         * (t_new, t_next) with t_new in gap j, and t_next following t_new
         * (s = 0) or preceding it (s = 1) along the forward direction; -1 if
         * gap j isn't a segment (or if the number of removed edges is
         * maximum).
         */
        std::array<std::array<int, 2>, 2 * maximum_number_of_removed_edges> next;

        /**
         * For each end, the gap of the segment containing it, and the side of
         * its tour edge in that segment (0: following it, 1: preceding it).
         */
        std::array<int, 2 * maximum_number_of_removed_edges> end_gaps;
        std::array<int, 2 * maximum_number_of_removed_edges> end_sides;

        /** Whether the move closed by (t[2m - 1], t[0]) is valid. */
        bool closable = false;

        /** Number of cycles (besides the path from t[0] to t[2m - 1]). */
        int number_of_cycles = 0;
    };

    /** Get the patterns (computed once). */
    static const SequentialMovePatterns& get()
    {
        static const SequentialMovePatterns patterns;
        return patterns;
    }

    /** Pattern of a single removed edge. */
    static constexpr int first_pattern_id = 0;

    /** Get a pattern. */
    inline const Pattern& pattern(int pattern_id) const { return patterns_[pattern_id]; }

    /** Number of patterns. */
    inline int number_of_patterns() const { return patterns_.size(); }

private:

    /** Constructor: enumerate the patterns reachable from the first one. */
    SequentialMovePatterns()
    {
        std::map<std::vector<int>, int> ids;
        auto add = [this, &ids](const std::vector<int>& ends)
        {
            auto it = ids.find(ends);
            if (it != ids.end())
                return it->second;
            int pattern_id = patterns_.size();
            ids[ends] = pattern_id;
            Pattern pattern;
            pattern.number_of_removed_edges = ends.size() / 2;
            pattern.ends = ends;
            analyze(pattern);
            patterns_.push_back(pattern);
            return pattern_id;
        };
        add({0, 1});
        // Breadth-first: 'patterns_' grows while it's traversed.
        for (int pattern_id = 0; pattern_id < (int)patterns_.size(); ++pattern_id) {
            std::vector<int> ends = patterns_[pattern_id].ends;
            int number_of_ends = ends.size();
            int m = number_of_ends / 2;
            for (int gap = 0; gap < 2 * maximum_number_of_removed_edges; ++gap)
                patterns_[pattern_id].next[gap] = {-1, -1};
            if (m == maximum_number_of_removed_edges)
                continue;
            for (int gap = 0; gap < number_of_ends; ++gap) {
                if (!is_segment(ends, gap))
                    continue;
                for (int side = 0; side < 2; ++side) {
                    std::vector<int> next_ends(ends.begin(), ends.begin() + gap + 1);
                    if (side == 0) {
                        next_ends.push_back(2 * m);
                        next_ends.push_back(2 * m + 1);
                    } else {
                        next_ends.push_back(2 * m + 1);
                        next_ends.push_back(2 * m);
                    }
                    next_ends.insert(next_ends.end(), ends.begin() + gap + 1, ends.end());
                    int next_id = add(next_ends);
                    patterns_[pattern_id].next[gap][side] = next_id;
                }
            }
        }
    }

    /** Are the ends of positions gap and gap + 1 not the ends of a removed edge? */
    static bool is_segment(
            const std::vector<int>& ends,
            int gap)
    {
        int number_of_ends = ends.size();
        // With a single removed edge, both gaps are between its ends: t[1]
        // follows t[0] along the removed edge.
        if (number_of_ends == 2)
            return gap == 1;
        int end_1 = ends[gap];
        int end_2 = ends[(gap + 1) % number_of_ends];
        return (end_1 ^ 1) != end_2;
    }

    /** Compute the fields of a pattern from its ends. */
    static void analyze(Pattern& pattern)
    {
        const std::vector<int>& ends = pattern.ends;
        int number_of_ends = ends.size();
        int m = number_of_ends / 2;

        // Neighbors of each end: along its segment, and along its added edge
        // (-1 for the ends of the path, t[0] and t[2m - 1]).
        std::vector<int> segment_neighbors(number_of_ends);
        for (int gap = 0; gap < number_of_ends; ++gap) {
            if (!is_segment(ends, gap))
                continue;
            int end_1 = ends[gap];
            int end_2 = ends[(gap + 1) % number_of_ends];
            segment_neighbors[end_1] = end_2;
            segment_neighbors[end_2] = end_1;
            pattern.end_gaps[end_1] = gap;
            pattern.end_sides[end_1] = 0;
            pattern.end_gaps[end_2] = gap;
            pattern.end_sides[end_2] = 1;
        }
        auto added_neighbor = [number_of_ends](int end)
        {
            if (end == 0 || end == number_of_ends - 1)
                return -1;
            return (end % 2 == 1)? end + 1: end - 1;
        };

        // Follow the path from end 0; the other ends are on cycles.
        std::vector<bool> visited(number_of_ends, false);
        int end = 0;
        for (;;) {
            visited[end] = true;
            int other_end = segment_neighbors[end];
            visited[other_end] = true;
            end = added_neighbor(other_end);
            if (end == -1)
                break;
        }
        pattern.number_of_cycles = 0;
        for (int start_end = 0; start_end < number_of_ends; ++start_end) {
            if (visited[start_end])
                continue;
            pattern.number_of_cycles++;
            int end = start_end;
            do {
                visited[end] = true;
                int other_end = segment_neighbors[end];
                visited[other_end] = true;
                end = added_neighbor(other_end);
            } while (end != start_end);
        }
        // Closing (t[2m - 1], t[0]) joins the ends of the path.
        pattern.closable = (pattern.number_of_cycles == 0);
        (void)m;
    }

    /** Patterns. */
    std::vector<Pattern> patterns_;

};

}
