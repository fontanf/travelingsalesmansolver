#pragma once

#include "travelingsalesmansolver/lin_kernighan/two_level_list.hpp"

#include <array>
#include <cstdint>
#include <unordered_map>
#include <vector>

namespace travelingsalesmansolver
{

/**
 * Sequential k-opt moves (Helsgaun, "General k-opt submoves for the
 * Lin-Kernighan TSP heuristic", Mathematical Programming Computation, 2009).
 *
 * A sequential k-opt move is given by 2k vertices t[0], ..., t[2k - 1]: it
 * removes the tour edges (t[2i], t[2i + 1]) and adds the edges
 * (t[2i + 1], t[2i + 2]), the last one being (t[2k - 1], t[0]). A vertex
 * may appear twice (when both its tour edges are removed).
 *
 * The removed edges cut the tour into k segments, and the added edges
 * reconnect their ends. The move is valid iff the result is a single cycle.
 */
class SequentialMove
{

public:

    /** Maximum value of k. */
    static constexpr int maximum_k = 6;

    /**
     * Analyze the move 't' (2k vertices) on the current tour. Return 'true'
     * iff it's valid (the removed edges must be tour edges).
     */
    bool analyze(
            const TwoLevelList& tour,
            const VertexId* t,
            int k)
    {
        start(tour, t);
        return is_valid(tour, t, k);
    }

    /**
     * Start building a move incrementally from its first removed edge
     * (t[0], t[1]).
     */
    void start(
            const TwoLevelList& tour,
            const VertexId* t)
    {
        origin_ = t[0];
        origin_position_ = tour.position(origin_);
        tour_size_ = tour.number_of_segments() * (tour.number_of_vertices() + 1);
        set_end(tour, t, 0);
        set_end(tour, t, 1);
        std::array<int, 2 * maximum_k>& order = orders_[0];
        order[0] = 0;
        order[1] = 1;
        if (before(tour, 1, 0))
            std::swap(order[0], order[1]);
        sorted_level_ = 0;
    }

    /**
     * Signal that the removed edges from 'level' on have changed (the sorted
     * ends are then recomputed from that level when needed).
     */
    inline void invalidate(int level)
    {
        if (sorted_level_ > level - 1)
            sorted_level_ = level - 1;
    }

private:

    /**
     * Add the removed edge (t[2 * level], t[2 * level + 1]) to the move
     * built incrementally (whose first 'level' removed edges are set).
     */
    void extend(
            const TwoLevelList& tour,
            const VertexId* t,
            int level)
    {
        const std::array<int, 2 * maximum_k>& previous_order = orders_[level - 1];
        std::array<int, 2 * maximum_k>& order = orders_[level];
        int number_of_ends = 2 * level;
        for (int pos = 0; pos < number_of_ends; ++pos)
            order[pos] = previous_order[pos];
        for (int end = 2 * level; end <= 2 * level + 1; ++end) {
            set_end(tour, t, end);
            int j = number_of_ends - 1;
            while (j >= 0 && before(tour, end, order[j])) {
                order[j + 1] = order[j];
                --j;
            }
            order[j + 1] = end;
            number_of_ends++;
        }
    }

public:

    /**
     * Is the move built incrementally (from 'start'; call 'invalidate'
     * when changing its removed edges), closed after its first k removed
     * edges, valid?
     */
    bool is_valid(
            const TwoLevelList& tour,
            const VertexId* t,
            int k)
    {
        if (!build_segments(tour, t, k))
            return false;

        // Follow the new cycle from the first end of segment 0: along the
        // segment, then through the added edge leaving its other end.
        int number_of_visited_segments = 0;
        int end = segments_[0].first;
        for (;;) {
            int other_end = (ends_[end].is_first)?
                segments_[ends_[end].segment_id].second:
                segments_[ends_[end].segment_id].first;
            number_of_visited_segments++;
            int next_end = added_partner(other_end);
            if (next_end == segments_[0].first)
                break;
            if (number_of_visited_segments >= k)
                return false;
            end = next_end;
        }
        return number_of_visited_segments == k;
    }

    /**
     * Count the cycles made by the move built incrementally, closed after its
     * first k removed edges (0 if its removed edges aren't tour edges).
     */
    int count_cycles(
            const TwoLevelList& tour,
            const VertexId* t,
            int k)
    {
        if (!build_segments(tour, t, k))
            return 0;
        for (int segment_id = 0; segment_id < k; ++segment_id)
            segment_cycles_[segment_id] = -1;
        int number_of_cycles = 0;
        for (int segment_id = 0; segment_id < k; ++segment_id) {
            if (segment_cycles_[segment_id] != -1)
                continue;
            int start_end = segments_[segment_id].first;
            int end = start_end;
            do {
                segment_cycles_[ends_[end].segment_id] = number_of_cycles;
                int other_end = (ends_[end].is_first)?
                    segments_[ends_[end].segment_id].second:
                    segments_[ends_[end].segment_id].first;
                end = added_partner(other_end);
            } while (end != start_end);
            number_of_cycles++;
        }
        return number_of_cycles;
    }

    /** Number of segments of the last analyzed move. */
    inline int number_of_segments() const { return k_; }

    /** First vertex (along 'next') of a segment of the last analyzed move. */
    inline VertexId segment_first_vertex(int segment_id) const
    {
        return ends_[segments_[segment_id].first].vertex_id;
    }

    /** Last vertex (along 'next') of a segment of the last analyzed move. */
    inline VertexId segment_last_vertex(int segment_id) const
    {
        return ends_[segments_[segment_id].second].vertex_id;
    }

    /** Cycle of a segment, after 'count_cycles'. */
    inline int segment_cycle(int segment_id) const
    {
        return segment_cycles_[segment_id];
    }

    /** Segment containing a vertex, for the last analyzed move. */
    int segment(
            const TwoLevelList& tour,
            VertexId vertex_id) const
    {
        int number_of_ends = 2 * k_;
        for (int end = 0; end < number_of_ends; ++end)
            if (ends_[end].vertex_id == vertex_id)
                return ends_[end].segment_id;
        // Number of ends before the vertex, along the tour from the origin.
        int64_t position = tour.position(vertex_id);
        if (position < origin_position_)
            position += tour_size_;
        int number_of_ends_before = 0;
        for (int end = 0; end < number_of_ends; ++end)
            if (!ends_[end].wraps && ends_[end].position < position)
                number_of_ends_before++;
        // The vertex is between the ends of ranks r and r + 1.
        int rank = (number_of_ends_before == 0)? number_of_ends - 1: number_of_ends_before - 1;
        int shifted_rank = rank - shift_;
        if (shifted_rank < 0)
            shifted_rank += number_of_ends;
        return shifted_rank / 2;
    }

    /**
     * Use explicit added edges instead of the sequential ones: end 'end' is
     * joined to end 'partners[end]', for the 2k ends of the next analyzed
     * moves, until 'clear_added_partners' is called.
     */
    void set_added_partners(
            const int* partners,
            int k)
    {
        for (int end = 0; end < 2 * k; ++end)
            partners_[end] = partners[end];
        use_partners_ = true;
    }

    /** Back to sequential moves. */
    inline void clear_added_partners() { use_partners_ = false; }

    /**
     * Apply the last analyzed move (which must be valid) to 'tour', as a
     * sequence of path reversals. 'reverse(u, v)' must reverse the path from
     * 'u' to 'v' following 'next' (it may reverse its complement instead).
     */
    template <typename Reverse>
    void apply(
            const TwoLevelList& tour,
            Reverse reverse)
    {
        int k = k_;
        if (k <= 1)
            return;
        // Target: the new cycle, read from segment 0 traversed forward, as
        // signed segments (+ forward, - reversed); segment 0 stays fixed.
        std::array<int, maximum_k> target;
        int end = segments_[0].first;
        for (int pos = 0; pos < k - 1; ++pos) {
            int other_end = (ends_[end].is_first)?
                segments_[ends_[end].segment_id].second:
                segments_[ends_[end].segment_id].first;
            int next_end = added_partner(other_end);
            int segment_id = ends_[next_end].segment_id;
            target[pos] = (ends_[next_end].is_first)? segment_id: -segment_id;
            end = next_end;
        }
        std::vector<std::array<int, 2>> reversals = find_reversals(target, k - 1);

        // Current arrangement of segments 1, ..., k - 1.
        std::array<int, maximum_k> arrangement;
        for (int pos = 0; pos < k - 1; ++pos)
            arrangement[pos] = pos + 1;
        VertexId segment_0_last = ends_[segments_[0].second].vertex_id;
        for (const auto& reversal: reversals) {
            int pos_1 = reversal[0];
            int pos_2 = reversal[1];
            VertexId first_vertex_id = first_vertex(arrangement[pos_1]);
            VertexId last_vertex_id = last_vertex(arrangement[pos_2]);
            // The tour may go through the arrangement in either direction.
            bool forward = (tour.next(segment_0_last) == first_vertex(arrangement[0]));
            if (forward) {
                reverse(first_vertex_id, last_vertex_id);
            } else {
                reverse(last_vertex_id, first_vertex_id);
            }
            for (int i = pos_1, j = pos_2; i <= j; ++i, --j) {
                int tmp = arrangement[i];
                arrangement[i] = -arrangement[j];
                arrangement[j] = -tmp;
            }
        }
    }

private:

    /**
     * Sort the ends of the first k removed edges of the move built
     * incrementally, and build its segments. Return 'false' if the removed
     * edges aren't consistent with the tour.
     */
    bool build_segments(
            const TwoLevelList& tour,
            const VertexId* t,
            int k,
            bool check = true)
    {
        for (int level = sorted_level_ + 1; level < k; ++level)
            extend(tour, t, level);
        if (sorted_level_ < k - 1)
            sorted_level_ = k - 1;
        k_ = k;
        int number_of_ends = 2 * k;
        const std::array<int, 2 * maximum_k>& order = orders_[k - 1];

        // Segments: the paths between consecutive removed edges. If the
        // first end's removed edge goes forward, the removed edges are the
        // pairs of ranks (0, 1), (2, 3), ... and the segments go from rank
        // 1 to 2, 3 to 4, ..., 2k - 1 to 0. Otherwise, the removed edges
        // are the pairs (2k - 1, 0), (1, 2), ... and the segments go from
        // rank 0 to 1, 2 to 3, ...
        shift_ = (ends_[order[0]].backward)? 0: 1;
        for (int segment_id = 0; segment_id < k; ++segment_id) {
            // (No modulo: divisions are slow.)
            int first_rank = 2 * segment_id + shift_;
            int last_rank = (first_rank + 1 == number_of_ends)? 0: first_rank + 1;
            int first_end = order[first_rank];
            int last_end = order[last_rank];
            segments_[segment_id] = {first_end, last_end};
            ends_[first_end].segment_id = segment_id;
            ends_[last_end].segment_id = segment_id;
            ends_[first_end].is_first = true;
            ends_[last_end].is_first = false;
        }
        // Check that the removed edges are the pairs expected above.
        for (int segment_id = 0; check && segment_id < k; ++segment_id) {
            int last_end = segments_[segment_id].second;
            int first_end_next = segments_[(segment_id + 1 == k)? 0: segment_id + 1].first;
            if ((last_end ^ 1) != first_end_next)
                return false;
        }
        return true;
    }

    struct End
    {
        VertexId vertex_id;
        int64_t position;
        bool backward;
        bool wraps;
        int segment_id;
        bool is_first;
    };

    /** Set the end 't[pos]'. */
    inline void set_end(
            const TwoLevelList& tour,
            const VertexId* t,
            int pos)
    {
        VertexId vertex_id = t[pos];
        VertexId partner_id = t[pos ^ 1];
        ends_[pos].vertex_id = vertex_id;
        // An end whose removed edge goes back in the tour comes first
        // among the ends of the same vertex; the one of the origin whose
        // removed edge goes back comes last of all (the tour wraps).
        ends_[pos].backward = (tour.previous(vertex_id) == partner_id);
        ends_[pos].wraps = (vertex_id == origin_ && ends_[pos].backward);
        // Position along the tour, from the origin.
        int64_t position = tour.position(vertex_id);
        ends_[pos].position = (position < origin_position_)?
            position + tour_size_:
            position;
    }

    /** Is end 'end_1' before end 'end_2' along the tour, from the origin? */
    inline bool before(
            const TwoLevelList& tour,
            int end_1,
            int end_2) const
    {
        VertexId origin = origin_;
        const End& e1 = ends_[end_1];
        const End& e2 = ends_[end_2];
        (void)tour;
        (void)origin;
        if (e1.wraps != e2.wraps)
            return e2.wraps;
        if (e1.vertex_id != e2.vertex_id)
            return e1.position < e2.position;
        return e1.backward && !e2.backward;
    }

    /** Get the end at the other end of the added edge at 'end'. */
    inline int added_partner(int end) const
    {
        if (use_partners_)
            return partners_[end];
        // Added edges: (t[2i + 1], t[2i + 2]), and (t[2k - 1], t[0]).
        int number_of_ends = 2 * k_;
        if (end & 1)
            return (end + 1 == number_of_ends)? 0: end + 1;
        return (end == 0)? number_of_ends - 1: end - 1;
    }

    /** First vertex of a signed segment, in arrangement order. */
    inline VertexId first_vertex(int signed_segment) const
    {
        return (signed_segment > 0)?
            ends_[segments_[signed_segment].first].vertex_id:
            ends_[segments_[-signed_segment].second].vertex_id;
    }

    /** Last vertex of a signed segment, in arrangement order. */
    inline VertexId last_vertex(int signed_segment) const
    {
        return (signed_segment > 0)?
            ends_[segments_[signed_segment].second].vertex_id:
            ends_[segments_[-signed_segment].first].vertex_id;
    }

    /**
     * Find a shortest sequence of block reversals (position ranges, each
     * reversing the order and the directions of its segments) turning
     * (1, 2, ..., m) into 'target', by breadth-first search. Results are
     * cached.
     */
    std::vector<std::array<int, 2>> find_reversals(
            const std::array<int, maximum_k>& target,
            int m)
    {
        auto encode = [m](const std::array<int, maximum_k>& state)
        {
            uint64_t code = m;
            for (int pos = 0; pos < m; ++pos)
                code = code * 16 + (uint64_t)(state[pos] + 8);
            return code;
        };
        uint64_t target_code = encode(target);
        auto it = cache_.find(target_code);
        if (it != cache_.end())
            return it->second;

        std::array<int, maximum_k> start;
        for (int pos = 0; pos < m; ++pos)
            start[pos] = pos + 1;
        // Breadth-first search, recording each state's parent and reversal.
        std::unordered_map<uint64_t, std::pair<uint64_t, std::array<int, 2>>> parents;
        std::vector<std::array<int, maximum_k>> queue = {start};
        uint64_t start_code = encode(start);
        parents[start_code] = {start_code, {-1, -1}};
        for (size_t queue_pos = 0; queue_pos < queue.size(); ++queue_pos) {
            std::array<int, maximum_k> state = queue[queue_pos];
            uint64_t code = encode(state);
            if (code == target_code)
                break;
            for (int pos_1 = 0; pos_1 < m; ++pos_1) {
                for (int pos_2 = pos_1; pos_2 < m; ++pos_2) {
                    std::array<int, maximum_k> next_state = state;
                    for (int i = pos_1, j = pos_2; i <= j; ++i, --j) {
                        next_state[i] = -state[j];
                        next_state[j] = -state[i];
                    }
                    uint64_t next_code = encode(next_state);
                    if (parents.count(next_code))
                        continue;
                    parents[next_code] = {code, {pos_1, pos_2}};
                    queue.push_back(next_state);
                }
            }
        }
        std::vector<std::array<int, 2>> reversals;
        for (uint64_t code = target_code; code != start_code; code = parents[code].first)
            reversals.push_back(parents[code].second);
        std::reverse(reversals.begin(), reversals.end());
        cache_[target_code] = reversals;
        return reversals;
    }

    /** k. */
    int k_ = 0;

    /** First vertex of the move, from which ends are sorted. */
    VertexId origin_ = -1;

    /** Position of 'origin_' along the tour. */
    int64_t origin_position_ = 0;

    /** Upper bound of the positions along the tour. */
    int64_t tour_size_ = 0;

    /** Levels whose sorted ends ('orders_') are up to date. */
    int sorted_level_ = -1;

    /** For each level, the ends of its first removed edges, sorted along the tour. */
    std::array<std::array<int, 2 * maximum_k>, maximum_k> orders_;

    /** Ends of the removed edges (entry 'pos' corresponds to t[pos]). */
    std::array<End, 2 * maximum_k> ends_;


    /** Segments: (first end, last end), in tour order. */
    std::array<std::pair<int, int>, maximum_k> segments_;


    /** 0 if the segments start at even ranks, 1 otherwise. */
    int shift_ = 0;

    /** Cycle of each segment, after 'count_cycles'. */
    std::array<int, maximum_k> segment_cycles_;

    /** Explicit added edges ('set_added_partners'). */
    std::array<int, 2 * maximum_k> partners_;

    /** Whether to use 'partners_' instead of the sequential added edges. */
    bool use_partners_ = false;

    /** Reversal sequences already computed, by target arrangement. */
    std::unordered_map<uint64_t, std::vector<std::array<int, 2>>> cache_;

};

}
