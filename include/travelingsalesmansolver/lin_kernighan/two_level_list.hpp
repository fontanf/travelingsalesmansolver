#pragma once

#include "travelingsalesmansolver/distances/commons.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace travelingsalesmansolver
{

/**
 * Tour represented as a two-level doubly-linked list (Fredman, Johnson &
 * McGeoch, "Data structures for traveling salesmen", 1995).
 *
 * The tour is cut into about sqrt(n) segments of consecutive vertices. Each
 * segment is a doubly-linked list of its vertices, in a "stored" order, plus
 * a reversal bit telling whether the tour goes through it in stored order or
 * in reverse. The segments are themselves linked in a cyclic list, in tour
 * order. Each vertex has a rank (its position in its segment's stored order,
 * consecutive within a segment) and each segment a rank (its position in the
 * list of segments).
 *
 * Complexity:
 * - 'next', 'previous' and 'between' are O(1);
 * - 'reverse' (and so a 2-opt move) is O(sqrt(n)): a path inside a segment
 *   is reversed vertex by vertex; otherwise the segments at both ends of the
 *   path are split so that the path consists of whole segments, and the
 *   shorter of the path and its complement is reversed segment by segment
 *   (their bits are flipped and their order in the list reversed).
 *
 * A split moves the smaller part of a segment into its neighbor, keeping the
 * number of segments constant; segment sizes drift over time, and 'rebuild'
 * resets them.
 *
 * Reversing a path or its complement gives the same cycle, but in opposite
 * directions: after a 'reverse' or a 'two_opt_move', the direction in which
 * the tour is traversed ('next' versus 'previous') is not preserved. Callers
 * must not assume it is.
 */
class TwoLevelList
{

public:

    /**
     * Build the list for a tour, given as the sequence of its vertices
     * (a permutation of 0, ..., n - 1).
     *
     * 'segment_size' is the number of vertices per segment; by default,
     * about sqrt(n).
     */
    explicit TwoLevelList(
            const std::vector<VertexId>& tour,
            VertexId segment_size = -1)
    {
        build(tour, segment_size);
    }

    /** Get the number of vertices. */
    inline VertexId number_of_vertices() const { return vertices_.size(); }

    /** Get the number of segments. */
    inline VertexId number_of_segments() const { return segments_.size(); }

    /** Get the vertex following 'vertex_id' in the tour. */
    inline VertexId next(VertexId vertex_id) const
    {
        const Vertex& vertex = vertices_[vertex_id];
        const Segment& segment = segments_[vertex.segment_id];
        VertexId next_vertex_id = (segment.reversed)?
            vertex.stored_previous:
            vertex.stored_next;
        if (next_vertex_id != -1)
            return next_vertex_id;
        return segment_first(segment.next);
    }

    /** Get the vertex preceding 'vertex_id' in the tour. */
    inline VertexId previous(VertexId vertex_id) const
    {
        const Vertex& vertex = vertices_[vertex_id];
        const Segment& segment = segments_[vertex.segment_id];
        VertexId previous_vertex_id = (segment.reversed)?
            vertex.stored_next:
            vertex.stored_previous;
        if (previous_vertex_id != -1)
            return previous_vertex_id;
        return segment_last(segment.previous);
    }

    /**
     * Return 'true' iff 'vertex_id_2' is on the path from 'vertex_id_1' to
     * 'vertex_id_3' following 'next' (ends included).
     */
    inline bool between(
            VertexId vertex_id_1,
            VertexId vertex_id_2,
            VertexId vertex_id_3) const
    {
        Key key_1 = key(vertex_id_1);
        Key key_2 = key(vertex_id_2);
        Key key_3 = key(vertex_id_3);
        if (!(key_3 < key_1))
            return !(key_2 < key_1) && !(key_3 < key_2);
        return !(key_2 < key_1) || !(key_3 < key_2);
    }

    /**
     * Reverse the path from 'vertex_id_1' to 'vertex_id_2' (following
     * 'next'): with 'a' the vertex before the path and 'd' the one after,
     * the edges (a, vertex_id_1) and (vertex_id_2, d) are replaced by
     * (a, vertex_id_2) and (vertex_id_1, d).
     *
     * Reversing the whole tour does nothing (it gives the same cycle).
     */
    void reverse(
            VertexId vertex_id_1,
            VertexId vertex_id_2)
    {
        if (vertex_id_1 == vertex_id_2)
            return;
        VertexId before_id = previous(vertex_id_1);
        VertexId after_id = next(vertex_id_2);
        if (after_id == vertex_id_1)
            return;

        // The path, or its complement, lies inside a single segment.
        if (within_segment(vertex_id_1, vertex_id_2)) {
            reverse_within_segment(vertex_id_1, vertex_id_2);
            return;
        }
        if (within_segment(after_id, before_id)) {
            reverse_within_segment(after_id, before_id);
            return;
        }

        // Split the segments so that the path consists of whole segments.
        split_before(vertex_id_1, -1);
        if (within_segment(vertex_id_1, vertex_id_2)) {
            reverse_within_segment(vertex_id_1, vertex_id_2);
            return;
        }
        // 'vertex_id_1' must stay the first vertex of its segment.
        split_before(after_id, vertices_[vertex_id_1].segment_id);
        if (within_segment(vertex_id_1, vertex_id_2)) {
            reverse_within_segment(vertex_id_1, vertex_id_2);
            return;
        }

        // Reverse the shorter of the path and its complement.
        SegmentId segment_id_1 = vertices_[vertex_id_1].segment_id;
        SegmentId segment_id_2 = vertices_[vertex_id_2].segment_id;
        SegmentId number_of_segments_path = run_size(segment_id_1, segment_id_2);
        if (2 * number_of_segments_path <= number_of_segments()) {
            reverse_segments(segment_id_1, segment_id_2);
        } else {
            reverse_segments(
                    vertices_[after_id].segment_id,
                    vertices_[before_id].segment_id);
        }
    }

    /**
     * Apply a 2-opt move, with Lin-Kernighan's naming: replace the edges
     * (t1, t2) and (t3, t4) by (t2, t3) and (t4, t1).
     *
     * For the result to be a tour, 't2' and 't4' must be on opposite sides:
     * either t2 = next(t1) and t4 = previous(t3), or t2 = previous(t1) and
     * t4 = next(t3).
     */
    void two_opt_move(
            VertexId t1,
            VertexId t2,
            VertexId t3,
            VertexId t4)
    {
        if (next(t1) == t2) {
            // t1 t2 ... t4 t3 -> t1 t4 ... t2 t3
            reverse(t2, t4);
        } else {
            // t2 t1 ... t3 t4 -> t2 t3 ... t1 t4
            reverse(t1, t3);
        }
    }

    /**
     * Get the position of a vertex along the tour, as an integer: positions
     * increase along 'next', from the first vertex of the first segment of
     * the list (a vertex whose position is not larger than the previous
     * one's starts the tour over). Valid until the tour changes.
     */
    inline int64_t position(VertexId vertex_id) const
    {
        return segments_[vertices_[vertex_id].segment_id].rank * (number_of_vertices() + 1)
            + offset(vertex_id);
    }

    /** Get the tour, as a sequence of vertices starting at 'vertex_id'. */
    std::vector<VertexId> tour(VertexId vertex_id = 0) const
    {
        std::vector<VertexId> vertex_ids;
        vertex_ids.reserve(number_of_vertices());
        VertexId current_vertex_id = vertex_id;
        for (VertexId pos = 0; pos < number_of_vertices(); ++pos) {
            vertex_ids.push_back(current_vertex_id);
            current_vertex_id = next(current_vertex_id);
        }
        return vertex_ids;
    }

    /** Rebuild the segments from the current tour (resetting their sizes). */
    void rebuild(VertexId segment_size = -1)
    {
        if (number_of_vertices() == 0)
            return;
        build(tour(), segment_size);
    }

private:

    using SegmentId = int64_t;

    /** Position of a vertex along the tour: (segment rank, offset). */
    struct Key
    {
        int64_t segment_rank;
        int64_t offset;

        inline bool operator<(const Key& key) const
        {
            return segment_rank < key.segment_rank
                || (segment_rank == key.segment_rank && offset < key.offset);
        }
    };

    struct Vertex
    {
        /** Previous vertex in its segment's stored order ('-1' if first). */
        VertexId stored_previous = -1;

        /** Next vertex in its segment's stored order ('-1' if last). */
        VertexId stored_next = -1;

        /** Segment of the vertex. */
        SegmentId segment_id = -1;

        /** Position in the segment's stored order (consecutive in a segment). */
        int64_t rank = 0;
    };

    struct Segment
    {
        /** First vertex, in stored order. */
        VertexId stored_first = -1;

        /** Last vertex, in stored order. */
        VertexId stored_last = -1;

        /** Previous segment in the list. */
        SegmentId previous = -1;

        /** Next segment in the list. */
        SegmentId next = -1;

        /** Position in the list of segments. */
        int64_t rank = 0;

        /** Number of vertices. */
        VertexId size = 0;

        /** Is the segment traversed in reverse stored order? */
        bool reversed = false;
    };

    /*
     * Private methods
     */

    /** Build the list from a tour. */
    void build(
            const std::vector<VertexId>& tour,
            VertexId segment_size)
    {
        VertexId number_of_vertices = tour.size();
        if (segment_size <= 0) {
            segment_size = (VertexId)std::sqrt((double)number_of_vertices);
            if (segment_size < 1)
                segment_size = 1;
        }
        vertices_ = std::vector<Vertex>(number_of_vertices);
        std::vector<bool> seen(number_of_vertices, false);
        for (VertexId vertex_id: tour) {
            if (vertex_id < 0
                    || vertex_id >= number_of_vertices
                    || seen[vertex_id]) {
                throw std::invalid_argument(
                        "travelingsalesmansolver::TwoLevelList: "
                        "the tour must be a permutation of 0, ..., n - 1.");
            }
            seen[vertex_id] = true;
        }
        SegmentId number_of_segments = (number_of_vertices + segment_size - 1) / segment_size;
        segments_ = std::vector<Segment>(number_of_segments);
        for (SegmentId segment_id = 0; segment_id < number_of_segments; ++segment_id) {
            Segment& segment = segments_[segment_id];
            segment.previous = (segment_id + number_of_segments - 1) % number_of_segments;
            segment.next = (segment_id + 1) % number_of_segments;
            segment.rank = segment_id;
            VertexId pos_first = segment_id * segment_size;
            VertexId pos_last = (std::min)(pos_first + segment_size, number_of_vertices) - 1;
            segment.stored_first = tour[pos_first];
            segment.stored_last = tour[pos_last];
            segment.size = pos_last - pos_first + 1;
            for (VertexId pos = pos_first; pos <= pos_last; ++pos) {
                Vertex& vertex = vertices_[tour[pos]];
                vertex.segment_id = segment_id;
                vertex.rank = pos - pos_first;
                vertex.stored_previous = (pos == pos_first)? -1: tour[pos - 1];
                vertex.stored_next = (pos == pos_last)? -1: tour[pos + 1];
            }
        }
    }

    /** Get the first vertex of a segment, in tour order. */
    inline VertexId segment_first(SegmentId segment_id) const
    {
        const Segment& segment = segments_[segment_id];
        return (segment.reversed)? segment.stored_last: segment.stored_first;
    }

    /** Get the last vertex of a segment, in tour order. */
    inline VertexId segment_last(SegmentId segment_id) const
    {
        const Segment& segment = segments_[segment_id];
        return (segment.reversed)? segment.stored_first: segment.stored_last;
    }

    /** Get the number of vertices before 'vertex_id' in its segment, in tour order. */
    inline int64_t offset(VertexId vertex_id) const
    {
        const Vertex& vertex = vertices_[vertex_id];
        const Segment& segment = segments_[vertex.segment_id];
        return (segment.reversed)?
            vertices_[segment.stored_last].rank - vertex.rank:
            vertex.rank - vertices_[segment.stored_first].rank;
    }

    /** Get the position of a vertex along the tour. */
    inline Key key(VertexId vertex_id) const
    {
        return {segments_[vertices_[vertex_id].segment_id].rank, offset(vertex_id)};
    }

    /**
     * Return 'true' iff the path from 'vertex_id_1' to 'vertex_id_2' lies
     * inside a single segment.
     */
    inline bool within_segment(
            VertexId vertex_id_1,
            VertexId vertex_id_2) const
    {
        return vertices_[vertex_id_1].segment_id == vertices_[vertex_id_2].segment_id
            && offset(vertex_id_1) <= offset(vertex_id_2);
    }

    /** Get the number of segments from 'segment_id_1' to 'segment_id_2' (following 'next'). */
    inline SegmentId run_size(
            SegmentId segment_id_1,
            SegmentId segment_id_2) const
    {
        SegmentId number_of_segments = segments_.size();
        return (segments_[segment_id_2].rank - segments_[segment_id_1].rank + number_of_segments)
            % number_of_segments + 1;
    }

    /** Append a vertex at the end of a segment, in tour order. */
    void push_back(
            SegmentId segment_id,
            VertexId vertex_id)
    {
        Segment& segment = segments_[segment_id];
        Vertex& vertex = vertices_[vertex_id];
        vertex.segment_id = segment_id;
        if (!segment.reversed) {
            Vertex& last = vertices_[segment.stored_last];
            vertex.stored_previous = segment.stored_last;
            vertex.stored_next = -1;
            vertex.rank = last.rank + 1;
            last.stored_next = vertex_id;
            segment.stored_last = vertex_id;
        } else {
            Vertex& first = vertices_[segment.stored_first];
            vertex.stored_next = segment.stored_first;
            vertex.stored_previous = -1;
            vertex.rank = first.rank - 1;
            first.stored_previous = vertex_id;
            segment.stored_first = vertex_id;
        }
        segment.size++;
    }

    /** Insert a vertex at the start of a segment, in tour order. */
    void push_front(
            SegmentId segment_id,
            VertexId vertex_id)
    {
        Segment& segment = segments_[segment_id];
        Vertex& vertex = vertices_[vertex_id];
        vertex.segment_id = segment_id;
        if (!segment.reversed) {
            Vertex& first = vertices_[segment.stored_first];
            vertex.stored_next = segment.stored_first;
            vertex.stored_previous = -1;
            vertex.rank = first.rank - 1;
            first.stored_previous = vertex_id;
            segment.stored_first = vertex_id;
        } else {
            Vertex& last = vertices_[segment.stored_last];
            vertex.stored_previous = segment.stored_last;
            vertex.stored_next = -1;
            vertex.rank = last.rank + 1;
            last.stored_next = vertex_id;
            segment.stored_last = vertex_id;
        }
        segment.size++;
    }

    /**
     * Split the segment of 'vertex_id' so that 'vertex_id' becomes the first
     * vertex of a segment (in tour order).
     *
     * The smaller part is moved into the neighboring segment: the part
     * before 'vertex_id' is appended to the previous segment, or the part
     * from 'vertex_id' on is prepended to the next segment. The latter is
     * never done when the next segment is 'protected_segment_id' (whose
     * first vertex must not change).
     */
    void split_before(
            VertexId vertex_id,
            SegmentId protected_segment_id)
    {
        SegmentId segment_id = vertices_[vertex_id].segment_id;
        if (segment_first(segment_id) == vertex_id)
            return;
        Segment& segment = segments_[segment_id];
        VertexId number_of_vertices_before = offset(vertex_id);
        VertexId number_of_vertices_after = segment.size - number_of_vertices_before;
        buffer_.clear();
        if (number_of_vertices_before <= number_of_vertices_after
                || segment.next == protected_segment_id) {
            // Move the part before 'vertex_id' to the previous segment.
            for (VertexId current_vertex_id = segment_first(segment_id);
                    current_vertex_id != vertex_id;
                    current_vertex_id = next(current_vertex_id)) {
                buffer_.push_back(current_vertex_id);
            }
            if (!segment.reversed) {
                segment.stored_first = vertex_id;
                vertices_[vertex_id].stored_previous = -1;
            } else {
                segment.stored_last = vertex_id;
                vertices_[vertex_id].stored_next = -1;
            }
            segment.size -= buffer_.size();
            SegmentId previous_segment_id = segment.previous;
            for (VertexId moved_vertex_id: buffer_)
                push_back(previous_segment_id, moved_vertex_id);
        } else {
            // Move the part from 'vertex_id' on to the next segment.
            VertexId last_vertex_id = segment_last(segment_id);
            for (VertexId current_vertex_id = vertex_id;;
                    current_vertex_id = next(current_vertex_id)) {
                buffer_.push_back(current_vertex_id);
                if (current_vertex_id == last_vertex_id)
                    break;
            }
            VertexId previous_vertex_id = previous(vertex_id);
            if (!segment.reversed) {
                segment.stored_last = previous_vertex_id;
                vertices_[previous_vertex_id].stored_next = -1;
            } else {
                segment.stored_first = previous_vertex_id;
                vertices_[previous_vertex_id].stored_previous = -1;
            }
            segment.size -= buffer_.size();
            SegmentId next_segment_id = segment.next;
            for (auto it = buffer_.rbegin(); it != buffer_.rend(); ++it)
                push_front(next_segment_id, *it);
        }
    }

    /**
     * Reverse the path from 'vertex_id_1' to 'vertex_id_2', which lies
     * inside a single segment, vertex by vertex.
     */
    void reverse_within_segment(
            VertexId vertex_id_1,
            VertexId vertex_id_2)
    {
        SegmentId segment_id = vertices_[vertex_id_1].segment_id;
        Segment& segment = segments_[segment_id];
        // Ends of the path in stored order.
        VertexId stored_first_id = (segment.reversed)? vertex_id_2: vertex_id_1;
        VertexId stored_last_id = (segment.reversed)? vertex_id_1: vertex_id_2;
        VertexId before_id = vertices_[stored_first_id].stored_previous;
        VertexId after_id = vertices_[stored_last_id].stored_next;
        int64_t rank = vertices_[stored_first_id].rank;
        buffer_.clear();
        for (VertexId current_vertex_id = stored_first_id;;
                current_vertex_id = vertices_[current_vertex_id].stored_next) {
            buffer_.push_back(current_vertex_id);
            if (current_vertex_id == stored_last_id)
                break;
        }
        std::reverse(buffer_.begin(), buffer_.end());
        VertexId size = buffer_.size();
        for (VertexId pos = 0; pos < size; ++pos) {
            Vertex& vertex = vertices_[buffer_[pos]];
            vertex.stored_previous = (pos == 0)? before_id: buffer_[pos - 1];
            vertex.stored_next = (pos == size - 1)? after_id: buffer_[pos + 1];
            vertex.rank = rank + pos;
        }
        if (before_id != -1) {
            vertices_[before_id].stored_next = buffer_.front();
        } else {
            segment.stored_first = buffer_.front();
        }
        if (after_id != -1) {
            vertices_[after_id].stored_previous = buffer_.back();
        } else {
            segment.stored_last = buffer_.back();
        }
    }

    /**
     * Reverse the run of whole segments from 'segment_id_1' to
     * 'segment_id_2' (following 'next'), which is not the whole list.
     */
    void reverse_segments(
            SegmentId segment_id_1,
            SegmentId segment_id_2)
    {
        SegmentId previous_segment_id = segments_[segment_id_1].previous;
        SegmentId next_segment_id = segments_[segment_id_2].next;
        segment_buffer_.clear();
        rank_buffer_.clear();
        for (SegmentId segment_id = segment_id_1;; segment_id = segments_[segment_id].next) {
            segment_buffer_.push_back(segment_id);
            rank_buffer_.push_back(segments_[segment_id].rank);
            if (segment_id == segment_id_2)
                break;
        }
        SegmentId size = segment_buffer_.size();
        for (SegmentId pos = 0; pos < size; ++pos) {
            Segment& segment = segments_[segment_buffer_[pos]];
            std::swap(segment.previous, segment.next);
            segment.reversed = !segment.reversed;
            // The run's positions in the list are reused in reverse order.
            segment.rank = rank_buffer_[size - 1 - pos];
        }
        // New order: previous, segment_id_2, ..., segment_id_1, next.
        segments_[segment_id_2].previous = previous_segment_id;
        segments_[previous_segment_id].next = segment_id_2;
        segments_[segment_id_1].next = next_segment_id;
        segments_[next_segment_id].previous = segment_id_1;
    }

    /*
     * Private attributes
     */

    /** Vertices. */
    std::vector<Vertex> vertices_;

    /** Segments. */
    std::vector<Segment> segments_;

    /** Buffer of vertices, used by the operations moving or reversing vertices. */
    std::vector<VertexId> buffer_;

    /** Buffers of segments and ranks, used by 'reverse_segments'. */
    std::vector<SegmentId> segment_buffer_;
    std::vector<int64_t> rank_buffer_;

};

}
