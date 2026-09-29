#pragma once

#include "travelingsalesmansolver/lin_kernighan/candidates.hpp"

#include <algorithm>
#include <array>
#include <numeric>
#include <random>
#include <tuple>
#include <vector>

namespace travelingsalesmansolver
{

/**
 * Minimum spanning forest of the candidate graph (Kruskal): for each vertex,
 * its neighbors in the forest.
 */
template <typename Distances>
std::vector<std::vector<VertexId>> minimum_spanning_forest(
        const Distances& distances,
        const CandidateLists& candidates)
{
    VertexId number_of_vertices = candidates.size();
    std::vector<std::tuple<Distance, VertexId, VertexId>> edges;
    for (VertexId vertex_id_1 = 0; vertex_id_1 < number_of_vertices; ++vertex_id_1) {
        for (VertexId vertex_id_2: candidates[vertex_id_1]) {
            VertexId u = (std::min)(vertex_id_1, vertex_id_2);
            VertexId v = (std::max)(vertex_id_1, vertex_id_2);
            edges.emplace_back(distances.distance(u, v), u, v);
        }
    }
    std::sort(edges.begin(), edges.end());
    edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
    std::vector<VertexId> parents(number_of_vertices);
    std::iota(parents.begin(), parents.end(), 0);
    auto find = [&parents](VertexId vertex_id)
    {
        while (parents[vertex_id] != vertex_id) {
            parents[vertex_id] = parents[parents[vertex_id]];
            vertex_id = parents[vertex_id];
        }
        return vertex_id;
    };
    std::vector<std::vector<VertexId>> neighbors(number_of_vertices);
    for (const auto& edge: edges) {
        VertexId root_1 = find(std::get<1>(edge));
        VertexId root_2 = find(std::get<2>(edge));
        if (root_1 == root_2)
            continue;
        parents[root_1] = root_2;
        neighbors[std::get<1>(edge)].push_back(std::get<2>(edge));
        neighbors[std::get<2>(edge)].push_back(std::get<1>(edge));
    }
    return neighbors;
}

/**
 * Merge two tours by partition crossover (GPX; Whitley, Hains & Howe,
 * "Tunneling between optima: partition crossover for the traveling
 * salesman problem", 2009).
 *
 * Removing the edges common to both tours splits the vertices into
 * components (connected by the other edges of either tour). A component
 * that the tours enter and leave through exactly two common edges is
 * crossed by each tour as a single path between the same two vertices, so
 * either path can be used, independently of the other components. The
 * result keeps 'tour_1' everywhere, except in such components where the
 * path of 'tour_2' is shorter. It is never longer than 'tour_1', nor than
 * 'tour_2' when all components qualify.
 *
 * Tours are sequences of vertices; so is the result.
 */
template <typename Distances>
std::vector<VertexId> partition_crossover(
        const Distances& distances,
        const std::vector<VertexId>& tour_1,
        const std::vector<VertexId>& tour_2)
{
    VertexId number_of_vertices = tour_1.size();
    if (number_of_vertices < 4)
        return tour_1;
    auto neighbors_in = [number_of_vertices](const std::vector<VertexId>& tour)
    {
        std::vector<std::array<VertexId, 2>> neighbors(number_of_vertices);
        for (VertexId pos = 0; pos < number_of_vertices; ++pos) {
            VertexId vertex_id = tour[pos];
            neighbors[vertex_id][0] = tour[(pos + number_of_vertices - 1) % number_of_vertices];
            neighbors[vertex_id][1] = tour[(pos + 1) % number_of_vertices];
        }
        return neighbors;
    };
    std::vector<std::array<VertexId, 2>> neighbors_1 = neighbors_in(tour_1);
    std::vector<std::array<VertexId, 2>> neighbors_2 = neighbors_in(tour_2);
    auto in_tour = [](const std::vector<std::array<VertexId, 2>>& neighbors, VertexId u, VertexId v)
    {
        return neighbors[u][0] == v || neighbors[u][1] == v;
    };

    // Components of the graph of the edges of either tour which aren't
    // common to both.
    std::vector<VertexId> parents(number_of_vertices);
    std::iota(parents.begin(), parents.end(), 0);
    auto find = [&parents](VertexId vertex_id)
    {
        while (parents[vertex_id] != vertex_id) {
            parents[vertex_id] = parents[parents[vertex_id]];
            vertex_id = parents[vertex_id];
        }
        return vertex_id;
    };
    for (VertexId vertex_id = 0; vertex_id < number_of_vertices; ++vertex_id) {
        for (int side = 0; side < 2; ++side) {
            VertexId u = neighbors_1[vertex_id][side];
            if (!in_tour(neighbors_2, vertex_id, u))
                parents[find(vertex_id)] = find(u);
            VertexId v = neighbors_2[vertex_id][side];
            if (!in_tour(neighbors_1, vertex_id, v))
                parents[find(vertex_id)] = find(v);
        }
    }

    // Number of common edges leaving each component, and length of each
    // tour inside it.
    std::vector<VertexId> number_of_cut_edges(number_of_vertices, 0);
    std::vector<Distance> length_1(number_of_vertices, 0);
    std::vector<Distance> length_2(number_of_vertices, 0);
    for (VertexId vertex_id = 0; vertex_id < number_of_vertices; ++vertex_id) {
        VertexId root = find(vertex_id);
        // Each edge counted from its tail (the vertex it leaves, following
        // the tour).
        VertexId next_1 = neighbors_1[vertex_id][1];
        if (find(next_1) != root) {
            number_of_cut_edges[root]++;
            number_of_cut_edges[find(next_1)]++;
        } else {
            length_1[root] += distances.distance(vertex_id, next_1);
        }
        VertexId next_2 = neighbors_2[vertex_id][1];
        if (find(next_2) == root)
            length_2[root] += distances.distance(vertex_id, next_2);
    }

    // Choose, for each component, the tour whose edges are used inside it.
    std::vector<bool> use_tour_2(number_of_vertices, false);
    bool any = false;
    for (VertexId vertex_id = 0; vertex_id < number_of_vertices; ++vertex_id) {
        if (find(vertex_id) != vertex_id)
            continue;
        if ((number_of_cut_edges[vertex_id] == 2 || number_of_cut_edges[vertex_id] == 0)
                && length_2[vertex_id] < length_1[vertex_id]) {
            use_tour_2[vertex_id] = true;
            any = true;
        }
    }
    if (!any)
        return tour_1;

    // Build the result, and check it is a single cycle.
    std::vector<std::array<VertexId, 2>> neighbors(number_of_vertices);
    // The edges leaving a component are common to both tours, so a vertex
    // simply keeps its neighbors in the tour chosen for its component.
    for (VertexId vertex_id = 0; vertex_id < number_of_vertices; ++vertex_id) {
        neighbors[vertex_id] = (use_tour_2[find(vertex_id)])?
            neighbors_2[vertex_id]:
            neighbors_1[vertex_id];
    }
    std::vector<VertexId> tour;
    tour.reserve(number_of_vertices);
    std::vector<bool> visited(number_of_vertices, false);
    VertexId previous_vertex_id = -1;
    VertexId vertex_id = 0;
    for (VertexId pos = 0; pos < number_of_vertices; ++pos) {
        if (visited[vertex_id])
            return tour_1;
        visited[vertex_id] = true;
        tour.push_back(vertex_id);
        VertexId next_vertex_id = (neighbors[vertex_id][0] != previous_vertex_id)?
            neighbors[vertex_id][0]:
            neighbors[vertex_id][1];
        // Both edges of a vertex must lead to its two tour neighbors.
        if (neighbors[next_vertex_id][0] != vertex_id && neighbors[next_vertex_id][1] != vertex_id)
            return tour_1;
        previous_vertex_id = vertex_id;
        vertex_id = next_vertex_id;
    }
    if (vertex_id != 0)
        return tour_1;
    return tour;
}


/**
 * Improve 'tour' (in place) with the shorter sub-paths of 'other': a
 * sub-path of 'tour' and one of 'other' (in either direction) which contain
 * the same vertices and have the same two end vertices can replace one
 * another. Non-overlapping replacements found in one pass are applied.
 * Return 'true' iff 'tour' has been improved.
 *
 * Two sub-paths tour[p..q] and other[k..l] can only match if the offsets of
 * their end vertices are the same (forward: pos_other - pos_tour; backward:
 * pos_other + pos_tour), so positions are grouped by offset, and each
 * position is only compared with the next 'maximum_number_of_pairs'
 * positions of its group. Vertex sets are compared by hashing (sums of
 * random values), then exactly.
 */
struct PartialTranscriptionBuffers
{
    std::vector<VertexId> positions_other;
    std::vector<uint64_t> hashes_tour;
    std::vector<uint64_t> hashes_other;
    std::vector<Distance> lengths_tour;
    std::vector<Distance> lengths_other;
    std::vector<uint8_t> used;
    std::vector<uint8_t> marks;
    std::vector<VertexId> offsets;
    std::vector<VertexId> group_starts;
    std::vector<VertexId> group_positions;
    std::vector<std::tuple<VertexId, VertexId, VertexId, bool>> replacements;
};

template <typename Distances>
bool partial_transcription(
        const Distances& distances,
        std::vector<VertexId>& tour,
        const std::vector<VertexId>& other,
        std::mt19937_64& generator,
        const std::vector<uint64_t>& vertex_hashes,
        PartialTranscriptionBuffers& buffers,
        int maximum_number_of_pairs = 10)
{
    VertexId n = tour.size();
    std::vector<VertexId>& positions_other = buffers.positions_other;
    positions_other.resize(n);
    for (VertexId pos = 0; pos < n; ++pos)
        positions_other[other[pos]] = pos;
    // Prefix sums of hashes and lengths (linear, no wrap-around).
    std::vector<uint64_t>& hashes_tour = buffers.hashes_tour;
    std::vector<uint64_t>& hashes_other = buffers.hashes_other;
    std::vector<Distance>& lengths_tour = buffers.lengths_tour;
    std::vector<Distance>& lengths_other = buffers.lengths_other;
    hashes_tour.resize(n + 1);
    hashes_other.resize(n + 1);
    lengths_tour.resize(n);
    lengths_other.resize(n);
    hashes_tour[0] = 0;
    hashes_other[0] = 0;
    lengths_tour[0] = 0;
    lengths_other[0] = 0;
    for (VertexId pos = 0; pos < n; ++pos) {
        hashes_tour[pos + 1] = hashes_tour[pos] + vertex_hashes[tour[pos]];
        hashes_other[pos + 1] = hashes_other[pos] + vertex_hashes[other[pos]];
    }
    for (VertexId pos = 1; pos < n; ++pos) {
        Distance distance_tour = distances.distance(tour[pos - 1], tour[pos]);
        Distance distance_other = distances.distance(other[pos - 1], other[pos]);
        lengths_tour[pos] = lengths_tour[pos - 1] + distance_tour;
        lengths_other[pos] = lengths_other[pos - 1] + distance_other;
    }

    // Replacements: (p, q, k, forward).
    auto& replacements = buffers.replacements;
    replacements.clear();
    std::vector<uint8_t>& used = buffers.used;
    std::vector<uint8_t>& marks = buffers.marks;
    used.assign(n, 0);
    marks.assign(n, 0);
    std::vector<VertexId>& offsets = buffers.offsets;
    std::vector<VertexId>& group_starts = buffers.group_starts;
    std::vector<VertexId>& group_positions = buffers.group_positions;
    offsets.resize(n);
    group_positions.resize(n);
    for (bool forward: {true, false}) {
        // Group the positions by offset (counting sort, positions in
        // increasing order within a group).
        group_starts.assign(n + 1, 0);
        for (VertexId pos = 0; pos < n; ++pos) {
            VertexId offset = (forward)?
                positions_other[tour[pos]] - pos:
                positions_other[tour[pos]] + pos;
            if (offset < 0)
                offset += n;
            if (offset >= n)
                offset -= n;
            offsets[pos] = offset;
            group_starts[offset + 1]++;
        }
        for (VertexId offset = 0; offset < n; ++offset)
            group_starts[offset + 1] += group_starts[offset];
        for (VertexId pos = 0; pos < n; ++pos)
            group_positions[group_starts[offsets[pos]]++] = pos;
        // (group_starts[offset] is now the end of the group 'offset'.)
        for (VertexId offset = 0; offset < n; ++offset) {
            VertexId group_start = (offset == 0)? 0: group_starts[offset - 1];
            VertexId group_end = group_starts[offset];
            for (VertexId t = group_start; t < group_end; ++t) {
                for (int r = 1; r <= maximum_number_of_pairs && t + r < group_end; ++r) {
                    VertexId p = group_positions[t];
                    VertexId q = group_positions[t + r];
                    VertexId length = q - p;
                    if (length < 2)
                        continue;
                    VertexId k = positions_other[tour[p]];
                    VertexId l = (forward)? k + length: k - length;
                    if (l < 0 || l >= n)
                        continue;
                    uint64_t hash_tour = hashes_tour[q + 1] - hashes_tour[p];
                    uint64_t hash_other = (forward)?
                        hashes_other[l + 1] - hashes_other[k]:
                        hashes_other[k + 1] - hashes_other[l];
                    if (hash_tour != hash_other)
                        continue;
                    Distance length_tour = lengths_tour[q] - lengths_tour[p];
                    Distance length_other = (forward)?
                        lengths_other[l] - lengths_other[k]:
                        lengths_other[k] - lengths_other[l];
                    if (length_other >= length_tour)
                        continue;
                    bool free = true;
                    for (VertexId pos = p; pos <= q && free; ++pos)
                        if (used[pos])
                            free = false;
                    if (!free)
                        continue;
                    // Exact check of the vertex sets.
                    for (VertexId pos = p; pos <= q; ++pos)
                        marks[tour[pos]] = 1;
                    bool same = true;
                    VertexId first = (forward)? k: l;
                    VertexId last = (forward)? l: k;
                    for (VertexId pos = first; pos <= last; ++pos)
                        if (!marks[other[pos]])
                            same = false;
                    for (VertexId pos = p; pos <= q; ++pos)
                        marks[tour[pos]] = 0;
                    if (!same)
                        continue;
                    for (VertexId pos = p; pos <= q; ++pos)
                        used[pos] = 1;
                    replacements.emplace_back(p, q, k, forward);
                }
            }
        }
    }
    (void)generator;
    for (const auto& replacement: replacements) {
        VertexId p = std::get<0>(replacement);
        VertexId q = std::get<1>(replacement);
        VertexId k = std::get<2>(replacement);
        bool forward = std::get<3>(replacement);
        for (VertexId pos = p; pos <= q; ++pos)
            tour[pos] = (forward)? other[k + (pos - p)]: other[k - (pos - p)];
    }
    return !replacements.empty();
}

template <typename Distances>
bool partial_transcription(
        const Distances& distances,
        std::vector<VertexId>& tour,
        const std::vector<VertexId>& other,
        std::mt19937_64& generator,
        const std::vector<uint64_t>& vertex_hashes,
        int maximum_number_of_pairs = 10)
{
    PartialTranscriptionBuffers buffers;
    return partial_transcription(distances, tour, other, generator, vertex_hashes, buffers, maximum_number_of_pairs);
}

/**
 * Merge two tours by iterative partial transcription (Möbius, Freund, et
 * al., "Combinatorial optimization by iterative partial transcription",
 * 2007): each tour is repeatedly improved with the shorter sub-paths of the
 * other (see 'partial_transcription'), with the tours rotated randomly
 * between rounds so that sub-paths crossing their start are found too.
 * Returns the shorter of the two resulting tours.
 */
template <typename Distances>
std::vector<VertexId> iterative_partial_transcription(
        const Distances& distances,
        const std::vector<VertexId>& tour_1,
        const std::vector<VertexId>& tour_2,
        std::mt19937_64& generator,
        int maximum_number_of_rounds = 20)
{
    VertexId n = tour_1.size();
    if (n < 4)
        return tour_1;
    std::vector<uint64_t> vertex_hashes(n);
    for (VertexId vertex_id = 0; vertex_id < n; ++vertex_id)
        vertex_hashes[vertex_id] = generator();
    std::vector<VertexId> tour_a = tour_1;
    std::vector<VertexId> tour_b = tour_2;
    std::uniform_int_distribution<VertexId> distribution(0, n - 1);
    PartialTranscriptionBuffers buffers;
    int number_of_rounds_without_improvement = 0;
    for (int round = 0;
            round < maximum_number_of_rounds && number_of_rounds_without_improvement < 2;
            ++round) {
        std::rotate(tour_a.begin(), tour_a.begin() + distribution(generator), tour_a.end());
        std::rotate(tour_b.begin(), tour_b.begin() + distribution(generator), tour_b.end());
        bool improved_a = partial_transcription(distances, tour_a, tour_b, generator, vertex_hashes, buffers);
        bool improved_b = partial_transcription(distances, tour_b, tour_a, generator, vertex_hashes, buffers);
        if (improved_a || improved_b) {
            number_of_rounds_without_improvement = 0;
        } else {
            number_of_rounds_without_improvement++;
        }
    }
    auto length = [&distances, n](const std::vector<VertexId>& tour)
    {
        Distance total = 0;
        for (VertexId pos = 0; pos < n; ++pos)
            total += distances.distance(tour[pos], tour[(pos + 1) % n]);
        return total;
    };
    return (length(tour_a) <= length(tour_b))? tour_a: tour_b;
}

}
