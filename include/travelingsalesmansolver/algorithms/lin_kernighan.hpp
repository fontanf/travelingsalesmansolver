#pragma once

#include "travelingsalesmansolver/algorithm_formatter.hpp"
#include "travelingsalesmansolver/solution_builder.hpp"
#include "travelingsalesmansolver/lin_kernighan/two_level_list.hpp"
#include "travelingsalesmansolver/lin_kernighan/candidates.hpp"
#include "travelingsalesmansolver/lin_kernighan/greedy.hpp"
#include "travelingsalesmansolver/lin_kernighan/k_opt_move.hpp"
#include "travelingsalesmansolver/lin_kernighan/sequential_move_patterns.hpp"
#include "travelingsalesmansolver/lin_kernighan/tour_merging.hpp"

#include <algorithm>
#include <cassert>
#include <limits>
#include <array>
#include <random>
#include <unordered_set>
#include <vector>

namespace travelingsalesmansolver
{

struct LinKernighanParameters: Parameters
{
    /** Number of candidate edges per vertex (nearest neighbors). */
    VertexId number_of_candidates = 10;

    /**
     * Size of the steps of a chain: 2 for 2-opt steps (chained LK, with
     * backtracking at the first levels, see 'breadth'); 3 to 5 for steps
     * searching exhaustively for a move of up to that many edges (as LKH,
     * whose default is 5).
     */
    int move_type = 2;

    /** Maximum number of steps in a chain. */
    int maximum_depth = 50;

    /**
     * Number of alternatives tried at each of the first levels of a chain;
     * a single one at the following levels.
     */
    std::vector<int> breadth = {5, 3};

    /**
     * How the search continues once the tour is a local optimum:
     * - "kicks": perturb the current tour locally (a random segment swap)
     *   and re-optimize around it; keep the result if it isn't longer;
     * - "trials" (as LKH): build a new starting tour by a random walk
     *   following the best tour's edges where they are minimum spanning
     *   tree edges, re-optimize the vertices whose edges differ from the
     *   best tour, then merge the result with the best tour.
     */
    std::string perturbation = "kicks";

    /**
     * Restricted search (as LKH): the first edge removed by a chain must
     * not be in the minimum spanning forest of the candidate graph (first
     * local search), then not in the best tour (trials).
     */
    bool restricted_search = true;

    /**
     * Non-sequential moves (as LKH's 'GAIN23'): at each new local optimum,
     * search for an improving move made of a sequential 2- or 3-opt move of
     * positive gain splitting the tour into two cycles, followed by a 2-opt
     * move joining them.
     */
    bool non_sequential_moves = true;

    /**
     * Trials: the starting tour's walk follows the best tour's edges where
     * they are spanning forest edges (otherwise, as LKH with nearest
     * neighbor candidates: a random walk on candidate edges, which include
     * the best tour's).
     */
    bool trial_walk_follows_best_tour = false;

    /**
     * Maximum number of kicks or trials ('-1': no limit, until the time
     * limit).
     */
    int64_t maximum_number_of_kicks = -1;

    /** Maximum length of each of the two segments swapped by a kick. */
    VertexId kick_segment_length = 50;

    /** Seed. */
    int seed = 0;


    virtual nlohmann::json to_json() const override
    {
        nlohmann::json json = Parameters::to_json();
        json.merge_patch({
                {"NumberOfCandidates", number_of_candidates},
                {"MoveType", move_type},
                {"MaximumDepth", maximum_depth},
                {"Breadth", breadth},
                {"Perturbation", perturbation},
                {"RestrictedSearch", restricted_search},
                {"NonSequentialMoves", non_sequential_moves},
                {"MaximumNumberOfKicks", maximum_number_of_kicks},
                {"KickSegmentLength", kick_segment_length},
                {"Seed", seed},
                });
        return json;
    }

    virtual int format_width() const override { return 37; }

    virtual void format(std::ostream& os) const override
    {
        Parameters::format(os);
        int width = format_width();
        os
            << std::setw(width) << std::left << "Number of candidates: " << number_of_candidates << std::endl
            << std::setw(width) << std::left << "Move type: " << move_type << std::endl
            << std::setw(width) << std::left << "Maximum depth: " << maximum_depth << std::endl
            << std::setw(width) << std::left << "Perturbation: " << perturbation << std::endl
            << std::setw(width) << std::left << "Restricted search: " << restricted_search << std::endl
            << std::setw(width) << std::left << "Non-sequential moves: " << non_sequential_moves << std::endl
            << std::setw(width) << std::left << "Maximum number of kicks: " << maximum_number_of_kicks << std::endl
            << std::setw(width) << std::left << "Kick segment length: " << kick_segment_length << std::endl
            << std::setw(width) << std::left << "Seed: " << seed << std::endl
            ;
    }
};

const Output lin_kernighan(
        const Instance& instance,
        const LinKernighanParameters& parameters = {});

template <typename Distances>
const Output lin_kernighan(
        const Distances& distances,
        const Instance& instance,
        const LinKernighanParameters& parameters = {});

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/**
 * State of the chained Lin-Kernighan algorithm (see 'lin_kernighan()').
 */
template <typename Distances>
struct LinKernighanData
{
    LinKernighanData(
            const Distances& distances,
            const Instance& instance,
            const LinKernighanParameters& parameters,
            const std::vector<VertexId>& initial_tour):
        distances(distances),
        instance(instance),
        parameters(parameters),
        tour(initial_tour),
        generator(parameters.seed),
        in_queue(instance.number_of_vertices(), false)
    {
        VertexId number_of_vertices = initial_tour.size();
        for (VertexId pos = 0; pos < number_of_vertices; ++pos) {
            length += distances.distance(
                    initial_tour[pos],
                    initial_tour[(pos + 1) % number_of_vertices]);
        }
    }

    const Distances& distances;

    const Instance& instance;

    const LinKernighanParameters& parameters;

    /** Candidate lists. */
    CandidateLists candidates;

    /** Lengths of the candidate edges (parallel to 'candidates'). */
    std::vector<std::vector<Distance>> candidate_distances;

    /**
     * For each vertex, its three shortest candidate edges, shortest first, as
     * (length, other end) ((-, -1) if there are fewer): at most two of them
     * are tour edges.
     */
    std::vector<std::array<std::pair<Distance, VertexId>, 3>> cheapest_candidates;

    /** Current tour. */
    TwoLevelList tour;

    /** Length of the current tour. */
    Distance length = 0;

    /** Random values of the vertices, for 'edge_hash'. */
    std::vector<uint64_t> vertex_hashes;

    /** Hash of the current tour (XOR of the hashes of its edges). */
    uint64_t tour_hash = 0;

    /**
     * Hashes of the local optima reached so far: a local search reaching
     * one of them again stops, since it would end at the same tour.
     */
    std::unordered_set<uint64_t> local_optima;

    /**
     * 2-opt moves applied since the last checkpoint (t1, t2, t3, t4, see
     * 'TwoLevelList::two_opt_move'), to be able to undo them.
     */
    std::vector<std::array<VertexId, 4>> moves;

    /** Random number generator. */
    std::mt19937_64 generator;

    /** Queue of the active vertices, to be used as 't1'. */
    std::vector<VertexId> queue;

    /** Position of the first vertex of 'queue' still to be processed. */
    VertexId queue_start = 0;

    /** Is a vertex in the queue? */
    std::vector<bool> in_queue;

    /** Edges added by the current chain (they can't be removed again). */
    std::vector<std::array<VertexId, 2>> added_edges;

    /** Best tour length found by the current chain. */
    Distance chain_best_length = 0;

    /** Number of moves applied when the current chain found its best tour. */
    VertexId chain_best_number_of_moves = 0;

    /*
     * k-opt steps ('move_type' >= 3)
     */

    /** Analysis of sequential moves. */
    SequentialMove sequential_move;

    /** Vertices of the move being built. */
    std::array<VertexId, 2 * SequentialMove::maximum_k> t;

    /** Best non-improving valid move of the current step. */
    std::array<VertexId, 2 * SequentialMove::maximum_k> best_t;

    /** Gain of 'best_t' before its closing edge. */
    Distance best_gain = 0;

    /** Has a best non-improving move been found? */
    bool has_best = false;

    /** For each vertex, the number of times it appears in the move being built. */
    std::vector<uint8_t> in_move;

    /** Patterns of sequential moves, to enumerate them ('k_opt_search'). */
    const SequentialMovePatterns& patterns = SequentialMovePatterns::get();

    /** Whether the move being built goes along 'next' (t[1] = next(t[0])). */
    bool forward_is_next = true;

    /** Position of t[0] along the tour. */
    int64_t key_origin = 0;

    /** Upper bound of the positions along the tour. */
    int64_t key_size = 0;

    /**
     * Keys of the ends of the move being built: their positions along the
     * direction of the move, from t[0] (see 'key').
     */
    std::array<int64_t, 2 * SequentialMove::maximum_k> end_keys;

    /** Vertex from which the next non-sequential search starts. */
    VertexId non_sequential_start = 0;

    /** Edges added by the steps already applied in the chain (keys, see 'edge_key'). */
    std::vector<uint64_t> chain_added_edges;

    /**
     * Restricted search: for each vertex, the neighbors 'n' such that the
     * edge (vertex, n) can't be the first edge removed by a chain (empty:
     * no restriction).
     */
    std::vector<std::vector<VertexId>> restricted_neighbors;
};

/** Can (t1, t2) be the first edge removed by a chain? */
template <typename Distances>
inline bool is_restricted(
        const LinKernighanData<Distances>& data,
        VertexId t1,
        VertexId t2)
{
    if (data.restricted_neighbors.empty())
        return false;
    const std::vector<VertexId>& neighbors = data.restricted_neighbors[t1];
    return std::find(neighbors.begin(), neighbors.end(), t2) != neighbors.end();
}

/** Update 'cheapest_candidates' with a new candidate edge. */
template <typename Distances>
inline void add_cheapest_candidate(
        LinKernighanData<Distances>& data,
        VertexId from,
        VertexId to,
        Distance distance)
{
    std::array<std::pair<Distance, VertexId>, 3>& cheapest = data.cheapest_candidates[from];
    std::pair<Distance, VertexId> entry = {distance, to};
    for (int pos = 0; pos < 3; ++pos) {
        if (cheapest[pos].second == -1 || entry.first < cheapest[pos].first)
            std::swap(entry, cheapest[pos]);
        if (entry.second == -1)
            break;
    }
}

/** Add an edge to the candidate lists of its ends, in distance order, if missing. */
template <typename Distances>
void add_candidate_edge(
        LinKernighanData<Distances>& data,
        VertexId vertex_id_1,
        VertexId vertex_id_2)
{
    for (int side = 0; side < 2; ++side) {
        VertexId from = (side == 0)? vertex_id_1: vertex_id_2;
        VertexId to = (side == 0)? vertex_id_2: vertex_id_1;
        std::vector<VertexId>& candidates = data.candidates[from];
        if (std::find(candidates.begin(), candidates.end(), to) != candidates.end())
            continue;
        Distance distance = data.distances.distance(from, to);
        std::vector<Distance>& candidate_distances = data.candidate_distances[from];
        size_t pos = 0;
        while (pos < candidates.size() && candidate_distances[pos] <= distance)
            ++pos;
        candidates.insert(candidates.begin() + pos, to);
        candidate_distances.insert(candidate_distances.begin() + pos, distance);
        add_cheapest_candidate(data, from, to, distance);
    }
}

/**
 * Move the candidate edges of both the best tour and the previous best tour
 * to the front of the candidate lists (keeping the order otherwise), so that
 * they are tried first (as in LKH).
 */
template <typename Distances>
void reorder_candidates(
        LinKernighanData<Distances>& data,
        const std::vector<std::array<VertexId, 2>>& best_neighbors,
        const std::vector<std::array<VertexId, 2>>& previous_best_neighbors)
{
    VertexId number_of_vertices = data.candidates.size();
    for (VertexId vertex_id = 0; vertex_id < number_of_vertices; ++vertex_id) {
        std::vector<VertexId>& candidates = data.candidates[vertex_id];
        std::vector<Distance>& candidate_distances = data.candidate_distances[vertex_id];
        auto is_common = [&](VertexId neighbor_id)
        {
            const std::array<VertexId, 2>& best = best_neighbors[vertex_id];
            const std::array<VertexId, 2>& previous_best = previous_best_neighbors[vertex_id];
            return (neighbor_id == best[0] || neighbor_id == best[1])
                && (neighbor_id == previous_best[0] || neighbor_id == previous_best[1]);
        };
        // Stable partition of both lists, common edges first.
        size_t number_of_common = 0;
        for (size_t pos = 0; pos < candidates.size(); ++pos) {
            if (!is_common(candidates[pos]))
                continue;
            for (size_t pos_2 = pos; pos_2 > number_of_common; --pos_2) {
                std::swap(candidates[pos_2], candidates[pos_2 - 1]);
                std::swap(candidate_distances[pos_2], candidate_distances[pos_2 - 1]);
            }
            number_of_common++;
        }
    }
}

/** Hash of an undirected edge. */
template <typename Distances>
inline uint64_t edge_hash(
        const LinKernighanData<Distances>& data,
        VertexId vertex_id_1,
        VertexId vertex_id_2)
{
    return data.vertex_hashes[vertex_id_1] * data.vertex_hashes[vertex_id_2];
}

/** Recompute the hash of the current tour. */
template <typename Distances>
void compute_tour_hash(
        LinKernighanData<Distances>& data)
{
    data.tour_hash = 0;
    for (VertexId vertex_id = 0; vertex_id < data.tour.number_of_vertices(); ++vertex_id)
        data.tour_hash ^= edge_hash(data, vertex_id, data.tour.next(vertex_id));
}

/** Apply a 2-opt move, and record it. */
template <typename Distances>
inline void apply_move(
        LinKernighanData<Distances>& data,
        VertexId t1,
        VertexId t2,
        VertexId t3,
        VertexId t4)
{
    data.tour.two_opt_move(t1, t2, t3, t4);
    data.tour_hash ^= edge_hash(data, t1, t2) ^ edge_hash(data, t3, t4)
        ^ edge_hash(data, t2, t3) ^ edge_hash(data, t4, t1);
    data.length += data.distances.distance(t2, t3)
        + data.distances.distance(t4, t1)
        - data.distances.distance(t1, t2)
        - data.distances.distance(t3, t4);
    data.moves.push_back({t1, t2, t3, t4});
}

/** Undo the last moves, until only 'number_of_moves' remain. */
template <typename Distances>
inline void undo_moves(
        LinKernighanData<Distances>& data,
        VertexId number_of_moves)
{
    while ((VertexId)data.moves.size() > number_of_moves) {
        std::array<VertexId, 4> move = data.moves.back();
        data.moves.pop_back();
        // The move replaced (t1, t2), (t3, t4) by (t2, t3), (t4, t1); the
        // reverse move replaces (t1, t4), (t3, t2) by (t4, t3), (t2, t1).
        data.tour.two_opt_move(move[0], move[3], move[2], move[1]);
        data.tour_hash ^= edge_hash(data, move[0], move[1]) ^ edge_hash(data, move[2], move[3])
            ^ edge_hash(data, move[1], move[2]) ^ edge_hash(data, move[3], move[0]);
        data.length += data.distances.distance(move[0], move[1])
            + data.distances.distance(move[2], move[3])
            - data.distances.distance(move[1], move[2])
            - data.distances.distance(move[3], move[0]);
    }
}

/** Add a vertex to the queue of active vertices, if not already in it. */
template <typename Distances>
inline void activate(
        LinKernighanData<Distances>& data,
        VertexId vertex_id)
{
    if (data.in_queue[vertex_id])
        return;
    data.in_queue[vertex_id] = true;
    data.queue.push_back(vertex_id);
}

/** Return 'true' iff the edge (u, v) has been added by the current chain. */
template <typename Distances>
inline bool added_by_chain(
        const LinKernighanData<Distances>& data,
        VertexId vertex_id_1,
        VertexId vertex_id_2)
{
    for (const auto& edge: data.added_edges) {
        if ((edge[0] == vertex_id_1 && edge[1] == vertex_id_2)
                || (edge[0] == vertex_id_2 && edge[1] == vertex_id_1)) {
            return true;
        }
    }
    return false;
}

/**
 * Extend the chain from the edge (t1, t2) of the current tour, 'gain' being
 * the total length of the edges removed by the chain so far (including
 * (t1, t2)) minus the total length of the edges added.
 *
 * Each step removes (t1, t2) and an edge (t3, t4), and adds (t2, t3) and
 * (t4, t1) (a 2-opt move); the chain then continues from (t1, t4). A step is
 * only considered while the gain stays positive (the positive gain
 * criterion), and the alternatives are tried in decreasing order of
 * d(t3, t4) - d(t2, t3).
 *
 * The search stops as soon as a chain reaches a tour shorter than the one
 * it started from; the moves after the best tour of that chain are then
 * still applied (the caller undoes them).
 */
template <typename Distances>
void lin_kernighan_step(
        LinKernighanData<Distances>& data,
        Distance start_length,
        VertexId t1,
        VertexId t2,
        Distance gain,
        int depth)
{
    if (depth >= data.parameters.maximum_depth)
        return;
    bool forward = (data.tour.next(t1) == t2);

    // Alternatives: (score, t3, t4).
    std::vector<std::tuple<Distance, VertexId, VertexId>> alternatives;
    for (VertexId t3: data.candidates[t2]) {
        Distance d23 = data.distances.distance(t2, t3);
        if (gain - d23 <= 0)
            break;
        if (t3 == t1)
            continue;
        VertexId t4 = (forward)? data.tour.previous(t3): data.tour.next(t3);
        if (t4 == t2)
            continue;
        if (added_by_chain(data, t3, t4))
            continue;
        alternatives.emplace_back(data.distances.distance(t3, t4) - d23, t3, t4);
    }
    std::sort(
            alternatives.begin(),
            alternatives.end(),
            [](const std::tuple<Distance, VertexId, VertexId>& alternative_1,
                const std::tuple<Distance, VertexId, VertexId>& alternative_2)
            {
                return std::get<0>(alternative_1) > std::get<0>(alternative_2);
            });
    int breadth = (depth < (int)data.parameters.breadth.size())?
        data.parameters.breadth[depth]:
        1;
    if ((int)alternatives.size() > breadth)
        alternatives.resize(breadth);

    for (const auto& alternative: alternatives) {
        VertexId t3 = std::get<1>(alternative);
        VertexId t4 = std::get<2>(alternative);
        VertexId number_of_moves = data.moves.size();
        Distance d23 = data.distances.distance(t2, t3);
        Distance d34 = data.distances.distance(t3, t4);
        apply_move(data, t1, t2, t3, t4);
        data.added_edges.push_back({t2, t3});
        if (data.length < data.chain_best_length) {
            data.chain_best_length = data.length;
            data.chain_best_number_of_moves = data.moves.size();
        }
        lin_kernighan_step(data, start_length, t1, t4, gain - d23 + d34, depth + 1);
        if (data.chain_best_length < start_length)
            return;
        undo_moves(data, number_of_moves);
        data.added_edges.pop_back();
    }
}

/**
 * Try to improve the tour with a chain starting at 't1'. Return 'true' iff
 * the tour has been improved.
 */
template <typename Distances>
bool lin_kernighan_improve(
        LinKernighanData<Distances>& data,
        VertexId t1)
{
    Distance start_length = data.length;
    VertexId start_number_of_moves = data.moves.size();
    for (VertexId t2: {data.tour.next(t1), data.tour.previous(t1)}) {
        if (is_restricted(data, t1, t2))
            continue;
        data.added_edges.clear();
        data.chain_best_length = start_length;
        data.chain_best_number_of_moves = start_number_of_moves;
        lin_kernighan_step(data, start_length, t1, t2, data.distances.distance(t1, t2), 0);
        if (data.chain_best_length < start_length) {
            undo_moves(data, data.chain_best_number_of_moves);
            for (VertexId move_pos = start_number_of_moves;
                    move_pos < (VertexId)data.moves.size();
                    ++move_pos) {
                for (VertexId vertex_id: data.moves[move_pos])
                    activate(data, vertex_id);
            }
            return true;
        }
    }
    return false;
}

/** Key of an undirected edge. */
template <typename Distances>
inline uint64_t edge_key(
        const LinKernighanData<Distances>& data,
        VertexId vertex_id_1,
        VertexId vertex_id_2)
{
    VertexId n = data.tour.number_of_vertices();
    return (vertex_id_1 < vertex_id_2)?
        (uint64_t)vertex_id_1 * n + vertex_id_2:
        (uint64_t)vertex_id_2 * n + vertex_id_1;
}

/** Apply the move analyzed last by 'data.sequential_move', as recorded 2-opt moves. */
template <typename Distances>
void apply_sequential_move(
        LinKernighanData<Distances>& data)
{
    data.sequential_move.apply(
            data.tour,
            [&data](VertexId vertex_id_1, VertexId vertex_id_2)
            {
                // Reverse the path from 'vertex_id_1' to 'vertex_id_2'.
                if (vertex_id_1 == vertex_id_2)
                    return;
                VertexId before_id = data.tour.previous(vertex_id_1);
                VertexId after_id = data.tour.next(vertex_id_2);
                if (after_id == vertex_id_1)
                    return;
                apply_move(data, before_id, vertex_id_1, after_id, vertex_id_2);
            });
}

/**
 * Length of the shortest candidate edge at 'vertex_id' which isn't a tour
 * edge.
 */
template <typename Distances>
inline Distance cheapest_candidate_distance(
        const LinKernighanData<Distances>& data,
        VertexId vertex_id)
{
    VertexId next_vertex_id = data.tour.next(vertex_id);
    VertexId previous_vertex_id = data.tour.previous(vertex_id);
    for (const auto& entry: data.cheapest_candidates[vertex_id]) {
        if (entry.second == -1)
            return 0;
        if (entry.second != next_vertex_id && entry.second != previous_vertex_id)
            return entry.first;
    }
    return 0;
}

/**
 * Key of a vertex: its position along the direction of the move being built,
 * from t[0].
 */
template <typename Distances>
inline int64_t key(
        const LinKernighanData<Distances>& data,
        VertexId vertex_id)
{
    int64_t position = data.tour.position(vertex_id) - data.key_origin;
    if (!data.forward_is_next)
        position = -position;
    if (position < 0)
        position += data.key_size;
    return position;
}

/**
 * Search for a sequential move of up to k = 'move_type' edges, extending
 * the move t[0], ..., t[2 * level - 1] (whose removed minus added edges
 * total 'gain', and whose pattern is 'pattern_id'; see
 * 'SequentialMovePatterns').
 *
 * Every combination of candidates is tried, level by level, as long as the
 * partial gain is positive (the gain criterion), skipping the moves which
 * can't be closed validly anymore. The first valid move with a positive gain
 * is applied (return 'true'). Otherwise, the valid k-opt move with the
 * largest gain before its closing edge is recorded ('best_t'), to continue
 * the chain from, provided that the next step can still find a positive
 * gain and that its last removed edge wasn't added by the chain (as in
 * LKH; the other edges aren't restricted).
 *
 * The pattern of the move gives its validity, and the pattern of each
 * extension, from the gap of the new vertex among the ends of the move (as
 * LKH's case analysis): moves are only analyzed when applied.
 */
template <typename Distances>
bool k_opt_search(
        LinKernighanData<Distances>& data,
        int level,
        Distance gain,
        int pattern_id)
{
    int k = data.parameters.move_type;
    bool last_level = (level + 1 == k);
    int number_of_ends = 2 * level;
    const SequentialMovePatterns::Pattern& pattern = data.patterns.pattern(pattern_id);
    VertexId t_last = data.t[2 * level - 1];
    const std::vector<VertexId>& candidates = data.candidates[t_last];
    const std::vector<Distance>& candidate_distances = data.candidate_distances[t_last];
    for (size_t candidate_pos = 0; candidate_pos < candidates.size(); ++candidate_pos) {
        VertexId t_new = candidates[candidate_pos];
        // Added edge (t_last, t_new).
        Distance gain_1 = gain - candidate_distances[candidate_pos];
        // Candidate lists aren't necessarily sorted by distance (see
        // 'reorder_candidates').
        if (gain_1 <= 0)
            continue;
        if (t_new == data.tour.next(t_last) || t_new == data.tour.previous(t_last))
            continue;
        if (t_new == data.t[0])
            continue;
        // Gap of t_new among the ends of the move, and the possible sides of
        // t_next.
        int gap = 0;
        int sides = 0;
        int64_t key_new = 0;
        if (data.in_move[t_new] == 0) {
            key_new = key(data, t_new);
            int number_of_ends_before = 0;
            for (int end = 0; end < number_of_ends; ++end)
                number_of_ends_before += (data.end_keys[end] < key_new);
            gap = number_of_ends_before - 1;
            sides = 3;
        } else if (data.in_move[t_new] == 1) {
            // An end of a removed edge: its other tour edge can be removed.
            int end = 0;
            while (data.t[end] != t_new)
                ++end;
            gap = pattern.end_gaps[end];
            sides = 1 << pattern.end_sides[end];
            key_new = data.end_keys[end];
        } else {
            continue;
        }
        for (int side_pos = 0; side_pos < 2; ++side_pos) {
            // Side 0: t_next follows t_new along the direction of the move.
            // (t_next = next(t_new) is tried first.)
            int side = (data.forward_is_next)? side_pos: 1 - side_pos;
            if (!(sides & (1 << side)))
                continue;
            int next_pattern_id = pattern.next[gap][side];
            const SequentialMovePatterns::Pattern& next_pattern = data.patterns.pattern(next_pattern_id);
            // Removed edge (t_new, t_next).
            VertexId t_next = (side_pos == 0)?
                data.tour.next(t_new):
                data.tour.previous(t_new);
            if (t_next == data.t[0])
                continue;
            data.t[2 * level] = t_new;
            data.t[2 * level + 1] = t_next;
            assert(data.sequential_move.analyze(data.tour, data.t.data(), level + 1)
                    == next_pattern.closable);
            assert(data.sequential_move.count_cycles(data.tour, data.t.data(), level + 1) - 1
                    == next_pattern.number_of_cycles);
            // Skip the moves which can't be closed validly anymore.
            if (last_level) {
                if (!next_pattern.closable)
                    continue;
            } else if (next_pattern.number_of_cycles > k - level - 1) {
                continue;
            }
            Distance gain_2 = gain_1 + data.distances.distance(t_new, t_next);
            if (next_pattern.closable
                    && gain_2 - data.distances.distance(t_next, data.t[0]) > 0) {
                data.sequential_move.analyze(data.tour, data.t.data(), level + 1);
                apply_sequential_move(data);
                return true;
            }
            if (!last_level) {
                data.end_keys[2 * level] = key_new;
                data.end_keys[2 * level + 1] = key(data, t_next);
                data.in_move[t_new]++;
                data.in_move[t_next]++;
                bool found = k_opt_search(data, level + 1, gain_2, next_pattern_id);
                data.in_move[t_new]--;
                data.in_move[t_next]--;
                if (found)
                    return true;
            } else if ((!data.has_best || gain_2 > data.best_gain)
                    // The next step must be able to add an edge at
                    // 't_next' with a positive gain.
                    && gain_2 > cheapest_candidate_distance(data, t_next)
                    // The last removed edge must be an original tour edge
                    // (not added by the chain), so that the chain can't
                    // cycle.
                    && std::find(
                        data.chain_added_edges.begin(),
                        data.chain_added_edges.end(),
                        edge_key(data, t_new, t_next)) == data.chain_added_edges.end()) {
                data.has_best = true;
                data.best_gain = gain_2;
                data.best_t = data.t;
            }
        }
    }
    return false;
}

/**
 * Try to improve the tour with a chain of k-opt steps starting at 't1'.
 * Return 'true' iff the tour has been improved.
 *
 * Each step searches for a move of up to k edges ('k_opt_search'). If none
 * improves the tour, the best one is applied anyway and the chain continues
 * from its closing edge.
 */
template <typename Distances>
bool lin_kernighan_improve_k_opt(
        LinKernighanData<Distances>& data,
        VertexId t1)
{
    int k = data.parameters.move_type;
    VertexId start_number_of_moves = data.moves.size();
    for (VertexId t2: {data.tour.next(t1), data.tour.previous(t1)}) {
        if (is_restricted(data, t1, t2))
            continue;
        data.chain_added_edges.clear();
        VertexId t2_current = t2;
        Distance gain = data.distances.distance(t1, t2);
        bool improved = false;
        for (int step = 0; step < data.parameters.maximum_depth; ++step) {
            data.t[0] = t1;
            data.t[1] = t2_current;
            data.has_best = false;
            data.forward_is_next = (data.tour.next(t1) == t2_current);
            data.key_origin = data.tour.position(t1);
            data.key_size = data.tour.number_of_segments() * (data.tour.number_of_vertices() + 1);
            data.end_keys[0] = 0;
            data.end_keys[1] = key(data, t2_current);
            data.in_move[t1]++;
            data.in_move[t2_current]++;
            bool found = k_opt_search(data, 1, gain, SequentialMovePatterns::first_pattern_id);
            data.in_move[t1]--;
            data.in_move[t2_current]--;
            if (found) {
                improved = true;
                break;
            }
            if (!data.has_best)
                break;
            // Apply the best non-improving move, and continue from it.
            data.sequential_move.analyze(data.tour, data.best_t.data(), k);
            apply_sequential_move(data);
            for (int i = 0; i < k - 1; ++i)
                data.chain_added_edges.push_back(edge_key(data, data.best_t[2 * i + 1], data.best_t[2 * i + 2]));
            gain = data.best_gain;
            t2_current = data.best_t[2 * k - 1];
        }
        if (improved) {
            for (VertexId move_pos = start_number_of_moves;
                    move_pos < (VertexId)data.moves.size();
                    ++move_pos) {
                for (VertexId vertex_id: data.moves[move_pos])
                    activate(data, vertex_id);
            }
            return true;
        }
        undo_moves(data, start_number_of_moves);
    }
    return false;
}

/** Is (u, v) one of the first k removed edges of the move 't'? */
inline bool is_removed(
        const VertexId* t,
        int k,
        VertexId u,
        VertexId v)
{
    for (int i = 0; i < k; ++i)
        if ((t[2 * i] == u && t[2 * i + 1] == v) || (t[2 * i] == v && t[2 * i + 1] == u))
            return true;
    return false;
}

/**
 * Apply the move last analyzed by 'data.sequential_move', of gain 'gain',
 * and activate its ends.
 */
template <typename Distances>
void apply_improving_move(
        LinKernighanData<Distances>& data,
        const VertexId* t,
        int k,
        Distance gain)
{
    Distance length = data.length;
    apply_sequential_move(data);
    (void)length;
    (void)gain;
    assert(data.length == length - gain);
    for (int pos = 0; pos < 2 * k; ++pos)
        activate(data, t[pos]);
}

/**
 * Try to join the two cycles made by the sequential move t[0], ...,
 * t[2k - 1] (of gain 'gain', whose cycles have been computed by
 * 'count_cycles') into a shorter tour, with a 2-opt move: remove a tour edge
 * (s1, s2) of the smaller cycle, add a candidate edge (s2, s3) to the other
 * cycle, remove a tour edge (s3, s4) there, and add (s4, s1). Apply the
 * first improving one found.
 */
template <typename Distances>
bool join_cycles(
        LinKernighanData<Distances>& data,
        std::array<VertexId, 2 * SequentialMove::maximum_k>& t,
        int k,
        Distance gain)
{
    SequentialMove& move = data.sequential_move;
    const TwoLevelList& tour = data.tour;
    std::array<std::array<int, SequentialMove::maximum_k>, 2> cycle_segments;
    std::array<int, 2> number_of_cycle_segments = {0, 0};
    for (int segment_id = 0; segment_id < k; ++segment_id) {
        int cycle_id = move.segment_cycle(segment_id);
        cycle_segments[cycle_id][number_of_cycle_segments[cycle_id]++] = segment_id;
    }

    // Find the smaller cycle, by walking both at the same pace.
    std::array<int, 2> positions = {0, 0};
    std::array<VertexId, 2> vertex_ids = {
        move.segment_first_vertex(cycle_segments[0][0]),
        move.segment_first_vertex(cycle_segments[1][0])};
    int smaller_cycle_id = -1;
    while (smaller_cycle_id == -1) {
        for (int cycle_id = 0; cycle_id < 2; ++cycle_id) {
            int segment_id = cycle_segments[cycle_id][positions[cycle_id]];
            if (vertex_ids[cycle_id] != move.segment_last_vertex(segment_id)) {
                vertex_ids[cycle_id] = tour.next(vertex_ids[cycle_id]);
                continue;
            }
            positions[cycle_id]++;
            if (positions[cycle_id] == number_of_cycle_segments[cycle_id]) {
                smaller_cycle_id = cycle_id;
                break;
            }
            vertex_ids[cycle_id] = move.segment_first_vertex(
                    cycle_segments[cycle_id][positions[cycle_id]]);
        }
    }

    // Added edges of the whole move: those of the sequential move, then
    // (s2, s3) and (s4, s1).
    std::array<int, 2 * SequentialMove::maximum_k> partners;
    for (int end = 0; end < 2 * k; ++end)
        partners[end] = (end % 2 == 1)? (end + 1) % (2 * k): (end + 2 * k - 1) % (2 * k);
    partners[2 * k] = 2 * k + 3;
    partners[2 * k + 3] = 2 * k;
    partners[2 * k + 1] = 2 * k + 2;
    partners[2 * k + 2] = 2 * k + 1;

    for (int pos = 0; pos < number_of_cycle_segments[smaller_cycle_id]; ++pos) {
        int segment_id = cycle_segments[smaller_cycle_id][pos];
        VertexId last_vertex_id = move.segment_last_vertex(segment_id);
        for (VertexId vertex_id = move.segment_first_vertex(segment_id);
                vertex_id != last_vertex_id;
                vertex_id = tour.next(vertex_id)) {
            VertexId next_vertex_id = tour.next(vertex_id);
            for (int side = 0; side < 2; ++side) {
                VertexId s1 = (side == 0)? vertex_id: next_vertex_id;
                VertexId s2 = (side == 0)? next_vertex_id: vertex_id;
                Distance gain_1 = gain + data.distances.distance(s1, s2);
                const std::vector<VertexId>& candidates = data.candidates[s2];
                const std::vector<Distance>& candidate_distances = data.candidate_distances[s2];
                for (size_t candidate_pos = 0; candidate_pos < candidates.size(); ++candidate_pos) {
                    VertexId s3 = candidates[candidate_pos];
                    Distance gain_2 = gain_1 - candidate_distances[candidate_pos];
                    if (gain_2 <= 0)
                        continue;
                    if (is_removed(t.data(), k, s2, s3))
                        continue;
                    if (move.segment_cycle(move.segment(tour, s3)) == smaller_cycle_id)
                        continue;
                    for (VertexId s4: {tour.next(s3), tour.previous(s3)}) {
                        if (is_removed(t.data(), k, s3, s4))
                            continue;
                        Distance total_gain = gain_2
                            + data.distances.distance(s3, s4)
                            - data.distances.distance(s4, s1);
                        if (total_gain <= 0)
                            continue;
                        t[2 * k] = s1;
                        t[2 * k + 1] = s2;
                        t[2 * k + 2] = s3;
                        t[2 * k + 3] = s4;
                        move.set_added_partners(partners.data(), k + 2);
                        bool valid = move.analyze(tour, t.data(), k + 2);
                        if (valid)
                            apply_improving_move(data, t.data(), k + 2, total_gain);
                        move.clear_added_partners();
                        if (valid)
                            return true;
                        // Restore the sequential move's analysis.
                        move.start(tour, t.data());
                        move.count_cycles(tour, t.data(), k);
                    }
                }
            }
        }
    }
    return false;
}

/**
 * Search for an improving non-sequential move (as LKH's 'GAIN23'): a
 * sequential 2- or 3-opt move of positive gain which splits the tour into
 * two cycles, followed by a 2-opt move joining them ('join_cycles'). Apply
 * the first one found.
 */
template <typename Distances>
bool non_sequential_search(
        LinKernighanData<Distances>& data)
{
    SequentialMove& move = data.sequential_move;
    const TwoLevelList& tour = data.tour;
    std::array<VertexId, 2 * SequentialMove::maximum_k> t;
    // Resume from where the previous search stopped.
    VertexId number_of_vertices = tour.number_of_vertices();
    for (VertexId count = 0; count < number_of_vertices; ++count) {
        VertexId t1 = data.non_sequential_start;
        data.non_sequential_start = (t1 + 1) % number_of_vertices;
        for (VertexId t2: {tour.next(t1), tour.previous(t1)}) {
            Distance gain_0 = data.distances.distance(t1, t2);
            const std::vector<VertexId>& candidates_2 = data.candidates[t2];
            for (size_t candidate_pos_2 = 0; candidate_pos_2 < candidates_2.size(); ++candidate_pos_2) {
                VertexId t3 = candidates_2[candidate_pos_2];
                Distance gain_1 = gain_0 - data.candidate_distances[t2][candidate_pos_2];
                if (gain_1 <= 0)
                    continue;
                if (t3 == tour.next(t2) || t3 == tour.previous(t2))
                    continue;
                for (VertexId t4: {tour.next(t3), tour.previous(t3)}) {
                    if (t4 == t1)
                        continue;
                    Distance gain_2 = gain_1 + data.distances.distance(t3, t4);
                    t[0] = t1;
                    t[1] = t2;
                    t[2] = t3;
                    t[3] = t4;
                    Distance gain = gain_2 - data.distances.distance(t4, t1);
                    if (gain > 0) {
                        move.start(tour, t.data());
                        int number_of_cycles = move.count_cycles(tour, t.data(), 2);
                        if (number_of_cycles == 1) {
                            apply_improving_move(data, t.data(), 2, gain);
                            return true;
                        }
                        if (number_of_cycles == 2 && join_cycles(data, t, 2, gain))
                            return true;
                    }
                    const std::vector<VertexId>& candidates_4 = data.candidates[t4];
                    for (size_t candidate_pos_4 = 0; candidate_pos_4 < candidates_4.size(); ++candidate_pos_4) {
                        VertexId t5 = candidates_4[candidate_pos_4];
                        Distance gain_3 = gain_2 - data.candidate_distances[t4][candidate_pos_4];
                        if (gain_3 <= 0)
                            continue;
                        if (t5 == tour.next(t4) || t5 == tour.previous(t4))
                            continue;
                        for (VertexId t6: {tour.next(t5), tour.previous(t5)}) {
                            if (t6 == t1 || is_removed(t.data(), 2, t5, t6))
                                continue;
                            Distance gain = gain_3
                                + data.distances.distance(t5, t6)
                                - data.distances.distance(t6, t1);
                            if (gain <= 0)
                                continue;
                            t[4] = t5;
                            t[5] = t6;
                            move.start(tour, t.data());
                            int number_of_cycles = move.count_cycles(tour, t.data(), 3);
                            if (number_of_cycles == 1) {
                                apply_improving_move(data, t.data(), 3, gain);
                                return true;
                            }
                            if (number_of_cycles == 2 && join_cycles(data, t, 3, gain))
                                return true;
                        }
                    }
                }
            }
        }
    }
    return false;
}

/** Run Lin-Kernighan steps until no active vertex remains. */
template <typename Distances>
void lin_kernighan_local_search(
        LinKernighanData<Distances>& data)
{
    for (;;) {
        while (data.queue_start < (VertexId)data.queue.size()) {
            VertexId t1 = data.queue[data.queue_start];
            data.queue_start++;
            data.in_queue[t1] = false;
            bool improved = (data.parameters.move_type <= 2)?
                lin_kernighan_improve(data, t1):
                lin_kernighan_improve_k_opt(data, t1);
            if (improved) {
                activate(data, t1);
                // Back to a local optimum already reached: this search would
                // end there again.
                if (data.local_optima.count(data.tour_hash)) {
                    for (VertexId pos = data.queue_start; pos < (VertexId)data.queue.size(); ++pos)
                        data.in_queue[data.queue[pos]] = false;
                    data.queue_start = data.queue.size();
                    break;
                }
            }
            if (data.queue_start > 4096 && 2 * data.queue_start > (VertexId)data.queue.size()) {
                data.queue.erase(data.queue.begin(), data.queue.begin() + data.queue_start);
                data.queue_start = 0;
            }
        }
        data.queue.clear();
        data.queue_start = 0;
        // A new local optimum: try to escape it with a non-sequential move.
        if (data.local_optima.count(data.tour_hash))
            break;
        data.local_optima.insert(data.tour_hash);
        if (!data.parameters.non_sequential_moves
                || !non_sequential_search(data)) {
            break;
        }
    }
}

/**
 * Reverse the path between 'vertex_id_1' and 'vertex_id_2', of
 * 'number_of_vertices' vertices, whichever direction the tour goes through
 * it, as a recorded 2-opt move.
 */
template <typename Distances>
void reverse_path(
        LinKernighanData<Distances>& data,
        VertexId vertex_id_1,
        VertexId vertex_id_2,
        VertexId number_of_vertices)
{
    if (number_of_vertices <= 1)
        return;
    VertexId vertex_id = vertex_id_1;
    for (VertexId pos = 1; pos < number_of_vertices; ++pos)
        vertex_id = data.tour.next(vertex_id);
    if (vertex_id != vertex_id_2)
        std::swap(vertex_id_1, vertex_id_2);
    // The path goes from 'vertex_id_1' to 'vertex_id_2' following 'next'.
    VertexId before_id = data.tour.previous(vertex_id_1);
    VertexId after_id = data.tour.next(vertex_id_2);
    apply_move(data, before_id, vertex_id_1, after_id, vertex_id_2);
}

/**
 * Kick: swap two consecutive segments of random lengths, starting at a
 * random vertex (a segment insertion, which a chain of 2-opt moves can't
 * easily undo). The vertices at the ends of the segments are activated.
 */
template <typename Distances>
void lin_kernighan_kick(
        LinKernighanData<Distances>& data)
{
    VertexId number_of_vertices = data.tour.number_of_vertices();
    VertexId maximum_length = (std::min)(
            data.parameters.kick_segment_length,
            (number_of_vertices - 2) / 2);
    std::uniform_int_distribution<VertexId> distribution_vertex(0, number_of_vertices - 1);
    std::uniform_int_distribution<VertexId> distribution_length(1, maximum_length);
    VertexId length_1 = distribution_length(data.generator);
    VertexId length_2 = distribution_length(data.generator);

    // a1 [a2 ... b1] [b2 ... c1] c2 -> a1 [b2 ... c1] [a2 ... b1] c2
    VertexId a1 = distribution_vertex(data.generator);
    VertexId a2 = data.tour.next(a1);
    VertexId b1 = a2;
    for (VertexId pos = 1; pos < length_1; ++pos)
        b1 = data.tour.next(b1);
    VertexId b2 = data.tour.next(b1);
    VertexId c1 = b2;
    for (VertexId pos = 1; pos < length_2; ++pos)
        c1 = data.tour.next(c1);
    VertexId c2 = data.tour.next(c1);

    // Reverse both segments together, then each of them.
    reverse_path(data, a2, c1, length_1 + length_2);
    reverse_path(data, c1, b2, length_2);
    reverse_path(data, b1, a2, length_1);

    for (VertexId vertex_id: {a1, a2, b1, b2, c1, c2})
        activate(data, vertex_id);
}

/**
 * Starting tour of a trial: a random walk. From each vertex, the next one
 * is (among unvisited vertices):
 * - if 'follow_best_tour', one of its neighbors in the best tour, if that
 *   edge is also in the minimum spanning forest of the candidate graph;
 * - otherwise, a random candidate neighbor;
 * - otherwise, the next vertex along the best tour.
 */
template <typename Distances>
std::vector<VertexId> trial_start_tour(
        LinKernighanData<Distances>& data,
        const std::vector<VertexId>& best_tour,
        const std::vector<std::array<VertexId, 2>>& best_neighbors,
        const std::vector<std::vector<VertexId>>& spanning_forest,
        bool follow_best_tour)
{
    VertexId number_of_vertices = best_tour.size();
    std::vector<VertexId> best_positions(number_of_vertices);
    for (VertexId pos = 0; pos < number_of_vertices; ++pos)
        best_positions[best_tour[pos]] = pos;
    // Next (and previous) unvisited position along the best tour
    // (union-find over positions: a visited position points to the next
    // (previous) one).
    std::vector<VertexId> next_unvisited(number_of_vertices);
    std::iota(next_unvisited.begin(), next_unvisited.end(), 0);
    std::vector<VertexId> previous_unvisited = next_unvisited;
    auto find = [](std::vector<VertexId>& links, VertexId pos)
    {
        VertexId root = pos;
        while (links[root] != root)
            root = links[root];
        while (links[pos] != root) {
            VertexId next_pos = links[pos];
            links[pos] = root;
            pos = next_pos;
        }
        return root;
    };
    std::vector<bool> visited(number_of_vertices, false);
    std::vector<VertexId> tour;
    tour.reserve(number_of_vertices);
    std::vector<VertexId> options;
    VertexId vertex_id = std::uniform_int_distribution<VertexId>(0, number_of_vertices - 1)(data.generator);
    for (;;) {
        tour.push_back(vertex_id);
        visited[vertex_id] = true;
        VertexId pos = best_positions[vertex_id];
        next_unvisited[pos] = (pos + 1) % number_of_vertices;
        previous_unvisited[pos] = (pos + number_of_vertices - 1) % number_of_vertices;
        if ((VertexId)tour.size() == number_of_vertices)
            break;
        options.clear();
        if (follow_best_tour) {
            for (VertexId neighbor_id: best_neighbors[vertex_id]) {
                if (visited[neighbor_id])
                    continue;
                const std::vector<VertexId>& forest_neighbors = spanning_forest[vertex_id];
                if (std::find(forest_neighbors.begin(), forest_neighbors.end(), neighbor_id) != forest_neighbors.end())
                    options.push_back(neighbor_id);
            }
        }
        if (options.empty()) {
            for (VertexId candidate_id: data.candidates[vertex_id])
                if (!visited[candidate_id])
                    options.push_back(candidate_id);
        }
        if (options.empty()) {
            // The closest of the next and previous unvisited vertices along
            // the best tour.
            VertexId next_vertex_id = best_tour[find(next_unvisited, pos)];
            VertexId previous_vertex_id = best_tour[find(previous_unvisited, pos)];
            vertex_id = (data.distances.distance(vertex_id, previous_vertex_id)
                    < data.distances.distance(vertex_id, next_vertex_id))?
                previous_vertex_id: next_vertex_id;
        } else {
            vertex_id = options[std::uniform_int_distribution<VertexId>(0, options.size() - 1)(data.generator)];
        }
    }
    return tour;
}

/** Length of a tour given as a sequence of vertices. */
template <typename Distances>
Distance tour_length(
        const Distances& distances,
        const std::vector<VertexId>& tour)
{
    Distance length = 0;
    for (size_t pos = 0; pos < tour.size(); ++pos)
        length += distances.distance(tour[pos], tour[(pos + 1) % tour.size()]);
    return length;
}

/** Build the solution corresponding to the current tour. */
template <typename Distances>
Solution lin_kernighan_solution(
        const LinKernighanData<Distances>& data)
{
    SolutionBuilder solution_builder(data.instance);
    for (VertexId vertex_id: data.tour.tour(0))
        solution_builder.add_vertex(vertex_id);
    Solution solution = solution_builder.build();
    // The length maintained incrementally must match the tour.
    assert(solution.distance() == data.length);
    return solution;
}

/**
 * Chained Lin-Kernighan (Applegate, Cook & Rohe, "Chained Lin-Kernighan for
 * large traveling salesman problems", 2003), with 2-opt steps.
 *
 * - Candidate edges: the nearest neighbors of each vertex.
 * - Initial tour: greedy edge.
 * - Local search: from each active vertex t1, chains of 2-opt moves (see
 *   'lin_kernighan_step'), with an active-vertex queue (don't-look bits):
 *   the vertices touched by an improvement become active again.
 * - Kicks: when no active vertex remains, a random segment swap; the new
 *   tour is kept if it isn't longer than the current one, undone otherwise.
 */
template <typename Distances>
const Output lin_kernighan(
        const Distances& distances,
        const Instance& instance,
        const LinKernighanParameters& parameters)
{
    Output output(instance);
    AlgorithmFormatter algorithm_formatter(parameters, output);
    algorithm_formatter.start("Lin-Kernighan");
    algorithm_formatter.print_header();

    VertexId number_of_vertices = instance.number_of_vertices();
    CandidateLists candidates = nearest_neighbor_candidates(
            instance.distances(),
            parameters.number_of_candidates);
    std::vector<VertexId> initial_tour = greedy_edge_tour(distances, candidates);
    LinKernighanData<Distances> data(distances, instance, parameters, initial_tour);
    data.candidates = std::move(candidates);
    data.in_move.assign(number_of_vertices, 0);
    data.vertex_hashes.resize(number_of_vertices);
    for (VertexId vertex_id = 0; vertex_id < number_of_vertices; ++vertex_id)
        data.vertex_hashes[vertex_id] = data.generator() | 1;
    compute_tour_hash(data);
    data.candidate_distances.resize(number_of_vertices);
    data.cheapest_candidates.assign(number_of_vertices, {{{0, -1}, {0, -1}, {0, -1}}});
    for (VertexId vertex_id = 0; vertex_id < number_of_vertices; ++vertex_id) {
        for (VertexId candidate_id: data.candidates[vertex_id]) {
            Distance distance = distances.distance(vertex_id, candidate_id);
            data.candidate_distances[vertex_id].push_back(distance);
            add_cheapest_candidate(data, vertex_id, candidate_id, distance);
        }
    }
    algorithm_formatter.update_solution(lin_kernighan_solution(data), "greedy edge");
    if (number_of_vertices < 8) {
        algorithm_formatter.end();
        return output;
    }

    std::vector<std::vector<VertexId>> spanning_forest = minimum_spanning_forest(
            distances,
            data.candidates);
    if (parameters.restricted_search)
        data.restricted_neighbors = spanning_forest;
    for (VertexId vertex_id = 0; vertex_id < number_of_vertices; ++vertex_id)
        activate(data, vertex_id);
    lin_kernighan_local_search(data);
    data.moves.clear();
    data.restricted_neighbors.clear();
    Distance best_length = data.length;
    algorithm_formatter.update_solution(lin_kernighan_solution(data), "Lin-Kernighan");

    double last_report_time = parameters.timer.elapsed_time();
    bool improved_since_report = false;
    if (parameters.perturbation == "trials") {
        std::vector<VertexId> best_tour = data.tour.tour(0);
        std::vector<std::array<VertexId, 2>> best_neighbors(number_of_vertices, {-1, -1});
        std::vector<std::array<VertexId, 2>> previous_best_neighbors(number_of_vertices, {-1, -1});
        // When the best tour changes: its edges become candidate edges,
        // and (restricted search) can't be the first edge removed by a
        // chain.
        auto update_best_neighbors = [&data, &parameters, &best_tour, &best_neighbors, &previous_best_neighbors, number_of_vertices]()
        {
            previous_best_neighbors.swap(best_neighbors);
            if (parameters.restricted_search)
                data.restricted_neighbors.assign(number_of_vertices, {});
            for (VertexId pos = 0; pos < number_of_vertices; ++pos) {
                VertexId vertex_id = best_tour[pos];
                VertexId next_vertex_id = best_tour[(pos + 1) % number_of_vertices];
                best_neighbors[vertex_id] = {
                    best_tour[(pos + number_of_vertices - 1) % number_of_vertices],
                    next_vertex_id};
                add_candidate_edge(data, vertex_id, next_vertex_id);
                if (parameters.restricted_search) {
                    data.restricted_neighbors[vertex_id].push_back(next_vertex_id);
                    data.restricted_neighbors[next_vertex_id].push_back(vertex_id);
                }
            }
            reorder_candidates(data, best_neighbors, previous_best_neighbors);
        };
        update_best_neighbors();
        for (int64_t trial_id = 0;
                !parameters.timer.needs_to_end()
                && (parameters.maximum_number_of_kicks == -1
                    || trial_id < parameters.maximum_number_of_kicks);
                ++trial_id) {
            std::vector<VertexId> start_tour = trial_start_tour(
                    data,
                    best_tour,
                    best_neighbors,
                    spanning_forest,
                    parameters.trial_walk_follows_best_tour);
            data.tour = TwoLevelList(start_tour);
            data.length = tour_length(distances, start_tour);
            compute_tour_hash(data);
            data.moves.clear();
            // Only the vertices whose edges differ from the best tour.
            for (VertexId pos = 0; pos < number_of_vertices; ++pos) {
                VertexId vertex_id = start_tour[pos];
                VertexId next_vertex_id = start_tour[(pos + 1) % number_of_vertices];
                if (best_neighbors[vertex_id][0] != next_vertex_id
                        && best_neighbors[vertex_id][1] != next_vertex_id) {
                    activate(data, vertex_id);
                    activate(data, next_vertex_id);
                }
            }
            lin_kernighan_local_search(data);
            data.moves.clear();
            // Merge with the best tour (starting from the shorter one).
            std::vector<VertexId> trial_tour = data.tour.tour(0);
            std::vector<VertexId> merged_tour = iterative_partial_transcription(
                    distances,
                    best_tour,
                    trial_tour,
                    data.generator);
            Distance merged_length = tour_length(distances, merged_tour);
            if (merged_length < best_length) {
                best_tour = std::move(merged_tour);
                best_length = merged_length;
                update_best_neighbors();
                improved_since_report = true;
            }
            if (improved_since_report
                    && parameters.timer.elapsed_time() > last_report_time + 1) {
                data.tour = TwoLevelList(best_tour);
                data.length = best_length;
                algorithm_formatter.update_solution(
                        lin_kernighan_solution(data),
                        "trial " + std::to_string(trial_id));
                last_report_time = parameters.timer.elapsed_time();
                improved_since_report = false;
            }
        }
        data.tour = TwoLevelList(best_tour);
        data.length = best_length;
    } else {
        for (int64_t kick_id = 0;
                !parameters.timer.needs_to_end()
                && (parameters.maximum_number_of_kicks == -1
                    || kick_id < parameters.maximum_number_of_kicks);
                ++kick_id) {
            lin_kernighan_kick(data);
            lin_kernighan_local_search(data);
            if (data.length <= best_length) {
                if (data.length < best_length)
                    improved_since_report = true;
                best_length = data.length;
                data.moves.clear();
            } else {
                undo_moves(data, 0);
            }
            // Report improvements at most every second (building a solution
            // costs O(n)).
            if (improved_since_report
                    && parameters.timer.elapsed_time() > last_report_time + 1) {
                algorithm_formatter.update_solution(
                        lin_kernighan_solution(data),
                        "kick " + std::to_string(kick_id));
                last_report_time = parameters.timer.elapsed_time();
                improved_since_report = false;
            }
            // Segment sizes drift with splits: reset them from time to time.
            if (kick_id % 1000 == 999)
                data.tour.rebuild();
        }
    }
    algorithm_formatter.update_solution(lin_kernighan_solution(data), "end");

    algorithm_formatter.end();
    return output;
}

}
