#pragma once

#include "travelingsalesmansolver/lin_kernighan/two_level_list.hpp"
#include "travelingsalesmansolver/candidates/candidate_lists.hpp"
#include "travelingsalesmansolver/lin_kernighan/k_opt_move.hpp"
#include "travelingsalesmansolver/lin_kernighan/sequential_move_patterns.hpp"
#include "travelingsalesmansolver/lin_kernighan/tour_merging.hpp"

#include <algorithm>
#include <cassert>
#include <chrono>
#include <functional>
#include <limits>
#include <array>
#include <numeric>
#include <random>
#include <string>
#include <unordered_set>
#include <vector>

/**
 * Lin-Kernighan engine: an LKH-style Lin-Kernighan on a tour of vertices,
 * shared by the algorithms of different problems, each given by a "problem"
 * class. The engine owns the tour, the search (sequential moves of up to 5
 * edges enumerated with 'SequentialMovePatterns', chains of steps,
 * non-sequential moves, active vertices, local optima), the candidate lists,
 * and the trials (random walks or kicks, merged with the best tour). The
 * problem gives the costs of the edges, which guide the search, and how the
 * tours are evaluated.
 *
 * A problem class 'Problem' provides (the hooks are resolved at compile time):
 *
 * - 'static constexpr bool evaluates_moves': 'false' if the objective is the
 *   total cost of the edges of the tour (a move is applied iff its gain is
 *   positive; 'Problem::Objective' must then be 'Distance'); 'true' if every
 *   valid closed move is evaluated by the problem ('evaluate_closing'), and
 *   applied iff it improves the objective of the start of the chain.
 * - 'static constexpr bool tracks_direction': whether the direction of the
 *   tour matters (see 'Engine::flipped').
 * - 'static constexpr bool uses_weights': whether the directed weights of
 *   the ends of the moves ('Engine::end_weights') are needed, from the
 *   weights of the edges ('weight').
 * - 'static constexpr bool has_dynamic_candidates': whether some vertices
 *   have dynamic candidates ('dynamic_candidates').
 * - 'using Objective': the objective of a tour, with 'operator<' (better).
 * - 'VertexId number_of_vertices() const'.
 * - 'Distance cost(VertexId, VertexId) const': cost of an edge, guiding the
 *   search (gain criterion).
 * - 'Distance weight(VertexId, VertexId) const' (if 'uses_weights').
 * - 'bool dynamic_candidates(VertexId) const', and
 *   'void compute_dynamic_candidates(const Engine&, VertexId,
 *   std::vector<VertexId>&) const' (if 'has_dynamic_candidates'): the
 *   candidates of such vertices are computed when needed, with cost 0.
 * - 'bool candidate_edge(VertexId, VertexId) const': whether an edge of the
 *   best tour may become a candidate edge.
 * - If 'evaluates_moves': 'void start_step(Engine&)', 'Objective
 *   evaluate_closing(Engine&, const SequentialMovePatterns::Pattern&,
 *   bool& backward)' and 'Objective evaluate_current(Engine&)'.
 * - 'Objective evaluate_tour(const std::vector<VertexId>&) const': objective
 *   of a tour given in its direction.
 * - 'std::vector<VertexId> trial_start_tour(Engine&, const
 *   std::vector<VertexId>& best_tour)': starting tour of a random-walk
 *   trial ('random_walk_tour' is the default one).
 * - 'void report(const std::vector<VertexId>& tour, const std::string&
 *   comment)': a new best tour.
 */

namespace travelingsalesmansolver
{
namespace lin_kernighan_engine
{

struct EngineParameters
{
    /**
     * Size of the steps of a chain: moves of up to that many edges (3 to 5;
     * LKH's default is 5).
     */
    int move_type = 5;

    /**
     * Maximum number of steps in a chain (1: a step only applies improving
     * moves, as LKH's 'MAX_SWAPS = 0').
     */
    int maximum_depth = 50;

    /**
     * How the search continues once the tour is a local optimum:
     * - "walks" (as LKH): trials starting from a random walk
     *   ('Problem::trial_start_tour'), re-optimizing the vertices whose edges
     *   differ from the best tour, then merged with the best tour;
     * - "double-bridge" (as LKH's 'SPECIAL' settings): trials starting from
     *   the best tour perturbed by a random double bridge;
     * - "segment-swap": perturb the current tour locally (a random segment
     *   swap) and re-optimize around it; keep the result if it isn't worse.
     */
    std::string perturbation = "walks";

    /**
     * Restricted search (as LKH): the first edge removed by a chain must not
     * be in the minimum spanning forest of the candidate graph (first local
     * search), then (walks only) not in the best tour.
     */
    bool restricted_search = true;

    /**
     * Non-sequential moves (as LKH's 'GAIN23'; problems not evaluating moves
     * only): at each new local optimum, search for an improving move made
     * of a sequential 2- or 3-opt move of positive gain splitting the tour
     * into two cycles, followed by a 2-opt move joining them.
     */
    bool non_sequential_moves = false;

    /** Maximum number of trials or kicks ('-1': no limit). */
    int64_t maximum_number_of_trials = -1;

    /** Maximum length of each of the two segments swapped by a segment swap. */
    VertexId kick_segment_length = 50;

    /** Should the search stop (time limit)? Checked between trials. */
    std::function<bool()> needs_to_end = []() { return false; };
};

/** State of the engine. */
template <typename Problem>
struct Engine
{
    using Objective = typename Problem::Objective;

    Engine(
            const Problem& problem,
            std::mt19937_64& generator,
            const EngineParameters& parameters):
        problem(problem),
        parameters(parameters),
        tour(std::vector<VertexId>{0}),
        generator(generator)
    {
    }

    /**
     * Problem (a copy: one indirection less than a reference in the hot
     * loops).
     */
    Problem problem;

    const EngineParameters& parameters;

    /** Number of vertices. */
    VertexId number_of_vertices = 0;

    /** Candidate lists (empty for the vertices with dynamic candidates). */
    CandidateLists candidates;

    /** Costs of the candidate edges (parallel to 'candidates'). */
    std::vector<std::vector<Distance>> candidate_distances;

    /**
     * For each vertex, its three cheapest candidate edges, cheapest first, as
     * (cost, other end) ((-, -1) if there are fewer): at most two of them
     * are tour edges.
     */
    std::vector<std::array<std::pair<Distance, VertexId>, 3>> cheapest_candidates;

    /** Current tour. */
    TwoLevelList tour;

    /**
     * Whether the direction of the tour is the direction of 'previous' in
     * 'tour' ('tracks_direction').
     */
    bool flipped = false;

    /** Total cost of the edges of the current tour (if not 'evaluates_moves'). */
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

    /** Random number generator (the caller's). */
    std::mt19937_64& generator;

    /** Queue of the active vertices, to be used as 't1'. */
    std::vector<VertexId> queue;

    /** Position of the first vertex of 'queue' still to be processed. */
    VertexId queue_start = 0;

    /** Is a vertex in the queue? */
    std::vector<bool> in_queue;

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

    /*
     * Used by some problems only (kept at the end, not to slow down the
     * others).
     */

    /** Objective of the current tour (if 'evaluates_moves'). */
    Objective objective = Objective();

    /** Objective at the start of the current chain (if 'evaluates_moves'). */
    Objective chain_objective = Objective();

    /**
     * Weights of the ends of the move being built (if 'uses_weights'): the
     * weight of the path from t[0] to them, along the direction of the move
     * (see 'directed_weight').
     */
    std::array<Distance, 2 * SequentialMove::maximum_k> end_weights;

    /** Total weight of the tour (if 'uses_weights'). */
    Distance total_weight = 0;

    /** Buffers of dynamic candidates, for each level of 'k_opt_search'. */
    std::array<std::vector<VertexId>, SequentialMove::maximum_k> dynamic_candidates;
    std::array<std::vector<Distance>, SequentialMove::maximum_k> dynamic_candidate_distances;
};

/** The costs of a problem, as distances (for the tour merging functions). */
template <typename Problem>
struct CostDistances
{
    const Problem& problem;

    inline Distance distance(VertexId vertex_id_1, VertexId vertex_id_2) const
    {
        return problem.cost(vertex_id_1, vertex_id_2);
    }
};

/*
 * Candidates
 */

/** Does a vertex have dynamic candidates? */
template <typename Problem>
inline bool has_dynamic_candidates(
        const Engine<Problem>& engine,
        VertexId vertex_id)
{
    if constexpr (Problem::has_dynamic_candidates) {
        return engine.problem.dynamic_candidates(vertex_id);
    } else {
        (void)engine;
        (void)vertex_id;
        return false;
    }
}

/** Update 'cheapest_candidates' with a new candidate edge. */
template <typename Problem>
inline void add_cheapest_candidate(
        Engine<Problem>& engine,
        VertexId from,
        VertexId to,
        Distance distance)
{
    std::array<std::pair<Distance, VertexId>, 3>& cheapest = engine.cheapest_candidates[from];
    std::pair<Distance, VertexId> entry = {distance, to};
    for (int pos = 0; pos < 3; ++pos) {
        if (cheapest[pos].second == -1 || entry.first < cheapest[pos].first)
            std::swap(entry, cheapest[pos]);
        if (entry.second == -1)
            break;
    }
}

/** Add a candidate edge from 'from' to 'to', in cost order, if missing. */
template <typename Problem>
void add_candidate(
        Engine<Problem>& engine,
        VertexId from,
        VertexId to)
{
    std::vector<VertexId>& candidates = engine.candidates[from];
    if (std::find(candidates.begin(), candidates.end(), to) != candidates.end())
        return;
    Distance distance = engine.problem.cost(from, to);
    std::vector<Distance>& candidate_distances = engine.candidate_distances[from];
    size_t pos = 0;
    while (pos < candidates.size() && candidate_distances[pos] <= distance)
        ++pos;
    candidates.insert(candidates.begin() + pos, to);
    candidate_distances.insert(candidate_distances.begin() + pos, distance);
    add_cheapest_candidate(engine, from, to, distance);
}

/**
 * Add a candidate edge to the candidate lists of its ends (except for the
 * ends with dynamic candidates), if missing.
 */
template <typename Problem>
void add_candidate_edge(
        Engine<Problem>& engine,
        VertexId vertex_id_1,
        VertexId vertex_id_2)
{
    if (!has_dynamic_candidates(engine, vertex_id_1))
        add_candidate(engine, vertex_id_1, vertex_id_2);
    if (!has_dynamic_candidates(engine, vertex_id_2))
        add_candidate(engine, vertex_id_2, vertex_id_1);
}

/** Set the candidate lists (and their costs). */
template <typename Problem>
void set_candidates(
        Engine<Problem>& engine,
        CandidateLists candidates)
{
    VertexId number_of_vertices = engine.number_of_vertices;
    engine.candidates = std::move(candidates);
    engine.candidates.resize(number_of_vertices);
    engine.candidate_distances.assign(number_of_vertices, {});
    engine.cheapest_candidates.assign(number_of_vertices, {{{0, -1}, {0, -1}, {0, -1}}});
    for (VertexId vertex_id = 0; vertex_id < number_of_vertices; ++vertex_id) {
        for (VertexId candidate_id: engine.candidates[vertex_id]) {
            Distance distance = engine.problem.cost(vertex_id, candidate_id);
            engine.candidate_distances[vertex_id].push_back(distance);
            add_cheapest_candidate(engine, vertex_id, candidate_id, distance);
        }
    }
}

/**
 * Move the candidate edges of both the best tour and the previous best tour
 * to the front of the candidate lists (keeping the order otherwise), so that
 * they are tried first (as in LKH).
 */
template <typename Problem>
void reorder_candidates(
        Engine<Problem>& engine,
        const std::vector<std::array<VertexId, 2>>& best_neighbors,
        const std::vector<std::array<VertexId, 2>>& previous_best_neighbors)
{
    for (VertexId vertex_id = 0; vertex_id < engine.number_of_vertices; ++vertex_id) {
        std::vector<VertexId>& candidates = engine.candidates[vertex_id];
        std::vector<Distance>& candidate_distances = engine.candidate_distances[vertex_id];
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

/**
 * Cost of the cheapest candidate edge at a vertex which isn't a tour edge.
 */
template <typename Problem>
inline Distance cheapest_candidate_distance(
        const Engine<Problem>& engine,
        VertexId vertex_id)
{
    if (has_dynamic_candidates(engine, vertex_id))
        return 0;
    VertexId next_vertex_id = engine.tour.next(vertex_id);
    VertexId previous_vertex_id = engine.tour.previous(vertex_id);
    for (const auto& entry: engine.cheapest_candidates[vertex_id]) {
        if (entry.second == -1)
            return 0;
        if (entry.second != next_vertex_id && entry.second != previous_vertex_id)
            return entry.first;
    }
    return 0;
}

/*
 * Tour
 */

/** Hash of an undirected edge. */
template <typename Problem>
inline uint64_t edge_hash(
        const Engine<Problem>& engine,
        VertexId vertex_id_1,
        VertexId vertex_id_2)
{
    return engine.vertex_hashes[vertex_id_1] * engine.vertex_hashes[vertex_id_2];
}

/** Key of an undirected edge. */
template <typename Problem>
inline uint64_t edge_key(
        const Engine<Problem>& engine,
        VertexId vertex_id_1,
        VertexId vertex_id_2)
{
    return (vertex_id_1 < vertex_id_2)?
        (uint64_t)vertex_id_1 * engine.number_of_vertices + vertex_id_2:
        (uint64_t)vertex_id_2 * engine.number_of_vertices + vertex_id_1;
}

/** Total cost of the edges of a tour given as a sequence of vertices. */
template <typename Problem>
Distance tour_cost(
        const Problem& problem,
        const std::vector<VertexId>& tour)
{
    Distance length = 0;
    for (size_t pos = 0; pos < tour.size(); ++pos)
        length += problem.cost(tour[pos], tour[(pos + 1) % tour.size()]);
    return length;
}

/** Set the current tour, given as a sequence of vertices (in its direction). */
template <typename Problem>
void set_tour(
        Engine<Problem>& engine,
        const std::vector<VertexId>& tour)
{
    engine.tour = TwoLevelList(tour);
    if constexpr (Problem::uses_weights) {
        Problem& problem = engine.problem;
        engine.tour.set_weight([&problem](VertexId vertex_id_1, VertexId vertex_id_2)
                {
                    return problem.weight(vertex_id_1, vertex_id_2);
                });
    }
    engine.flipped = false;
    engine.tour_hash = 0;
    for (VertexId vertex_id = 0; vertex_id < engine.tour.number_of_vertices(); ++vertex_id)
        engine.tour_hash ^= edge_hash(engine, vertex_id, engine.tour.next(vertex_id));
    if constexpr (Problem::evaluates_moves) {
        engine.objective = engine.problem.evaluate_current(engine);
    } else {
        engine.length = tour_cost(engine.problem, tour);
    }
    engine.moves.clear();
}

/** Sequence of the vertices of the current tour, in its direction. */
template <typename Problem>
std::vector<VertexId> current_tour(
        const Engine<Problem>& engine)
{
    std::vector<VertexId> tour = engine.tour.tour(0);
    if (engine.flipped)
        std::reverse(tour.begin(), tour.end());
    return tour;
}

/** Objective of the current tour. */
template <typename Problem>
inline typename Problem::Objective current_objective(
        const Engine<Problem>& engine)
{
    if constexpr (Problem::evaluates_moves) {
        return engine.objective;
    } else {
        return engine.length;
    }
}

/** Apply a 2-opt move, and record it. */
template <typename Problem>
inline void apply_move(
        Engine<Problem>& engine,
        VertexId t1,
        VertexId t2,
        VertexId t3,
        VertexId t4)
{
    bool flipped = engine.tour.two_opt_move(t1, t2, t3, t4);
    if constexpr (Problem::tracks_direction) {
        if (flipped)
            engine.flipped = !engine.flipped;
    } else {
        (void)flipped;
    }
    engine.tour_hash ^= edge_hash(engine, t1, t2) ^ edge_hash(engine, t3, t4)
        ^ edge_hash(engine, t2, t3) ^ edge_hash(engine, t4, t1);
    if constexpr (!Problem::evaluates_moves) {
        const Problem& problem = engine.problem;
        engine.length += problem.cost(t2, t3)
            + problem.cost(t4, t1)
            - problem.cost(t1, t2)
            - problem.cost(t3, t4);
    }
    engine.moves.push_back({t1, t2, t3, t4});
}

/** Undo the last moves, until only 'number_of_moves' remain. */
template <typename Problem>
inline void undo_moves(
        Engine<Problem>& engine,
        VertexId number_of_moves)
{
    while ((VertexId)engine.moves.size() > number_of_moves) {
        std::array<VertexId, 4> move = engine.moves.back();
        engine.moves.pop_back();
        // The move replaced (t1, t2), (t3, t4) by (t2, t3), (t4, t1); the
        // reverse move replaces (t1, t4), (t3, t2) by (t4, t3), (t2, t1).
        bool flipped = engine.tour.two_opt_move(move[0], move[3], move[2], move[1]);
        if constexpr (Problem::tracks_direction) {
            if (flipped)
                engine.flipped = !engine.flipped;
        } else {
            (void)flipped;
        }
        engine.tour_hash ^= edge_hash(engine, move[0], move[1]) ^ edge_hash(engine, move[2], move[3])
            ^ edge_hash(engine, move[1], move[2]) ^ edge_hash(engine, move[3], move[0]);
        if constexpr (!Problem::evaluates_moves) {
            const Problem& problem = engine.problem;
            engine.length += problem.cost(move[0], move[1])
                + problem.cost(move[2], move[3])
                - problem.cost(move[1], move[2])
                - problem.cost(move[3], move[0]);
        }
    }
}

/**
 * Apply the move last analyzed by 'engine.sequential_move' (a sequence of
 * path reversals keeping the direction of its segment 0).
 */
template <typename Problem>
void apply_sequential_move(
        Engine<Problem>& engine)
{
    engine.sequential_move.apply(
            engine.tour,
            [&engine](VertexId vertex_id_1, VertexId vertex_id_2)
            {
                // Reverse the path from 'vertex_id_1' to 'vertex_id_2'.
                if (vertex_id_1 == vertex_id_2)
                    return;
                VertexId before_id = engine.tour.previous(vertex_id_1);
                VertexId after_id = engine.tour.next(vertex_id_2);
                if (after_id == vertex_id_1)
                    return;
                apply_move(engine, before_id, vertex_id_1, after_id, vertex_id_2);
            });
}

/*
 * Sequential moves
 */

/** Can (t1, t2) be the first edge removed by a chain? */
template <typename Problem>
inline bool is_restricted(
        const Engine<Problem>& engine,
        VertexId t1,
        VertexId t2)
{
    if (engine.restricted_neighbors.empty())
        return false;
    const std::vector<VertexId>& neighbors = engine.restricted_neighbors[t1];
    return std::find(neighbors.begin(), neighbors.end(), t2) != neighbors.end();
}

/** Add a vertex to the queue of active vertices, if not already in it. */
template <typename Problem>
inline void activate(
        Engine<Problem>& engine,
        VertexId vertex_id)
{
    if (engine.in_queue[vertex_id])
        return;
    engine.in_queue[vertex_id] = true;
    engine.queue.push_back(vertex_id);
}

/**
 * Key of a vertex: its position along the direction of the move being built,
 * from t[0].
 */
template <typename Problem>
inline int64_t key(
        const Engine<Problem>& engine,
        VertexId vertex_id)
{
    int64_t position = engine.tour.position(vertex_id) - engine.key_origin;
    if (!engine.forward_is_next)
        position = -position;
    if (position < 0)
        position += engine.key_size;
    return position;
}

/**
 * Weight of the path from t[0] to a vertex, along the direction of the move
 * being built (t[0] must be set). (Not modulo the total weight: the weight of
 * a path can be 0 or the total weight if some edges weigh 0.)
 */
template <typename Problem>
inline Distance directed_weight(
        const Engine<Problem>& engine,
        VertexId vertex_id)
{
    return (engine.forward_is_next)?
        engine.tour.path_weight(engine.t[0], vertex_id):
        engine.tour.path_weight(vertex_id, engine.t[0]);
}

/** Initialize the move built by a step from its first removed edge (t1, t2). */
template <typename Problem>
inline void start_step(
        Engine<Problem>& engine,
        VertexId t1,
        VertexId t2)
{
    engine.t[0] = t1;
    engine.t[1] = t2;
    engine.forward_is_next = (engine.tour.next(t1) == t2);
    engine.key_origin = engine.tour.position(t1);
    engine.key_size = engine.tour.positions_size();
    engine.end_keys[0] = 0;
    engine.end_keys[1] = key(engine, t2);
    if constexpr (Problem::uses_weights) {
        engine.total_weight = engine.tour.total_weight();
        engine.end_weights[0] = 0;
        engine.end_weights[1] = directed_weight(engine, t2);
    }
    if constexpr (Problem::evaluates_moves)
        engine.problem.start_step(engine);
}

/**
 * Search for a sequential move of up to k = 'move_type' edges, extending the
 * move t[0], ..., t[2 * level - 1] (whose removed minus added edges total
 * 'gain', and whose pattern is 'pattern_id'; see 'SequentialMovePatterns').
 *
 * Every combination of candidates is tried, level by level, as long as the
 * partial gain is positive (the gain criterion), skipping the moves which
 * can't be closed validly anymore. A valid closed move is applied (return
 * 'true') if its gain is positive (not 'evaluates_moves'), or if it improves
 * the objective of the start of the chain ('evaluates_moves'). Otherwise, the
 * valid k-opt move with the largest gain before its closing edge is recorded
 * ('best_t'), to continue the chain from, provided that the next step can
 * still find a positive gain and that its last removed edge wasn't added by
 * the chain (as in LKH; the other edges aren't restricted).
 *
 * The pattern of the move gives its validity, and the pattern of each
 * extension, from the gap of the new vertex among the ends of the move (as
 * LKH's case analysis): moves are only analyzed when applied.
 */
template <typename Problem>
bool k_opt_search(
        Engine<Problem>& engine,
        int level,
        Distance gain,
        int pattern_id)
{
    const Problem& problem = engine.problem;
    int k = engine.parameters.move_type;
    bool last_level = (level + 1 == k);
    int number_of_ends = 2 * level;
    const SequentialMovePatterns::Pattern& pattern = engine.patterns.pattern(pattern_id);
    VertexId t_last = engine.t[2 * level - 1];
    const std::vector<VertexId>* candidates = &engine.candidates[t_last];
    const std::vector<Distance>* candidate_distances = &engine.candidate_distances[t_last];
    if constexpr (Problem::has_dynamic_candidates) {
        if (problem.dynamic_candidates(t_last)) {
            problem.compute_dynamic_candidates(engine, t_last, engine.dynamic_candidates[level]);
            engine.dynamic_candidate_distances[level].assign(engine.dynamic_candidates[level].size(), 0);
            candidates = &engine.dynamic_candidates[level];
            candidate_distances = &engine.dynamic_candidate_distances[level];
        }
    }
    for (size_t candidate_pos = 0; candidate_pos < candidates->size(); ++candidate_pos) {
        VertexId t_new = (*candidates)[candidate_pos];
        // Added edge (t_last, t_new).
        Distance gain_1 = gain - (*candidate_distances)[candidate_pos];
        // Candidate lists aren't necessarily sorted by cost (see
        // 'reorder_candidates').
        if (gain_1 <= 0)
            continue;
        if (t_new == engine.tour.next(t_last) || t_new == engine.tour.previous(t_last))
            continue;
        if (t_new == engine.t[0])
            continue;
        // Gap of t_new among the ends of the move, and the possible sides of
        // t_next.
        int gap = 0;
        int sides = 0;
        int64_t key_new = 0;
        Distance weight_new = 0;
        if (engine.in_move[t_new] == 0) {
            key_new = key(engine, t_new);
            if constexpr (Problem::uses_weights)
                weight_new = directed_weight(engine, t_new);
            int number_of_ends_before = 0;
            for (int end = 0; end < number_of_ends; ++end)
                number_of_ends_before += (engine.end_keys[end] < key_new);
            gap = number_of_ends_before - 1;
            sides = 3;
        } else if (engine.in_move[t_new] == 1) {
            // An end of a removed edge: its other tour edge can be removed.
            int end = 0;
            while (engine.t[end] != t_new)
                ++end;
            gap = pattern.end_gaps[end];
            sides = 1 << pattern.end_sides[end];
            key_new = engine.end_keys[end];
            if constexpr (Problem::uses_weights)
                weight_new = engine.end_weights[end];
        } else {
            continue;
        }
        (void)weight_new;
        for (int side_pos = 0; side_pos < 2; ++side_pos) {
            // Side 0: t_next follows t_new along the direction of the move.
            // (t_next = next(t_new) is tried first.)
            int side = (engine.forward_is_next)? side_pos: 1 - side_pos;
            if (!(sides & (1 << side)))
                continue;
            int next_pattern_id = pattern.next[gap][side];
            const SequentialMovePatterns::Pattern& next_pattern = engine.patterns.pattern(next_pattern_id);
            // Skip the moves which can't be closed validly anymore.
            if (last_level) {
                if (!next_pattern.closable)
                    continue;
            } else if (next_pattern.number_of_cycles > k - level - 1) {
                continue;
            }
            // Removed edge (t_new, t_next).
            VertexId t_next = (side_pos == 0)?
                engine.tour.next(t_new):
                engine.tour.previous(t_new);
            if (t_next == engine.t[0])
                continue;
            engine.t[2 * level] = t_new;
            engine.t[2 * level + 1] = t_next;
            assert(engine.sequential_move.analyze(engine.tour, engine.t.data(), level + 1)
                    == next_pattern.closable);
            Distance gain_2 = gain_1 + problem.cost(t_new, t_next);
            if constexpr (Problem::evaluates_moves) {
                if (next_pattern.closable) {
                    engine.end_keys[2 * level] = key_new;
                    engine.end_keys[2 * level + 1] = key(engine, t_next);
                    if constexpr (Problem::uses_weights) {
                        engine.end_weights[2 * level] = weight_new;
                        engine.end_weights[2 * level + 1] = directed_weight(engine, t_next);
                    }
                    bool backward = false;
                    typename Problem::Objective objective
                        = engine.problem.evaluate_closing(engine, next_pattern, backward);
                    if (objective < engine.chain_objective) {
                        bool valid = engine.sequential_move.analyze(engine.tour, engine.t.data(), level + 1);
                        (void)valid;
                        assert(valid);
                        apply_sequential_move(engine);
                        if constexpr (Problem::tracks_direction) {
                            // Direction of the new tour: t[0] follows
                            // t[2m - 1] iff not 'backward'.
                            bool t0_follows = (engine.tour.next(t_next) == engine.t[0]);
                            engine.flipped = (t0_follows == backward);
                        }
                        engine.objective = objective;
                        assert(!(engine.objective < engine.problem.evaluate_current(engine))
                                && !(engine.problem.evaluate_current(engine) < engine.objective));
                        return true;
                    }
                }
            } else {
                if (next_pattern.closable
                        && gain_2 - problem.cost(t_next, engine.t[0]) > 0) {
                    engine.sequential_move.analyze(engine.tour, engine.t.data(), level + 1);
                    apply_sequential_move(engine);
                    return true;
                }
            }
            if (!last_level) {
                engine.end_keys[2 * level] = key_new;
                engine.end_keys[2 * level + 1] = key(engine, t_next);
                if constexpr (Problem::uses_weights) {
                    engine.end_weights[2 * level] = weight_new;
                    engine.end_weights[2 * level + 1] = directed_weight(engine, t_next);
                }
                engine.in_move[t_new]++;
                engine.in_move[t_next]++;
                bool found = k_opt_search(engine, level + 1, gain_2, next_pattern_id);
                engine.in_move[t_new]--;
                engine.in_move[t_next]--;
                if (found)
                    return true;
            } else if ((!engine.has_best || gain_2 > engine.best_gain)
                    // The next step must be able to add an edge at
                    // 't_next' with a positive gain.
                    && gain_2 > cheapest_candidate_distance(engine, t_next)
                    // The last removed edge must be an original tour edge
                    // (not added by the chain), so that the chain can't
                    // cycle.
                    && std::find(
                        engine.chain_added_edges.begin(),
                        engine.chain_added_edges.end(),
                        edge_key(engine, t_new, t_next)) == engine.chain_added_edges.end()) {
                engine.has_best = true;
                engine.best_gain = gain_2;
                engine.best_t = engine.t;
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
 * from its closing edge (up to 'maximum_depth' steps).
 */
template <typename Problem>
bool improve(
        Engine<Problem>& engine,
        VertexId t1)
{
    int k = engine.parameters.move_type;
    VertexId start_number_of_moves = engine.moves.size();
    for (VertexId t2: {engine.tour.next(t1), engine.tour.previous(t1)}) {
        if (is_restricted(engine, t1, t2))
            continue;
        engine.chain_added_edges.clear();
        if constexpr (Problem::evaluates_moves)
            engine.chain_objective = engine.objective;
        VertexId t2_current = t2;
        Distance gain = engine.problem.cost(t1, t2);
        bool improved = false;
        for (int step = 0; step < engine.parameters.maximum_depth; ++step) {
            start_step(engine, t1, t2_current);
            engine.has_best = false;
            engine.in_move[t1]++;
            engine.in_move[t2_current]++;
            bool found = k_opt_search(engine, 1, gain, SequentialMovePatterns::first_pattern_id);
            engine.in_move[t1]--;
            engine.in_move[t2_current]--;
            if (found) {
                improved = true;
                break;
            }
            // (Not at the last step: the move would be undone.)
            if (!engine.has_best || step + 1 == engine.parameters.maximum_depth)
                break;
            // Apply the best non-improving move, and continue from it.
            engine.sequential_move.analyze(engine.tour, engine.best_t.data(), k);
            apply_sequential_move(engine);
            for (int i = 0; i < k - 1; ++i)
                engine.chain_added_edges.push_back(edge_key(engine, engine.best_t[2 * i + 1], engine.best_t[2 * i + 2]));
            gain = engine.best_gain;
            t2_current = engine.best_t[2 * k - 1];
        }
        if (improved) {
            for (VertexId move_pos = start_number_of_moves;
                    move_pos < (VertexId)engine.moves.size();
                    ++move_pos) {
                for (VertexId vertex_id: engine.moves[move_pos])
                    activate(engine, vertex_id);
            }
            return true;
        }
        undo_moves(engine, start_number_of_moves);
    }
    return false;
}

/*
 * Non-sequential moves
 */

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
 * Apply the move last analyzed by 'engine.sequential_move', of gain 'gain',
 * and activate its ends.
 */
template <typename Problem>
void apply_improving_move(
        Engine<Problem>& engine,
        const VertexId* t,
        int k,
        Distance gain)
{
    Distance length = engine.length;
    apply_sequential_move(engine);
    (void)length;
    (void)gain;
    assert(engine.length == length - gain);
    for (int pos = 0; pos < 2 * k; ++pos)
        activate(engine, t[pos]);
}

/**
 * Try to join the two cycles made by the sequential move t[0], ...,
 * t[2k - 1] (of gain 'gain', whose cycles have been computed by
 * 'count_cycles') into a cheaper tour, with a 2-opt move: remove a tour edge
 * (s1, s2) of the smaller cycle, add a candidate edge (s2, s3) to the other
 * cycle, remove a tour edge (s3, s4) there, and add (s4, s1). Apply the
 * first improving one found.
 */
template <typename Problem>
bool join_cycles(
        Engine<Problem>& engine,
        std::array<VertexId, 2 * SequentialMove::maximum_k>& t,
        int k,
        Distance gain)
{
    const Problem& problem = engine.problem;
    SequentialMove& move = engine.sequential_move;
    const TwoLevelList& tour = engine.tour;
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
                Distance gain_1 = gain + problem.cost(s1, s2);
                const std::vector<VertexId>& candidates = engine.candidates[s2];
                const std::vector<Distance>& candidate_distances = engine.candidate_distances[s2];
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
                            + problem.cost(s3, s4)
                            - problem.cost(s4, s1);
                        if (total_gain <= 0)
                            continue;
                        t[2 * k] = s1;
                        t[2 * k + 1] = s2;
                        t[2 * k + 2] = s3;
                        t[2 * k + 3] = s4;
                        move.set_added_partners(partners.data(), k + 2);
                        bool valid = move.analyze(tour, t.data(), k + 2);
                        if (valid)
                            apply_improving_move(engine, t.data(), k + 2, total_gain);
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
template <typename Problem>
bool non_sequential_search(
        Engine<Problem>& engine)
{
    const Problem& problem = engine.problem;
    SequentialMove& move = engine.sequential_move;
    const TwoLevelList& tour = engine.tour;
    std::array<VertexId, 2 * SequentialMove::maximum_k> t;
    // Resume from where the previous search stopped.
    VertexId number_of_vertices = tour.number_of_vertices();
    for (VertexId count = 0; count < number_of_vertices; ++count) {
        VertexId t1 = engine.non_sequential_start;
        engine.non_sequential_start = (t1 + 1) % number_of_vertices;
        for (VertexId t2: {tour.next(t1), tour.previous(t1)}) {
            Distance gain_0 = problem.cost(t1, t2);
            const std::vector<VertexId>& candidates_2 = engine.candidates[t2];
            for (size_t candidate_pos_2 = 0; candidate_pos_2 < candidates_2.size(); ++candidate_pos_2) {
                VertexId t3 = candidates_2[candidate_pos_2];
                Distance gain_1 = gain_0 - engine.candidate_distances[t2][candidate_pos_2];
                if (gain_1 <= 0)
                    continue;
                if (t3 == tour.next(t2) || t3 == tour.previous(t2))
                    continue;
                for (VertexId t4: {tour.next(t3), tour.previous(t3)}) {
                    if (t4 == t1)
                        continue;
                    Distance gain_2 = gain_1 + problem.cost(t3, t4);
                    t[0] = t1;
                    t[1] = t2;
                    t[2] = t3;
                    t[3] = t4;
                    Distance gain = gain_2 - problem.cost(t4, t1);
                    if (gain > 0) {
                        move.start(tour, t.data());
                        int number_of_cycles = move.count_cycles(tour, t.data(), 2);
                        if (number_of_cycles == 1) {
                            apply_improving_move(engine, t.data(), 2, gain);
                            return true;
                        }
                        if (number_of_cycles == 2 && join_cycles(engine, t, 2, gain))
                            return true;
                    }
                    const std::vector<VertexId>& candidates_4 = engine.candidates[t4];
                    for (size_t candidate_pos_4 = 0; candidate_pos_4 < candidates_4.size(); ++candidate_pos_4) {
                        VertexId t5 = candidates_4[candidate_pos_4];
                        Distance gain_3 = gain_2 - engine.candidate_distances[t4][candidate_pos_4];
                        if (gain_3 <= 0)
                            continue;
                        if (t5 == tour.next(t4) || t5 == tour.previous(t4))
                            continue;
                        for (VertexId t6: {tour.next(t5), tour.previous(t5)}) {
                            if (t6 == t1 || is_removed(t.data(), 2, t5, t6))
                                continue;
                            Distance gain = gain_3
                                + problem.cost(t5, t6)
                                - problem.cost(t6, t1);
                            if (gain <= 0)
                                continue;
                            t[4] = t5;
                            t[5] = t6;
                            move.start(tour, t.data());
                            int number_of_cycles = move.count_cycles(tour, t.data(), 3);
                            if (number_of_cycles == 1) {
                                apply_improving_move(engine, t.data(), 3, gain);
                                return true;
                            }
                            if (number_of_cycles == 2 && join_cycles(engine, t, 3, gain))
                                return true;
                        }
                    }
                }
            }
        }
    }
    return false;
}

/*
 * Local search
 */

/**
 * Run Lin-Kernighan steps until no active vertex remains; then, at a new
 * local optimum, try to escape it with a non-sequential move.
 */
template <typename Problem>
void local_search(
        Engine<Problem>& engine)
{
    for (;;) {
        while (engine.queue_start < (VertexId)engine.queue.size()) {
            VertexId t1 = engine.queue[engine.queue_start];
            engine.queue_start++;
            engine.in_queue[t1] = false;
            if (improve(engine, t1)) {
                activate(engine, t1);
                // Back to a local optimum already reached: this search would
                // end there again.
                if (engine.local_optima.count(engine.tour_hash)) {
                    for (VertexId pos = engine.queue_start; pos < (VertexId)engine.queue.size(); ++pos)
                        engine.in_queue[engine.queue[pos]] = false;
                    engine.queue_start = engine.queue.size();
                    break;
                }
            }
            if (engine.queue_start > 4096 && 2 * engine.queue_start > (VertexId)engine.queue.size()) {
                engine.queue.erase(engine.queue.begin(), engine.queue.begin() + engine.queue_start);
                engine.queue_start = 0;
            }
        }
        engine.queue.clear();
        engine.queue_start = 0;
        if (engine.local_optima.count(engine.tour_hash))
            break;
        engine.local_optima.insert(engine.tour_hash);
        if constexpr (Problem::evaluates_moves) {
            break;
        } else {
            if (!engine.parameters.non_sequential_moves
                    || !non_sequential_search(engine)) {
                break;
            }
        }
    }
}

/*
 * Perturbations
 */

/**
 * Reverse the path between 'vertex_id_1' and 'vertex_id_2', of
 * 'number_of_vertices' vertices, whichever direction the tour goes through
 * it, as a recorded 2-opt move.
 */
template <typename Problem>
void reverse_path(
        Engine<Problem>& engine,
        VertexId vertex_id_1,
        VertexId vertex_id_2,
        VertexId number_of_vertices)
{
    if (number_of_vertices <= 1)
        return;
    VertexId vertex_id = vertex_id_1;
    for (VertexId pos = 1; pos < number_of_vertices; ++pos)
        vertex_id = engine.tour.next(vertex_id);
    if (vertex_id != vertex_id_2)
        std::swap(vertex_id_1, vertex_id_2);
    // The path goes from 'vertex_id_1' to 'vertex_id_2' following 'next'.
    VertexId before_id = engine.tour.previous(vertex_id_1);
    VertexId after_id = engine.tour.next(vertex_id_2);
    apply_move(engine, before_id, vertex_id_1, after_id, vertex_id_2);
}

/**
 * Kick: swap two consecutive segments of random lengths, starting at a
 * random vertex (a segment insertion, which a chain of 2-opt moves can't
 * easily undo). The vertices at the ends of the segments are activated.
 */
template <typename Problem>
void segment_swap_kick(
        Engine<Problem>& engine)
{
    VertexId number_of_vertices = engine.tour.number_of_vertices();
    VertexId maximum_length = (std::min)(
            engine.parameters.kick_segment_length,
            (number_of_vertices - 2) / 2);
    std::uniform_int_distribution<VertexId> distribution_vertex(0, number_of_vertices - 1);
    std::uniform_int_distribution<VertexId> distribution_length(1, maximum_length);
    VertexId length_1 = distribution_length(engine.generator);
    VertexId length_2 = distribution_length(engine.generator);

    // a1 [a2 ... b1] [b2 ... c1] c2 -> a1 [b2 ... c1] [a2 ... b1] c2
    VertexId a1 = distribution_vertex(engine.generator);
    VertexId a2 = engine.tour.next(a1);
    VertexId b1 = a2;
    for (VertexId pos = 1; pos < length_1; ++pos)
        b1 = engine.tour.next(b1);
    VertexId b2 = engine.tour.next(b1);
    VertexId c1 = b2;
    for (VertexId pos = 1; pos < length_2; ++pos)
        c1 = engine.tour.next(c1);
    VertexId c2 = engine.tour.next(c1);

    // Reverse both segments together, then each of them.
    reverse_path(engine, a2, c1, length_1 + length_2);
    reverse_path(engine, c1, b2, length_2);
    reverse_path(engine, b1, a2, length_1);

    for (VertexId vertex_id: {a1, a2, b1, b2, c1, c2})
        activate(engine, vertex_id);
}

/**
 * Starting tour of a trial: a random walk, from each vertex to a random
 * unvisited candidate, else to the closest (in cost) of the next and
 * previous unvisited vertices along the best tour.
 */
template <typename Problem>
std::vector<VertexId> random_walk_tour(
        Engine<Problem>& engine,
        const std::vector<VertexId>& best_tour)
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
    VertexId vertex_id = std::uniform_int_distribution<VertexId>(0, number_of_vertices - 1)(engine.generator);
    for (;;) {
        tour.push_back(vertex_id);
        visited[vertex_id] = true;
        VertexId pos = best_positions[vertex_id];
        next_unvisited[pos] = (pos + 1) % number_of_vertices;
        previous_unvisited[pos] = (pos + number_of_vertices - 1) % number_of_vertices;
        if ((VertexId)tour.size() == number_of_vertices)
            break;
        options.clear();
        for (VertexId candidate_id: engine.candidates[vertex_id])
            if (!visited[candidate_id])
                options.push_back(candidate_id);
        if (options.empty()) {
            VertexId next_vertex_id = best_tour[find(next_unvisited, pos)];
            VertexId previous_vertex_id = best_tour[find(previous_unvisited, pos)];
            vertex_id = (engine.problem.cost(vertex_id, previous_vertex_id)
                    < engine.problem.cost(vertex_id, next_vertex_id))?
                previous_vertex_id: next_vertex_id;
        } else {
            vertex_id = options[std::uniform_int_distribution<VertexId>(0, options.size() - 1)(engine.generator)];
        }
    }
    return tour;
}

/**
 * Starting tour of a trial with kicks: the best tour, perturbed by a double
 * bridge on four random edges (as LKH's 'KICK_TYPE = 4').
 */
template <typename Problem>
std::vector<VertexId> double_bridge_tour(
        Engine<Problem>& engine,
        const std::vector<VertexId>& best_tour)
{
    VertexId number_of_vertices = best_tour.size();
    // Positions of the four removed edges (pos, pos + 1).
    std::array<VertexId, 4> positions;
    std::uniform_int_distribution<VertexId> distribution(0, number_of_vertices - 1);
    for (int i = 0; i < 4; ++i) {
        for (;;) {
            positions[i] = distribution(engine.generator);
            bool duplicate = false;
            for (int j = 0; j < i; ++j)
                if (positions[j] == positions[i])
                    duplicate = true;
            if (!duplicate)
                break;
        }
    }
    std::sort(positions.begin(), positions.end());
    // A = (p0 + 1 .. p1), B = (p1 + 1 .. p2), C = (p2 + 1 .. p3),
    // D = (p3 + 1 .. p0): A B C D -> A C B D.
    std::vector<VertexId> tour;
    tour.reserve(number_of_vertices);
    auto append = [&](VertexId pos_from, VertexId pos_to)
    {
        for (VertexId pos = pos_from; pos <= pos_to; ++pos)
            tour.push_back(best_tour[pos % number_of_vertices]);
    };
    append(positions[0] + 1, positions[1]);
    append(positions[2] + 1, positions[3]);
    append(positions[1] + 1, positions[2]);
    append(positions[3] + 1, positions[0] + number_of_vertices);
    return tour;
}

/*
 * Run
 */

/**
 * Run the Lin-Kernighan on a problem, from an initial tour and candidate
 * lists: a first local search (restricted by the minimum spanning forest of
 * the candidate graph), then trials or kicks ('perturbation') until
 * 'maximum_number_of_trials' or 'needs_to_end'. The best tours are reported
 * to the problem ('Problem::report', at most every second); the best tour is
 * returned. The engine works on a copy of 'input_problem'.
 */
template <typename Problem>
std::vector<VertexId> run(
        const Problem& input_problem,
        std::mt19937_64& generator,
        const EngineParameters& parameters,
        CandidateLists candidates,
        const std::vector<VertexId>& initial_tour)
{
    using Objective = typename Problem::Objective;
    Engine<Problem> engine(input_problem, generator, parameters);
    Problem& problem = engine.problem;
    VertexId number_of_vertices = problem.number_of_vertices();
    engine.number_of_vertices = number_of_vertices;
    engine.in_queue.assign(number_of_vertices, false);
    engine.in_move.assign(number_of_vertices, 0);
    engine.vertex_hashes.resize(number_of_vertices);
    for (VertexId vertex_id = 0; vertex_id < number_of_vertices; ++vertex_id)
        engine.vertex_hashes[vertex_id] = engine.generator() | 1;
    set_candidates(engine, std::move(candidates));
    set_tour(engine, initial_tour);
    problem.report(initial_tour, "initial tour");
    if (number_of_vertices < 8)
        return initial_tour;

    // First local search, restricted by the minimum spanning forest of the
    // candidate graph.
    CostDistances<Problem> costs{problem};
    if (parameters.restricted_search)
        engine.restricted_neighbors = minimum_spanning_forest(costs, engine.candidates);
    for (VertexId vertex_id = 0; vertex_id < number_of_vertices; ++vertex_id)
        activate(engine, vertex_id);
    local_search(engine);
    engine.moves.clear();
    engine.restricted_neighbors.clear();
    std::vector<VertexId> best_tour = current_tour(engine);
    Objective best_objective = current_objective(engine);
    problem.report(best_tour, "Lin-Kernighan");

    auto start = std::chrono::steady_clock::now();
    auto elapsed_time = [&start]()
    {
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    };
    double last_report_time = elapsed_time();
    bool improved_since_report = false;
    auto report = [&](const std::string& comment)
    {
        if (improved_since_report && elapsed_time() > last_report_time + 1) {
            problem.report(best_tour, comment);
            last_report_time = elapsed_time();
            improved_since_report = false;
        }
    };

    if (parameters.perturbation == "segment-swap") {
        for (int64_t kick_id = 0;
                !parameters.needs_to_end()
                && (parameters.maximum_number_of_trials == -1
                    || kick_id < parameters.maximum_number_of_trials);
                ++kick_id) {
            segment_swap_kick(engine);
            local_search(engine);
            Objective objective = current_objective(engine);
            if (!(best_objective < objective)) {
                if (objective < best_objective) {
                    best_objective = objective;
                    improved_since_report = true;
                }
                engine.moves.clear();
            } else {
                undo_moves(engine, 0);
                if constexpr (Problem::evaluates_moves)
                    engine.objective = best_objective;
            }
            if (improved_since_report
                    && elapsed_time() > last_report_time + 1) {
                best_tour = current_tour(engine);
            }
            report("kick " + std::to_string(kick_id));
            // Segment sizes drift with splits: reset them from time to time.
            if (kick_id % 1000 == 999)
                engine.tour.rebuild();
        }
        best_tour = current_tour(engine);
    } else {
        bool walks = (parameters.perturbation == "walks");
        std::vector<std::array<VertexId, 2>> best_neighbors(number_of_vertices, {-1, -1});
        std::vector<std::array<VertexId, 2>> previous_best_neighbors(number_of_vertices, {-1, -1});
        // When the best tour changes: its edges become candidate edges,
        // and (restricted search, walks) can't be the first edge removed by
        // a chain.
        auto update_best_neighbors = [&]()
        {
            previous_best_neighbors.swap(best_neighbors);
            bool restricted_search = parameters.restricted_search && walks;
            engine.restricted_neighbors.clear();
            if (restricted_search)
                engine.restricted_neighbors.assign(number_of_vertices, {});
            for (VertexId pos = 0; pos < number_of_vertices; ++pos) {
                VertexId vertex_id = best_tour[pos];
                VertexId next_vertex_id = best_tour[(pos + 1) % number_of_vertices];
                best_neighbors[vertex_id] = {
                    best_tour[(pos + number_of_vertices - 1) % number_of_vertices],
                    next_vertex_id};
                if (problem.candidate_edge(vertex_id, next_vertex_id))
                    add_candidate_edge(engine, vertex_id, next_vertex_id);
                if (restricted_search) {
                    engine.restricted_neighbors[vertex_id].push_back(next_vertex_id);
                    engine.restricted_neighbors[next_vertex_id].push_back(vertex_id);
                }
            }
            reorder_candidates(engine, best_neighbors, previous_best_neighbors);
        };
        update_best_neighbors();
        for (int64_t trial_id = 0;
                !parameters.needs_to_end()
                && (parameters.maximum_number_of_trials == -1
                    || trial_id < parameters.maximum_number_of_trials);
                ++trial_id) {
            std::vector<VertexId> start_tour = (walks)?
                problem.trial_start_tour(engine, best_tour):
                double_bridge_tour(engine, best_tour);
            set_tour(engine, start_tour);
            // Only the vertices whose edges differ from the best tour.
            for (VertexId pos = 0; pos < number_of_vertices; ++pos) {
                VertexId vertex_id = start_tour[pos];
                VertexId next_vertex_id = start_tour[(pos + 1) % number_of_vertices];
                if (best_neighbors[vertex_id][0] != next_vertex_id
                        && best_neighbors[vertex_id][1] != next_vertex_id) {
                    activate(engine, vertex_id);
                    activate(engine, next_vertex_id);
                }
            }
            local_search(engine);
            engine.moves.clear();
            // Merge with the best tour.
            std::vector<VertexId> trial_tour = current_tour(engine);
            std::vector<VertexId> merged_tour = iterative_partial_transcription(
                    costs,
                    best_tour,
                    trial_tour,
                    engine.generator);
            Objective merged_objective = problem.evaluate_tour(merged_tour);
            Objective trial_objective = current_objective(engine);
            if (trial_objective < merged_objective) {
                merged_tour = std::move(trial_tour);
                merged_objective = trial_objective;
            }
            if (merged_objective < best_objective) {
                best_tour = std::move(merged_tour);
                best_objective = merged_objective;
                update_best_neighbors();
                improved_since_report = true;
            }
            report("trial " + std::to_string(trial_id));
        }
    }
    problem.report(best_tour, "end");
    return best_tour;
}

}
}
