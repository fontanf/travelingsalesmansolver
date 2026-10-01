#include "travelingsalesmansolver/algorithms/eax.hpp"

#include <cassert>

using namespace travelingsalesmansolver;
using namespace travelingsalesmansolver::eax_internal;

namespace
{

/** Selection strategy of the AB-cycles of an E-set. */
enum class Strategy
{
    /** A single AB-cycle (localized EAX, stage I). */
    single,

    /** Block2 strategy (stage II). */
    block2,
};

struct NeighborList
{
    std::array<VertexId, 2> vertices;
    int size = 0;
};

template <typename Distances>
struct EaxData
{
    EaxData(
            const Distances& distances,
            const Instance& instance,
            std::mt19937_64& generator,
            const EaxParameters& parameters,
            AlgorithmFormatter& algorithm_formatter,
            const Output& output):
        distances(distances),
        instance(instance),
        generator(generator),
        parameters(parameters),
        algorithm_formatter(algorithm_formatter),
        output(output),
        number_of_vertices(instance.number_of_vertices()),
        edge_frequencies(number_of_vertices),
        remaining_a(number_of_vertices),
        remaining_b(number_of_vertices),
        trace_stamps(number_of_vertices, 0),
        occurrences(number_of_vertices),
        is_cut(number_of_vertices, 0),
        links(number_of_vertices, {-1, -1}),
        in_subtour(number_of_vertices, 0),
        cached_neighbors(number_of_vertices),
        cached_neighbor_distances(number_of_vertices),
        neighbor_stamps(number_of_vertices, 0),
        ineffective(number_of_vertices, 0),
        visits(number_of_vertices),
        is_linked(number_of_vertices, 0)
    {
    }

    const Distances& distances;

    const Instance& instance;

    std::mt19937_64& generator;

    const EaxParameters& parameters;

    AlgorithmFormatter& algorithm_formatter;

    const Output& output;

    VertexId number_of_vertices;

    /** Nearest neighbors, for the connection of the subtours, and their distances. */
    std::vector<std::vector<VertexId>> near_vertices;
    std::vector<std::vector<Distance>> near_distances;

    /**
     * Neighbors in the intermediate solution of the vertices near the
     * subtour being connected, and the lengths of their edges (valid if
     * 'neighbor_stamps' is 'subtour_stamp').
     */
    std::vector<std::array<VertexId, 2>> cached_neighbors;
    std::vector<std::array<Distance, 2>> cached_neighbor_distances;
    std::vector<int64_t> neighbor_stamps;

    std::vector<Tour> population;

    EdgeFrequencies edge_frequencies;

    /** Order of the vertices along pA ('Tour::order' of pA). */
    const int32_t* order = nullptr;

    /** Position of each vertex along pA ('Tour::positions' of pA). */
    const int32_t* positions = nullptr;

    /** Edges of pA (by position) in ineffective AB-cycles. */
    std::vector<int64_t> ineffective;
    int64_t ineffective_stamp = 0;

    /** Edges not yet in an AB-cycle. */
    std::vector<NeighborList> remaining_a;
    std::vector<NeighborList> remaining_b;

    /** Tracing: parents, and vertices reached (stamps). */
    const Tour* tracing_tour_a = nullptr;
    const Tour* tracing_tour_b = nullptr;
    std::vector<int64_t> trace_stamps;
    int64_t trace_stamp = 0;
    VertexId number_of_remaining_edges = 0;

    /** Traced path. */
    std::vector<VertexId> path;

    /** Indices of each vertex in the traced path. */
    std::vector<NeighborList> occurrences;

    /** Effective AB-cycles ('number_of_ab_cycles' first ones). */
    std::vector<std::vector<VertexId>> ab_cycles;
    VertexId number_of_ab_cycles = 0;

    /** Random bits, drawn 64 at a time. */
    uint64_t random_bits = 0;
    int number_of_random_bits = 0;

    /** Intermediate solution. */
    std::vector<char> is_cut;
    std::vector<VertexId> cuts;
    std::vector<std::array<VertexId, 2>> links;
    /** Vertices which have had links (each once). */
    std::vector<VertexId> linked_vertices;
    std::vector<char> is_linked;
    std::vector<VertexId> segment_subtours;
    std::vector<VertexId> subtour_sizes;

    /** Marks of the vertices of the subtour being connected. */
    std::vector<int64_t> in_subtour;
    int64_t subtour_stamp = 0;

    Offspring offspring;
    Offspring best_offspring;

    /** Block2 strategy. */
    std::vector<NeighborList> visits;
    std::vector<VertexId> visited_vertices;
    std::vector<int64_t> c_vertices;
    std::vector<std::vector<std::pair<VertexId, int64_t>>> contacts;
    std::vector<VertexId> ab_cycles_by_size;
    std::vector<char> in_e_set;
    std::vector<int64_t> contact_sums;
    std::vector<int64_t> tabu;
    std::vector<VertexId> candidates;
    std::vector<VertexId> candidate_positions;
    std::vector<VertexId> e_set_members;
    int64_t number_of_c_vertices = 0;
};

/*
 * Tours
 */

template <typename Distances>
Distance distance(
        const EaxData<Distances>& data,
        VertexId vertex_id_1,
        VertexId vertex_id_2)
{
    return data.distances.distance(vertex_id_1, vertex_id_2);
}

template <typename Distances>
Tour tour_from_solution(
        const EaxData<Distances>& data,
        const Solution& solution)
{
    Tour tour;
    tour.neighbors.resize(data.number_of_vertices);
    VertexId n = data.number_of_vertices;
    for (VertexId pos = 0; pos < n; ++pos) {
        VertexId vertex_id = solution.vertex_id(pos);
        tour.neighbors[vertex_id] = {
            solution.vertex_id((pos + n - 1) % n),
            solution.vertex_id((pos + 1) % n)};
    }
    tour.length = solution.distance();
    return tour;
}

/**
 * Make 'order' and 'positions' the order of the vertices along a tour
 * (computed if the tour has changed since the last time).
 */
template <typename Distances>
void compute_order(
        EaxData<Distances>& data,
        Tour& tour)
{
    if (!tour.order_valid) {
        tour.order.resize(data.number_of_vertices);
        tour.positions.resize(data.number_of_vertices);
        VertexId previous_vertex_id = tour.neighbors[0][0];
        VertexId vertex_id = 0;
        for (VertexId pos = 0; pos < data.number_of_vertices; ++pos) {
            tour.order[pos] = vertex_id;
            tour.positions[vertex_id] = pos;
            VertexId next_vertex_id = (tour.neighbors[vertex_id][0] != previous_vertex_id)?
                tour.neighbors[vertex_id][0]:
                tour.neighbors[vertex_id][1];
            previous_vertex_id = vertex_id;
            vertex_id = next_vertex_id;
        }
        tour.order_valid = true;
    }
    data.order = tour.order.data();
    data.positions = tour.positions.data();
}

template <typename Distances>
Solution solution_from_tour(
        EaxData<Distances>& data,
        const Tour& tour)
{
    SolutionBuilder solution_builder(data.instance);
    VertexId previous_vertex_id = tour.neighbors[0][0];
    VertexId vertex_id = 0;
    for (VertexId pos = 0; pos < data.number_of_vertices; ++pos) {
        solution_builder.add_vertex(vertex_id);
        VertexId next_vertex_id = (tour.neighbors[vertex_id][0] != previous_vertex_id)?
            tour.neighbors[vertex_id][0]:
            tour.neighbors[vertex_id][1];
        previous_vertex_id = vertex_id;
        vertex_id = next_vertex_id;
    }
    return solution_builder.build();
}

/*
 * Initial population
 */

template <typename Distances>
void initialize_population(EaxData<Distances>& data)
{
    // Nearest neighbors, for the 2-opt and the merging of the subtours.
    VertexId number_of_two_opt_candidates = 50;
    CandidateLists two_opt_candidates = nearest_neighbor_candidates(
            data.instance.distances(),
            number_of_two_opt_candidates);
    data.near_vertices.resize(data.number_of_vertices);
    data.near_distances.resize(data.number_of_vertices);
    for (VertexId vertex_id = 0; vertex_id < data.number_of_vertices; ++vertex_id) {
        const auto& list = two_opt_candidates[vertex_id];
        data.near_vertices[vertex_id] = list;
        data.near_distances[vertex_id].resize(list.size());
        for (size_t pos = 0; pos < list.size(); ++pos)
            data.near_distances[vertex_id][pos] = distance(data, vertex_id, list[pos]);
    }

    const std::string& initial_tours = data.parameters.initial_tours;
    const std::string& initial_local_search = data.parameters.initial_local_search;
    if (initial_tours != "random-permutation" && initial_tours != "random-walk") {
        throw std::invalid_argument(
                "travelingsalesmansolver::eax: "
                "unknown initial tours \"" + initial_tours + "\".");
    }
    if (initial_local_search != "2-opt"
            && initial_local_search != "lin-kernighan"
            && initial_local_search != "none") {
        throw std::invalid_argument(
                "travelingsalesmansolver::eax: "
                "unknown initial local search \"" + initial_local_search + "\".");
    }
    // Alpha-nearness candidates, for the random walks and the
    // Lin-Kernighan.
    CandidateLists alpha_nearness;
    if (initial_tours == "random-walk" || initial_local_search == "lin-kernighan") {
        AlphaNearnessParameters alpha_nearness_parameters;
        alpha_nearness_parameters.number_of_candidates = data.parameters.number_of_candidates;
        if (data.number_of_vertices >= 10000)
            alpha_nearness_parameters.initial_period = 100;
        alpha_nearness = alpha_nearness_candidates(
                data.instance.distances(),
                alpha_nearness_parameters).candidates;
    }

    for (int i = 0; i < data.parameters.population_size; ++i) {
        if (i > 0 && data.parameters.timer.needs_to_end())
            break;
        Solution solution(data.instance);
        if (initial_tours == "random-permutation") {
            RandomPermutationParameters random_permutation_parameters;
            random_permutation_parameters.verbosity_level = 0;
            solution = random_permutation(data.instance, data.generator, random_permutation_parameters).solution;
        } else {
            RandomWalkParameters random_walk_parameters;
            random_walk_parameters.verbosity_level = 0;
            solution = random_walk(data.instance, data.generator, random_walk_parameters, &alpha_nearness).solution;
        }
        if (initial_local_search == "2-opt") {
            TwoOptParameters two_opt_parameters;
            two_opt_parameters.verbosity_level = 0;
            solution = two_opt(data.distances, data.instance, data.generator, two_opt_parameters, &solution, &two_opt_candidates).solution;
        } else if (initial_local_search == "lin-kernighan") {
            LinKernighanParameters lin_kernighan_parameters;
            lin_kernighan_parameters.verbosity_level = 0;
            lin_kernighan_parameters.maximum_number_of_trials = 0;
            solution = lin_kernighan(data.distances, data.instance, data.generator, lin_kernighan_parameters, &solution, &alpha_nearness).solution;
        }
        data.population.push_back(tour_from_solution(data, solution));
        if (!data.output.solution.feasible()
                || solution.distance() < data.output.solution.distance()) {
            data.algorithm_formatter.update_solution(solution, "initial solution " + std::to_string(i));
        }
    }
    // Each edge is counted from its smaller end.
    for (const Tour& tour: data.population)
        for (VertexId vertex_id = 0; vertex_id < data.number_of_vertices; ++vertex_id)
            for (VertexId neighbor_id: tour.neighbors[vertex_id])
                if (vertex_id < neighbor_id)
                    data.edge_frequencies.add(vertex_id, neighbor_id, 1);
}

/*
 * AB-cycles
 */

/**
 * Compute the AB-cycles of parents A (whose order is in 'order') and B:
 * the edges of the multigraph G_AB (the edges of A and the edges of B,
 * common edges included twice) are partitioned into cycles alternating
 * edges of A and B, by random tracing.
 *
 * Each effective AB-cycle (with at least one edge of A not in B) is
 * stored in 'ab_cycles' as its vertices c[0], ..., c[2m - 1]: the edges
 * (c[2i], c[2i + 1]) are edges of A, the edges (c[2i + 1], c[2i + 2])
 * (indices modulo 2m) edges of B. The edges of A of the ineffective
 * AB-cycles (made of common edges) are marked in 'ineffective'.
 *
 * The tracing stops once 'maximum_number_of_ab_cycles' effective
 * AB-cycles have been formed ('-1': no limit). Each new tracing starting
 * from a random vertex, the AB-cycles formed are then biased towards the
 * largest ones.
 */
template <typename Distances>
void compute_ab_cycles(
        EaxData<Distances>& data,
        const Tour& tour_a,
        const Tour& tour_b,
        VertexId maximum_number_of_ab_cycles)
{
    data.number_of_ab_cycles = 0;
    data.ineffective_stamp++;
    // The remaining edges of the vertices are initialized when they're
    // first reached (the tracing may stop early).
    data.tracing_tour_a = &tour_a;
    data.tracing_tour_b = &tour_b;
    data.trace_stamp++;
    data.number_of_remaining_edges = 2 * data.number_of_vertices;

    std::vector<VertexId>& path = data.path;
    path.clear();
    while (data.number_of_remaining_edges > 0) {
        if (path.empty()) {
            // Random vertex with remaining edges.
            VertexId vertex_id = -1;
            std::uniform_int_distribution<VertexId> distribution(0, data.number_of_vertices - 1);
            for (int attempt = 0; attempt < 32 && vertex_id == -1; ++attempt) {
                VertexId candidate_id = distribution(data.generator);
                reach(data, candidate_id);
                if (data.remaining_a[candidate_id].size > 0)
                    vertex_id = candidate_id;
            }
            if (vertex_id == -1) {
                // Few vertices left: the first one from a random position.
                VertexId start = distribution(data.generator);
                for (VertexId offset = 0; offset < data.number_of_vertices; ++offset) {
                    VertexId candidate_id = (start + offset) % data.number_of_vertices;
                    reach(data, candidate_id);
                    if (data.remaining_a[candidate_id].size > 0) {
                        vertex_id = candidate_id;
                        break;
                    }
                }
            }
            path.push_back(vertex_id);
            data.occurrences[vertex_id].size = 0;
            data.occurrences[vertex_id].vertices[data.occurrences[vertex_id].size++] = 0;
        }
        VertexId vertex_id = path.back();
        // Edge index of the next edge: even for an edge of A.
        bool edge_a = ((path.size() - 1) % 2 == 0);
        NeighborList& remaining = (edge_a)? data.remaining_a[vertex_id]: data.remaining_b[vertex_id];
        if (remaining.size == 0) {
            // Only possible for a path reduced to its start vertex.
            clear_path(data);
            continue;
        }
        int pick = (remaining.size == 1)? 0: random_bit(data);
        VertexId next_vertex_id = remaining.vertices[pick];
        reach(data, next_vertex_id);
        remove_remaining_edge(data, vertex_id, next_vertex_id, edge_a);
        data.number_of_remaining_edges--;
        path.push_back(next_vertex_id);
        VertexId k = path.size() - 1;
        // Does the path close an AB-cycle at 'next_vertex_id'?
        NeighborList& occurrences = data.occurrences[next_vertex_id];
        bool closed = false;
        for (int o = 0; o < occurrences.size; ++o) {
            VertexId j = occurrences.vertices[o];
            if ((k - j) % 2 != 0)
                continue;
            // AB-cycle path[j], ..., path[k] (= path[j]); its edges of
            // A are the edges (path[i], path[i + 1]) for even i.
            VertexId first_a = (j % 2 == 0)? j: j + 1;
            bool effective = false;
            for (VertexId i = first_a; i < k; i += 2) {
                const auto& neighbors_b = tour_b.neighbors[path[i]];
                if (neighbors_b[0] != path[i + 1] && neighbors_b[1] != path[i + 1]) {
                    effective = true;
                    break;
                }
            }
            if (effective) {
                // Stored from an edge of A (vectors reused between pairs).
                if (data.number_of_ab_cycles == (VertexId)data.ab_cycles.size())
                    data.ab_cycles.emplace_back();
                std::vector<VertexId>& cycle = data.ab_cycles[data.number_of_ab_cycles++];
                cycle.clear();
                for (VertexId i = first_a; i < k; ++i)
                    cycle.push_back(path[i]);
                if (first_a != j)
                    cycle.push_back(path[j]);
            } else {
                for (VertexId i = first_a; i < k; i += 2)
                    data.ineffective[cut_position(data, path[i], path[i + 1])] = data.ineffective_stamp;
            }
            // Truncate the path to path[0], ..., path[j].
            for (VertexId i = j + 1; i < k; ++i)
                remove_occurrence(data, path[i], i);
            path.resize(j + 1);
            closed = true;
            break;
        }
        if (data.number_of_ab_cycles == maximum_number_of_ab_cycles && maximum_number_of_ab_cycles != -1)
            break;
        if (!closed)
            occurrences.vertices[occurrences.size++] = k;
    }
    clear_path(data);
}

template <typename Distances>
int random_bit(EaxData<Distances>& data)
{
    if (data.number_of_random_bits == 0) {
        data.random_bits = data.generator();
        data.number_of_random_bits = 64;
    }
    int bit = data.random_bits & 1;
    data.random_bits >>= 1;
    data.number_of_random_bits--;
    return bit;
}

template <typename Distances>
void remove_occurrence(
        EaxData<Distances>& data,
        VertexId vertex_id,
        VertexId index)
{
    NeighborList& occurrences = data.occurrences[vertex_id];
    for (int o = 0; o < occurrences.size; ++o) {
        if (occurrences.vertices[o] == index) {
            occurrences.vertices[o] = occurrences.vertices[occurrences.size - 1];
            occurrences.size--;
            return;
        }
    }
}

template <typename Distances>
void clear_path(EaxData<Distances>& data)
{
    for (VertexId vertex_id: data.path)
        data.occurrences[vertex_id].size = 0;
    data.path.clear();
}

template <typename Distances>
void remove_from_list(
        EaxData<Distances>& data,
        NeighborList& list,
        VertexId vertex_id)
{
    for (int i = 0; i < list.size; ++i) {
        if (list.vertices[i] == vertex_id) {
            list.vertices[i] = list.vertices[list.size - 1];
            list.size--;
            return;
        }
    }
}

template <typename Distances>
void remove_remaining_edge(
        EaxData<Distances>& data,
        VertexId vertex_id_1,
        VertexId vertex_id_2,
        bool edge_a)
{
    if (edge_a) {
        remove_from_list(data, data.remaining_a[vertex_id_1], vertex_id_2);
        remove_from_list(data, data.remaining_a[vertex_id_2], vertex_id_1);
    } else {
        remove_from_list(data, data.remaining_b[vertex_id_1], vertex_id_2);
        remove_from_list(data, data.remaining_b[vertex_id_2], vertex_id_1);
    }
}

/** Initialize the remaining edges of a vertex, the first time it's reached. */
template <typename Distances>
void reach(
        EaxData<Distances>& data,
        VertexId vertex_id)
{
    if (data.trace_stamps[vertex_id] == data.trace_stamp)
        return;
    data.trace_stamps[vertex_id] = data.trace_stamp;
    data.remaining_a[vertex_id].vertices = data.tracing_tour_a->neighbors[vertex_id];
    data.remaining_a[vertex_id].size = 2;
    data.remaining_b[vertex_id].vertices = data.tracing_tour_b->neighbors[vertex_id];
    data.remaining_b[vertex_id].size = 2;
    data.occurrences[vertex_id].size = 0;
}

/*
 * Intermediate solutions
 *
 * An intermediate solution is pA (whose order along the tour is
 * 'order') whose edges at the cut positions are removed ('cuts': the
 * edge between the positions c and c + 1), and with additional edges
 * ('links') between the ends of the resulting paths (the segments).
 */

template <typename Distances>
void reset_intermediate(EaxData<Distances>& data)
{
    for (VertexId cut: data.cuts)
        data.is_cut[cut] = 0;
    data.cuts.clear();
    for (VertexId vertex_id: data.linked_vertices) {
        data.links[vertex_id] = {-1, -1};
        data.is_linked[vertex_id] = 0;
    }
    data.linked_vertices.clear();
}

/** Position of the edge between consecutive vertices of pA (-1 if not consecutive). */
template <typename Distances>
VertexId cut_position(
        const EaxData<Distances>& data,
        VertexId vertex_id_1,
        VertexId vertex_id_2)
{
    VertexId pos_1 = data.positions[vertex_id_1];
    VertexId pos_2 = data.positions[vertex_id_2];
    if (pos_1 + 1 == pos_2 || (pos_1 == data.number_of_vertices - 1 && pos_2 == 0))
        return pos_1;
    if (pos_2 + 1 == pos_1 || (pos_2 == data.number_of_vertices - 1 && pos_1 == 0))
        return pos_2;
    return -1;
}

template <typename Distances>
void add_link(
        EaxData<Distances>& data,
        VertexId vertex_id_1,
        VertexId vertex_id_2)
{
    for (auto [a, b]: {std::pair{vertex_id_1, vertex_id_2}, std::pair{vertex_id_2, vertex_id_1}}) {
        auto& links = data.links[a];
        if (!data.is_linked[a]) {
            data.is_linked[a] = 1;
            data.linked_vertices.push_back(a);
        }
        if (links[0] == -1) {
            links[0] = b;
        } else {
            links[1] = b;
        }
    }
}

template <typename Distances>
void remove_link(
        EaxData<Distances>& data,
        VertexId vertex_id_1,
        VertexId vertex_id_2)
{
    for (auto [a, b]: {std::pair{vertex_id_1, vertex_id_2}, std::pair{vertex_id_2, vertex_id_1}}) {
        auto& links = data.links[a];
        if (links[0] == b) {
            links[0] = -1;
        } else {
            links[1] = -1;
        }
    }
}

/** Remove an edge of the intermediate solution. */
template <typename Distances>
void remove_edge(
        EaxData<Distances>& data,
        VertexId vertex_id_1,
        VertexId vertex_id_2)
{
    VertexId cut = cut_position(data, vertex_id_1, vertex_id_2);
    if (cut != -1 && !data.is_cut[cut]) {
        data.is_cut[cut] = 1;
        data.cuts.insert(std::upper_bound(data.cuts.begin(), data.cuts.end(), cut), cut);
    } else {
        remove_link(data, vertex_id_1, vertex_id_2);
    }
}

/** Get the two neighbors of a vertex in the intermediate solution. */
template <typename Distances>
std::array<VertexId, 2> neighbors(
        const EaxData<Distances>& data,
        VertexId vertex_id)
{
    std::array<VertexId, 2> result;
    int size = 0;
    VertexId pos = data.positions[vertex_id];
    VertexId previous_pos = (pos == 0)? data.number_of_vertices - 1: pos - 1;
    if (!data.is_cut[previous_pos])
        result[size++] = data.order[previous_pos];
    if (!data.is_cut[pos])
        result[size++] = data.order[(pos + 1 == data.number_of_vertices)? 0: pos + 1];
    for (VertexId link: data.links[vertex_id])
        if (link != -1 && size < 2)
            result[size++] = link;
    return result;
}

/** Segment containing a position: the index of the cut ending it. */
template <typename Distances>
VertexId segment(
        const EaxData<Distances>& data,
        VertexId pos)
{
    VertexId index = std::lower_bound(data.cuts.begin(), data.cuts.end(), pos) - data.cuts.begin();
    return (index == (VertexId)data.cuts.size())? 0: index;
}

template <typename Distances>
VertexId segment_start(
        const EaxData<Distances>& data,
        VertexId segment_id)
{
    VertexId pos = (segment_id == 0)? data.cuts.back() + 1: data.cuts[segment_id - 1] + 1;
    return (pos == data.number_of_vertices)? 0: pos;
}

template <typename Distances>
VertexId segment_end(
        const EaxData<Distances>& data,
        VertexId segment_id)
{
    return data.cuts[segment_id];
}

template <typename Distances>
VertexId segment_size(
        const EaxData<Distances>& data,
        VertexId segment_id)
{
    VertexId start = segment_start(data, segment_id);
    VertexId end = segment_end(data, segment_id);
    return (end >= start)? end - start + 1: end - start + 1 + data.number_of_vertices;
}

/** Compute the subtours of the intermediate solution. */
template <typename Distances>
void compute_subtours(EaxData<Distances>& data)
{
    VertexId number_of_segments = data.cuts.size();
    data.segment_subtours.assign(number_of_segments, -1);
    data.subtour_sizes.clear();
    for (VertexId segment_id_0 = 0; segment_id_0 < number_of_segments; ++segment_id_0) {
        if (data.segment_subtours[segment_id_0] != -1)
            continue;
        VertexId subtour_id = data.subtour_sizes.size();
        VertexId size = 0;
        VertexId segment_id = segment_id_0;
        VertexId entry_vertex_id = data.order[segment_start(data, segment_id)];
        VertexId previous_vertex_id = -1;
        for (;;) {
            data.segment_subtours[segment_id] = subtour_id;
            size += segment_size(data, segment_id);
            VertexId start_vertex_id = data.order[segment_start(data, segment_id)];
            VertexId end_vertex_id = data.order[segment_end(data, segment_id)];
            VertexId exit_vertex_id = (entry_vertex_id == start_vertex_id)?
                end_vertex_id:
                start_vertex_id;
            // Link leaving the segment.
            const auto& links = data.links[exit_vertex_id];
            VertexId next_vertex_id;
            if (start_vertex_id == end_vertex_id) {
                // Single vertex: leave by the link it wasn't entered by.
                next_vertex_id = (links[0] != previous_vertex_id)? links[0]: links[1];
            } else {
                next_vertex_id = (links[0] != -1)? links[0]: links[1];
            }
            previous_vertex_id = exit_vertex_id;
            entry_vertex_id = next_vertex_id;
            segment_id = segment(data, data.positions[next_vertex_id]);
            if (segment_id == segment_id_0)
                break;
        }
        data.subtour_sizes.push_back(size);
    }
}

/**
 * Connect the subtours of the intermediate solution into a tour (step 5
 * of EAX): repeatedly, the smallest subtour is connected to another one
 * by the best exchange of an edge of it and an edge of the other
 * subtour, at least one end of which is among the nearest neighbors of
 * an end of the first edge.
 */
template <typename Distances>
void merge_subtours(EaxData<Distances>& data)
{
    compute_subtours(data);
    while (data.subtour_sizes.size() > 1) {
        VertexId subtour_id = std::min_element(data.subtour_sizes.begin(), data.subtour_sizes.end()) - data.subtour_sizes.begin();
        Distance best_delta = std::numeric_limits<Distance>::max();
        std::array<VertexId, 4> best = {-1, -1, -1, -1};
        // Mark the vertices of the subtour.
        data.subtour_stamp++;
        for (VertexId segment_id = 0; segment_id < (VertexId)data.cuts.size(); ++segment_id) {
            if (data.segment_subtours[segment_id] != subtour_id)
                continue;
            VertexId pos = segment_start(data, segment_id);
            VertexId end = segment_end(data, segment_id);
            for (;;) {
                data.in_subtour[data.order[pos]] = data.subtour_stamp;
                if (pos == end)
                    break;
                pos = (pos + 1 == data.number_of_vertices)? 0: pos + 1;
            }
        }
        // The nearest neighbors first, then more of them if none is
        // outside of the subtour.
        for (VertexId near_limit: {data.parameters.number_of_near_vertices, (VertexId)50}) {
            for (VertexId segment_id = 0; segment_id < (VertexId)data.cuts.size(); ++segment_id) {
                if (data.segment_subtours[segment_id] != subtour_id)
                    continue;
                VertexId pos = segment_start(data, segment_id);
                VertexId end = segment_end(data, segment_id);
                for (;;) {
                    // Edges (u, u_1) and (u, u_2) of the subtour; edges (w,
                    // w_1) and (w, w_2) of another subtour, w being near u.
                    VertexId u = data.order[pos];
                    std::array<VertexId, 2> u_neighbors = neighbors(data, u);
                    Distance d_u_0 = distance(data, u, u_neighbors[0]);
                    Distance d_u_1 = distance(data, u, u_neighbors[1]);
                    const std::vector<VertexId>& near_vertices = data.near_vertices[u];
                    const std::vector<Distance>& near_distances = data.near_distances[u];
                    VertexId number_of_near_vertices = std::min(near_limit, (VertexId)near_vertices.size());
                    for (VertexId near_pos = 0; near_pos < number_of_near_vertices; ++near_pos) {
                        VertexId w = near_vertices[near_pos];
                        if (data.in_subtour[w] == data.subtour_stamp)
                            continue;
                        if (data.neighbor_stamps[w] != data.subtour_stamp) {
                            data.neighbor_stamps[w] = data.subtour_stamp;
                            data.cached_neighbors[w] = neighbors(data, w);
                            data.cached_neighbor_distances[w] = {
                                distance(data, w, data.cached_neighbors[w][0]),
                                distance(data, w, data.cached_neighbors[w][1])};
                        }
                        const std::array<VertexId, 2>& w_neighbors = data.cached_neighbors[w];
                        const std::array<Distance, 2>& w_distances = data.cached_neighbor_distances[w];
                        Distance d_uw = near_distances[near_pos];
                        for (int i = 0; i < 2; ++i) {
                            VertexId u_2 = u_neighbors[i];
                            Distance d_u = (i == 0)? d_u_0: d_u_1;
                            Distance d_u2_w = distance(data, u_2, w);
                            for (int j = 0; j < 2; ++j) {
                                VertexId w_2 = w_neighbors[j];
                                Distance d_removed = d_u + w_distances[j];
                                // Add (u, w) and (u_2, w_2). (The length of
                                // the second added edge is only computed
                                // if the move can be better than the best
                                // one, and the same below.)
                                if (best_delta > d_uw - d_removed) {
                                    Distance delta_1 = d_uw + distance(data, u_2, w_2) - d_removed;
                                    if (best_delta > delta_1) {
                                        best_delta = delta_1;
                                        best = {u, u_2, w, w_2};
                                    }
                                }
                                // Add (u, w_2) and (u_2, w).
                                if (best_delta > d_u2_w - d_removed) {
                                    Distance delta_2 = distance(data, u, w_2) + d_u2_w - d_removed;
                                    if (best_delta > delta_2) {
                                        best_delta = delta_2;
                                        best = {u, u_2, w_2, w};
                                    }
                                }
                            }
                        }
                    }
                    if (pos == end)
                        break;
                    pos = (pos + 1 == data.number_of_vertices)? 0: pos + 1;
                }
            }
            if (best[0] != -1)
                break;
        }
        if (best[0] == -1) {
            // No near vertex in another subtour: connect a vertex of the
            // subtour to the nearest vertex outside of it.
            VertexId vertex_id = -1;
            for (VertexId segment_id = 0; segment_id < (VertexId)data.cuts.size(); ++segment_id)
                if (data.segment_subtours[segment_id] == subtour_id)
                    vertex_id = data.order[segment_start(data, segment_id)];
            VertexId vertex_id_2 = neighbors(data, vertex_id)[0];
            for (VertexId vertex_id_3 = 0; vertex_id_3 < data.number_of_vertices; ++vertex_id_3) {
                if (data.in_subtour[vertex_id_3] == data.subtour_stamp)
                    continue;
                for (VertexId vertex_id_4: neighbors(data, vertex_id_3)) {
                    Distance delta = distance(data, vertex_id, vertex_id_3)
                        + distance(data, vertex_id_2, vertex_id_4)
                        - distance(data, vertex_id, vertex_id_2)
                        - distance(data, vertex_id_3, vertex_id_4);
                    if (best_delta > delta) {
                        best_delta = delta;
                        best = {vertex_id, vertex_id_2, vertex_id_3, vertex_id_4};
                    }
                }
            }
        }
        // Remove (0, 1) and (2, 3), add (0, 2) and (1, 3).
        remove_edge(data, best[0], best[1]);
        remove_edge(data, best[2], best[3]);
        add_link(data, best[0], best[2]);
        add_link(data, best[1], best[3]);
        compute_subtours(data);
    }
}

/*
 * Offspring
 */

template <typename Distances>
double entropy_term(
        const EaxData<Distances>& data,
        int64_t frequency)
{
    if (frequency <= 0)
        return 0;
    double x = (double)frequency / data.parameters.population_size;
    return -x * std::log(x);
}

/**
 * Build the offspring solution of the current intermediate solution
 * (after 'merge_subtours'): its edges removed from and added to pA, and
 * its evaluation.
 */
template <typename Distances>
void build_offspring(
        EaxData<Distances>& data,
        Offspring& offspring)
{
    offspring.removed_edges.clear();
    offspring.added_edges.clear();
    for (VertexId cut: data.cuts) {
        offspring.removed_edges.push_back(Edge(
                    data.order[cut],
                    data.order[(cut + 1 == data.number_of_vertices)? 0: cut + 1]));
    }
    for (VertexId vertex_id: data.linked_vertices)
        for (VertexId link: data.links[vertex_id])
            if (link != -1 && vertex_id < link)
                offspring.added_edges.push_back(Edge(vertex_id, link));
    // Remove the edges both removed and added.
    std::sort(offspring.removed_edges.begin(), offspring.removed_edges.end());
    std::sort(offspring.added_edges.begin(), offspring.added_edges.end());
    std::vector<Edge>& removed = offspring.removed_edges;
    std::vector<Edge>& added = offspring.added_edges;
    size_t i_removed = 0;
    size_t i_added = 0;
    size_t new_removed_size = 0;
    size_t new_added_size = 0;
    while (i_removed < removed.size() || i_added < added.size()) {
        if (i_added == added.size()
                || (i_removed < removed.size() && removed[i_removed] < added[i_added])) {
            removed[new_removed_size++] = removed[i_removed++];
        } else if (i_removed == removed.size() || added[i_added] < removed[i_removed]) {
            added[new_added_size++] = added[i_added++];
        } else {
            i_removed++;
            i_added++;
        }
    }
    removed.resize(new_removed_size);
    added.resize(new_added_size);

    Distance length_difference = 0;
    double entropy_difference = 0;
    for (const Edge& edge: removed) {
        length_difference -= distance(data, edge.vertex_id_1, edge.vertex_id_2);
        int64_t frequency = data.edge_frequencies.get(edge.vertex_id_1, edge.vertex_id_2);
        entropy_difference += entropy_term(data, frequency - 1) - entropy_term(data, frequency);
    }
    for (const Edge& edge: added) {
        length_difference += distance(data, edge.vertex_id_1, edge.vertex_id_2);
        int64_t frequency = data.edge_frequencies.get(edge.vertex_id_1, edge.vertex_id_2);
        entropy_difference += entropy_term(data, frequency + 1) - entropy_term(data, frequency);
    }
    offspring.length_difference = length_difference;
    // Entropy-preserving selection (the average length difference is
    // proportional to the length difference).
    double average_length_difference = (double)length_difference / data.parameters.population_size;
    if (length_difference >= 0) {
        offspring.evaluation = -average_length_difference;
    } else if (entropy_difference < 0) {
        offspring.evaluation = average_length_difference / entropy_difference;
    } else {
        offspring.evaluation = -average_length_difference / 1e-9;
    }
}

/** Replace a tour of the population by an offspring solution of it. */
template <typename Distances>
void apply_offspring(
        EaxData<Distances>& data,
        Tour& tour,
        const Offspring& offspring)
{
    for (const Edge& edge: offspring.removed_edges) {
        for (auto [a, b]: {std::pair{edge.vertex_id_1, edge.vertex_id_2}, std::pair{edge.vertex_id_2, edge.vertex_id_1}}) {
            auto& neighbors = tour.neighbors[a];
            if (neighbors[0] == b) {
                neighbors[0] = -1;
            } else {
                neighbors[1] = -1;
            }
        }
        data.edge_frequencies.add(edge.vertex_id_1, edge.vertex_id_2, -1);
    }
    for (const Edge& edge: offspring.added_edges) {
        for (auto [a, b]: {std::pair{edge.vertex_id_1, edge.vertex_id_2}, std::pair{edge.vertex_id_2, edge.vertex_id_1}}) {
            auto& neighbors = tour.neighbors[a];
            if (neighbors[0] == -1) {
                neighbors[0] = b;
            } else {
                neighbors[1] = b;
            }
        }
        data.edge_frequencies.add(edge.vertex_id_1, edge.vertex_id_2, 1);
    }
    tour.length += offspring.length_difference;
    tour.order_valid = false;
#ifndef NDEBUG
    // The tour must be a Hamiltonian cycle of the recorded length.
    Distance length = 0;
    VertexId previous_vertex_id = tour.neighbors[0][0];
    VertexId vertex_id = 0;
    for (VertexId pos = 0; pos < data.number_of_vertices; ++pos) {
        VertexId next_vertex_id = (tour.neighbors[vertex_id][0] != previous_vertex_id)?
            tour.neighbors[vertex_id][0]:
            tour.neighbors[vertex_id][1];
        length += distance(data, vertex_id, next_vertex_id);
        previous_vertex_id = vertex_id;
        vertex_id = next_vertex_id;
        assert(pos == data.number_of_vertices - 1 || vertex_id != 0);
    }
    assert(vertex_id == 0);
    assert(length == tour.length);
#endif
}

/*
 * Crossover
 */

/** Apply the AB-cycles of an E-set to pA, connect the subtours, and build the offspring. */
template <typename Distances>
void generate_offspring(
        EaxData<Distances>& data,
        const std::vector<VertexId>& e_set,
        Offspring& offspring)
{
    reset_intermediate(data);
    for (VertexId ab_cycle_id: e_set) {
        const std::vector<VertexId>& cycle = data.ab_cycles[ab_cycle_id];
        VertexId size = cycle.size();
        for (VertexId i = 0; i < size; i += 2) {
            VertexId cut = cut_position(data, cycle[i], cycle[i + 1]);
            data.is_cut[cut] = 1;
            data.cuts.push_back(cut);
            add_link(data, cycle[i + 1], cycle[(i + 2 == size)? 0: i + 2]);
        }
    }
    std::sort(data.cuts.begin(), data.cuts.end());
    merge_subtours(data);
    build_offspring(data, offspring);
}

/**
 * Generate offspring solutions of pA and pB, and replace pA by the best
 * one if it improves it (entropy-preserving selection), unless its
 * length is the one of pB (it's then very likely pB itself, and
 * replacing pA by it would decrease the diversity). Return 'true' iff pA
 * has been replaced.
 */
template <typename Distances>
bool crossover(
        EaxData<Distances>& data,
        Tour& tour_a,
        const Tour& tour_b,
        Strategy strategy)
{
    compute_order(data, tour_a);
    compute_ab_cycles(data, 
            tour_a,
            tour_b,
            (strategy == Strategy::single)? data.parameters.number_of_children: -1);
    if (data.number_of_ab_cycles == 0)
        return false;

    data.best_offspring.evaluation = 0;
    bool found = false;
    std::vector<VertexId> e_set;
    if (strategy == Strategy::single) {
        // Each offspring from a single AB-cycle, without repetition.
        std::vector<VertexId> ab_cycle_ids(data.number_of_ab_cycles);
        std::iota(ab_cycle_ids.begin(), ab_cycle_ids.end(), 0);
        std::shuffle(ab_cycle_ids.begin(), ab_cycle_ids.end(), data.generator);
        VertexId number_of_children = std::min(
                (VertexId)data.parameters.number_of_children,
                data.number_of_ab_cycles);
        for (VertexId child_id = 0; child_id < number_of_children; ++child_id) {
            e_set = {ab_cycle_ids[child_id]};
            generate_offspring(data, e_set, data.offspring);
            if (data.offspring.evaluation > data.best_offspring.evaluation
                    && tour_a.length + data.offspring.length_difference != tour_b.length) {
                std::swap(data.best_offspring, data.offspring);
                found = true;
            }
        }
    } else {
        for (VertexId child_id = 0; child_id < data.parameters.number_of_children; ++child_id) {
            if (!block2_e_set(data, child_id, e_set))
                break;
            generate_offspring(data, e_set, data.offspring);
            if (data.offspring.evaluation > data.best_offspring.evaluation
                    && tour_a.length + data.offspring.length_difference != tour_b.length) {
                std::swap(data.best_offspring, data.offspring);
                found = true;
            }
        }
    }
    reset_intermediate(data);
    if (found)
        apply_offspring(data, tour_a, data.best_offspring);
    return found;
}

/*
 * Block2 strategy
 */

/**
 * Compute, for the AB-cycles of the current pair of parents, the values
 * used by the block2 strategy: for each AB-cycle i, the number n_i of
 * C-vertices of the intermediate solution of the E-set {i}, and for each
 * pair of AB-cycles i and j, the number w_ij of their contact points (the
 * vertices visited by both).
 *
 * They're computed on "redefined" AB-cycles: each sequence of edges
 * common to both parents (an A-block whose edges are also edges of B) is
 * integrated into one of the two AB-cycles at its ends, so that it isn't
 * counted as an A-block between two C-vertices.
 */
template <typename Distances>
void prepare_block2(EaxData<Distances>& data)
{
    VertexId number_of_ab_cycles = data.number_of_ab_cycles;
    // Visits of the vertices by the AB-cycles (at most two).
    for (VertexId vertex_id: data.visited_vertices)
        data.visits[vertex_id].size = 0;
    data.visited_vertices.clear();
    auto add_visit = [&data](VertexId vertex_id, VertexId ab_cycle_id)
    {
        NeighborList& visits = data.visits[vertex_id];
        if (visits.size == 0)
            data.visited_vertices.push_back(vertex_id);
        visits.vertices[visits.size++] = ab_cycle_id;
    };
    for (VertexId ab_cycle_id = 0; ab_cycle_id < number_of_ab_cycles; ++ab_cycle_id)
        for (VertexId vertex_id: data.ab_cycles[ab_cycle_id])
            add_visit(vertex_id, ab_cycle_id);

    // Sequences of successive ineffective AB-cycles (of common edges),
    // along pA: each one is integrated into an effective AB-cycle
    // visiting one of its ends (each end is visited once by an effective
    // AB-cycle).
    auto is_common = [&data](VertexId pos)
    {
        return data.ineffective[pos] == data.ineffective_stamp;
    };
    // Start after a non-common edge.
    VertexId start_pos = 0;
    while (is_common(start_pos))
        start_pos++;
    start_pos = (start_pos + 1 == data.number_of_vertices)? 0: start_pos + 1;
    VertexId pos = start_pos;
    for (VertexId count = 0; count < data.number_of_vertices;) {
        if (!is_common(pos)) {
            pos = (pos + 1 == data.number_of_vertices)? 0: pos + 1;
            count++;
            continue;
        }
        // Sequence of common edges from the vertex at 'pos'.
        VertexId first_vertex_id = data.order[pos];
        while (count < data.number_of_vertices && is_common(pos)) {
            pos = (pos + 1 == data.number_of_vertices)? 0: pos + 1;
            count++;
        }
        VertexId last_vertex_id = data.order[pos];
        if (data.visits[first_vertex_id].size != 1 || data.visits[last_vertex_id].size != 1)
            continue;
        VertexId first_ab_cycle_id = data.visits[first_vertex_id].vertices[0];
        VertexId last_ab_cycle_id = data.visits[last_vertex_id].vertices[0];
        VertexId ab_cycle_id = (std::uniform_int_distribution<int>(0, 1)(data.generator) == 0)?
            first_ab_cycle_id:
            last_ab_cycle_id;
        add_visit(first_vertex_id, ab_cycle_id);
        add_visit(last_vertex_id, ab_cycle_id);
    }

    // n_i and w_ij.
    data.c_vertices.assign(number_of_ab_cycles, 0);
    data.contacts.resize(number_of_ab_cycles);
    for (VertexId ab_cycle_id = 0; ab_cycle_id < number_of_ab_cycles; ++ab_cycle_id)
        data.contacts[ab_cycle_id].clear();
    for (VertexId vertex_id: data.visited_vertices) {
        const NeighborList& visits = data.visits[vertex_id];
        if (visits.size == 1) {
            data.c_vertices[visits.vertices[0]]++;
        } else if (visits.vertices[0] != visits.vertices[1]) {
            VertexId ab_cycle_id_1 = visits.vertices[0];
            VertexId ab_cycle_id_2 = visits.vertices[1];
            data.c_vertices[ab_cycle_id_1]++;
            data.c_vertices[ab_cycle_id_2]++;
            add_contact(data, ab_cycle_id_1, ab_cycle_id_2);
            add_contact(data, ab_cycle_id_2, ab_cycle_id_1);
        }
    }

    // AB-cycles by decreasing size (the central AB-cycles).
    data.ab_cycles_by_size.resize(number_of_ab_cycles);
    std::iota(data.ab_cycles_by_size.begin(), data.ab_cycles_by_size.end(), 0);
    std::sort(
            data.ab_cycles_by_size.begin(),
            data.ab_cycles_by_size.end(),
            [&data](VertexId ab_cycle_id_1, VertexId ab_cycle_id_2)
            {
                return data.ab_cycles[ab_cycle_id_1].size() > data.ab_cycles[ab_cycle_id_2].size();
            });
    data.in_e_set.assign(number_of_ab_cycles, 0);
    data.contact_sums.assign(number_of_ab_cycles, 0);
    data.tabu.assign(number_of_ab_cycles, 0);
    data.candidate_positions.assign(number_of_ab_cycles, -1);
    data.candidates.clear();
    data.e_set_members.clear();
}

template <typename Distances>
void add_contact(
        EaxData<Distances>& data,
        VertexId ab_cycle_id_1,
        VertexId ab_cycle_id_2)
{
    for (auto& contact: data.contacts[ab_cycle_id_1]) {
        if (contact.first == ab_cycle_id_2) {
            contact.second++;
            return;
        }
    }
    data.contacts[ab_cycle_id_1].push_back({ab_cycle_id_2, 1});
}

/** Add an AB-cycle to the E-set of the tabu search (or remove it). */
template <typename Distances>
void block2_move(
        EaxData<Distances>& data,
        VertexId ab_cycle_id,
        bool add)
{
    int64_t sign = (add)? 1: -1;
    data.number_of_c_vertices += sign * (data.c_vertices[ab_cycle_id] - 2 * data.contact_sums[ab_cycle_id]);
    data.in_e_set[ab_cycle_id] = add;
    if (add) {
        data.e_set_members.push_back(ab_cycle_id);
        remove_candidate(data, ab_cycle_id);
    } else {
        data.e_set_members.erase(std::find(data.e_set_members.begin(), data.e_set_members.end(), ab_cycle_id));
        if (data.contact_sums[ab_cycle_id] > 0)
            add_candidate(data, ab_cycle_id);
    }
    for (const auto& contact: data.contacts[ab_cycle_id]) {
        VertexId other_ab_cycle_id = contact.first;
        data.contact_sums[other_ab_cycle_id] += sign * contact.second;
        if (!data.in_e_set[other_ab_cycle_id]) {
            if (data.contact_sums[other_ab_cycle_id] > 0) {
                add_candidate(data, other_ab_cycle_id);
            } else {
                remove_candidate(data, other_ab_cycle_id);
            }
        }
    }
}

template <typename Distances>
void add_candidate(
        EaxData<Distances>& data,
        VertexId ab_cycle_id)
{
    if (data.candidate_positions[ab_cycle_id] != -1)
        return;
    data.candidate_positions[ab_cycle_id] = data.candidates.size();
    data.candidates.push_back(ab_cycle_id);
}

template <typename Distances>
void remove_candidate(
        EaxData<Distances>& data,
        VertexId ab_cycle_id)
{
    VertexId pos = data.candidate_positions[ab_cycle_id];
    if (pos == -1)
        return;
    VertexId last_ab_cycle_id = data.candidates.back();
    data.candidates[pos] = last_ab_cycle_id;
    data.candidate_positions[last_ab_cycle_id] = pos;
    data.candidates.pop_back();
    data.candidate_positions[ab_cycle_id] = -1;
}

/**
 * E-set of the child 'child_id' with the block2 strategy: a tabu search
 * minimizing the number of C-vertices of the intermediate solution (an
 * estimate of its number of subtours), from an E-set made of the central
 * AB-cycle (the 'child_id'-th largest one) and of some of the smaller
 * AB-cycles in contact with it. Return 'false' if there's no AB-cycle
 * left to be central.
 */
template <typename Distances>
bool block2_e_set(
        EaxData<Distances>& data,
        VertexId child_id,
        std::vector<VertexId>& e_set)
{
    if (child_id == 0)
        prepare_block2(data);
    if (child_id >= (VertexId)data.ab_cycles_by_size.size())
        return false;
    // Reset the E-set.
    for (VertexId ab_cycle_id: std::vector<VertexId>(data.e_set_members))
        block2_move(data, ab_cycle_id, false);
    for (VertexId ab_cycle_id: std::vector<VertexId>(data.candidates))
        remove_candidate(data, ab_cycle_id);
    for (VertexId ab_cycle_id = 0; ab_cycle_id < data.number_of_ab_cycles; ++ab_cycle_id)
        data.tabu[ab_cycle_id] = 0;
    data.number_of_c_vertices = 0;

    VertexId central_ab_cycle_id = data.ab_cycles_by_size[child_id];
    VertexId central_size = data.ab_cycles[central_ab_cycle_id].size();
    block2_move(data, central_ab_cycle_id, true);
    for (const auto& contact: std::vector<std::pair<VertexId, int64_t>>(data.contacts[central_ab_cycle_id])) {
        VertexId ab_cycle_id = contact.first;
        if ((VertexId)data.ab_cycles[ab_cycle_id].size() < central_size
                && std::uniform_int_distribution<int>(0, 1)(data.generator) == 0) {
            block2_move(data, ab_cycle_id, true);
        }
    }

    int64_t best_number_of_c_vertices = data.number_of_c_vertices;
    std::vector<VertexId> best_e_set = data.e_set_members;
    int64_t iteration = 0;
    int64_t iterations_without_improvement = 0;
    while (iterations_without_improvement < data.parameters.tabu_maximum_number_of_iterations_without_improvement) {
        // Best move (random tie-breaking).
        VertexId best_ab_cycle_id = -1;
        bool best_add = false;
        int64_t best_delta = std::numeric_limits<int64_t>::max();
        int64_t number_of_ties = 0;
        auto consider = [&](VertexId ab_cycle_id, bool add)
        {
            int64_t delta = (add)?
                data.c_vertices[ab_cycle_id] - 2 * data.contact_sums[ab_cycle_id]:
                -data.c_vertices[ab_cycle_id] + 2 * data.contact_sums[ab_cycle_id];
            bool tabu = (iteration <= data.tabu[ab_cycle_id]);
            if (tabu && data.number_of_c_vertices + delta >= best_number_of_c_vertices)
                return;
            if (delta < best_delta) {
                best_delta = delta;
                best_ab_cycle_id = ab_cycle_id;
                best_add = add;
                number_of_ties = 1;
            } else if (delta == best_delta) {
                number_of_ties++;
                if (std::uniform_int_distribution<int64_t>(0, number_of_ties - 1)(data.generator) == 0) {
                    best_ab_cycle_id = ab_cycle_id;
                    best_add = add;
                }
            }
        };
        for (VertexId ab_cycle_id: data.candidates)
            consider(ab_cycle_id, true);
        for (VertexId ab_cycle_id: data.e_set_members)
            if (ab_cycle_id != central_ab_cycle_id)
                consider(ab_cycle_id, false);
        if (best_ab_cycle_id == -1)
            break;
        block2_move(data, best_ab_cycle_id, best_add);
        data.tabu[best_ab_cycle_id] = iteration
            + std::uniform_int_distribution<int64_t>(0, data.parameters.tabu_maximum_tenure)(data.generator);
        iteration++;
        if (data.number_of_c_vertices < best_number_of_c_vertices) {
            best_number_of_c_vertices = data.number_of_c_vertices;
            best_e_set = data.e_set_members;
            iterations_without_improvement = 0;
        } else {
            iterations_without_improvement++;
        }
    }
    e_set = best_e_set;
    return true;
}

/*
 * Generations
 */

template <typename Distances>
Distance best_length(const EaxData<Distances>& data)
{
    Distance best = data.population[0].length;
    for (const Tour& tour: data.population)
        best = std::min(best, tour.length);
    return best;
}

template <typename Distances>
void report_best(
        EaxData<Distances>& data,
        const std::string& comment)
{
    VertexId best_id = 0;
    for (VertexId i = 0; i < (VertexId)data.population.size(); ++i)
        if (data.population[i].length < data.population[best_id].length)
            best_id = i;
    if (!data.output.solution.feasible()
            || data.population[best_id].length < data.output.solution.distance()) {
        data.algorithm_formatter.update_solution(solution_from_tour(data, data.population[best_id]), comment);
    }
}

/**
 * Run a stage: generations until the best solution doesn't improve over
 * G / 10 generations, G being the number of generations at which it
 * first didn't improve over 1500 / N_ch generations.
 */
template <typename Distances>
void run_stage(
        EaxData<Distances>& data,
        Strategy strategy)
{
    VertexId population_size = data.population.size();
    std::vector<VertexId> permutation(population_size);
    std::iota(permutation.begin(), permutation.end(), 0);
    Distance best = best_length(data);
    int64_t generation = 0;
    int64_t stagnation = 0;
    int64_t g = -1;
    int64_t stagnation_limit = 1500 / data.parameters.number_of_children;
    std::string stage_name = (strategy == Strategy::single)? "stage 1": "stage 2";
    while (!data.parameters.timer.needs_to_end()) {
        std::shuffle(permutation.begin(), permutation.end(), data.generator);
        for (VertexId i = 0; i < population_size; ++i) {
            if (data.parameters.timer.needs_to_end())
                break;
            Tour& tour_a = data.population[permutation[i]];
            const Tour& tour_b = data.population[permutation[(i + 1) % population_size]];
            crossover(data, tour_a, tour_b, strategy);
        }
        generation++;
        Distance new_best = best_length(data);
        if (new_best < best) {
            best = new_best;
            stagnation = 0;
            report_best(data, stage_name + ", generation " + std::to_string(generation));
        } else {
            stagnation++;
        }
        if (g == -1 && stagnation >= stagnation_limit)
            g = generation;
        // Converged population.
        double average_length = 0;
        for (const Tour& tour: data.population)
            average_length += tour.length;
        average_length /= population_size;
        if (average_length - best < 0.001)
            break;
        if (g != -1 && stagnation >= g / 10)
            break;
    }
    report_best(data, stage_name + ", end");
}

template <typename Distances>
void run(EaxData<Distances>& data)
{
    initialize_population(data);
    if (data.number_of_vertices < 8)
        return;
    run_stage(data, Strategy::single);
    if (data.parameters.stage_2)
        run_stage(data, Strategy::block2);
}

template <typename Distances>
const Output eax_dispatch(
        const Distances& distances,
        const Instance& instance,
        std::mt19937_64& generator,
        const EaxParameters& parameters)
{
    Output output(instance);
    AlgorithmFormatter algorithm_formatter(parameters, output);
    algorithm_formatter.start("Genetic algorithm with EAX");
    algorithm_formatter.print_header();

    EaxData<Distances> data(distances, instance, generator, parameters, algorithm_formatter, output);
    run(data);

    algorithm_formatter.end();
    return output;
}

}

const Output travelingsalesmansolver::eax(
        const Instance& instance,
        std::mt19937_64& generator,
        const EaxParameters& parameters)
{
    return FUNCTION_WITH_DISTANCES(
            eax_dispatch,
            instance.distances(),
            instance,
            generator,
            parameters);
}
