#pragma once

#include "travelingsalesmansolver/solution.hpp"
#include "travelingsalesmansolver/algorithm_formatter.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <utility>
#include <vector>

namespace travelingsalesmansolver
{

struct LocalSearchParameters: Parameters
{
    /** Population size. */
    int population_size = 100;

    /** Number of children generated per generation. */
    int number_of_children = 30;

    /** Seed of the random number generator. */
    int seed = 0;
};

const Output local_search(
        const Instance& instance,
        const LocalSearchParameters& parameters = {});

template <typename Distances>
const Output local_search(
        const Distances& distances,
        const Instance& instance,
        const LocalSearchParameters& parameters = {});

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/**
 * Implementation of the Edge Assembly Crossover (EAX) genetic algorithm for
 * the TSP, adapted from Shujia Liu's C++ implementation
 * (https://github.com/Sugia/GA-for-TSP); see 'licenses/eax-ga/NOTICE.md' for
 * the full list of changes made to the original source.
 *
 * All mutable working state lives in 'LocalSearchData'; every operation is a
 * free function taking it (or a 'const&' to it for read-only queries) in
 * place of an implicit 'this', following the same 'LocalSearchData' +
 * free-function pattern as 'local_search_pfss_makespan.cpp'.
 */
namespace
{

/** A TSP tour, represented as a doubly-linked list over its vertices. */
struct Individual
{
    /** Create an individual with no vertices. */
    Individual() { }

    /** Create an individual for a given number of vertices. */
    Individual(VertexId number_of_vertices):
        neighbors(number_of_vertices) { }

    /** Checks if two individuals represent the same tour. */
    bool operator==(
            const Individual& individual) const
    {
        if (neighbors.size() != individual.neighbors.size())
            return false;
        if (length != individual.length)
            return false;

        VertexId current_vertex_id = 0;
        VertexId previous_vertex_id = -1;
        for (std::size_t i = 0; i < neighbors.size(); ++i) {
            VertexId next_vertex_id = (neighbors[current_vertex_id][0] == previous_vertex_id)? neighbors[current_vertex_id][1]: neighbors[current_vertex_id][0];
            if (individual.neighbors[current_vertex_id][0] != next_vertex_id && individual.neighbors[current_vertex_id][1] != next_vertex_id)
                return false;
            previous_vertex_id = current_vertex_id;
            current_vertex_id = next_vertex_id;
        }
        return true;
    }

    /** 'neighbors[vertex_id]' holds the two vertices adjacent to 'vertex_id'. */
    std::vector<std::array<VertexId, 2>> neighbors;

    /** Length of the tour. */
    Distance length = 0;
};

/** Seed the random number generator used throughout the algorithm. */
inline void seed_random(int seed)
{
    std::srand(seed);
}

/** Random integer in ['min', 'max']. */
inline int random_integer(int min, int max)
{
    return min + (std::rand() % (max - min + 1));
}

inline double random_double(double min, double max)
{
    return min + std::rand() % (int)(max - min);
}

/** Random number drawn from the normal distribution of mean 'mu' and standard deviation 'sigma'. */
inline double random_normal(double mu, double sigma)
{
    const double pi = 3.1415926;
    double u1;
    do {
        u1 = random_double(0.0, 1.0);
    } while (u1 == 0.0);
    double u2 = random_double(0.0, 1.0);
    double x = std::sqrt(-2.0 * std::log(u1)) * std::cos(2 * pi * u2);
    return mu + sigma * x;
}

/**
 * Fill 'array' with 'number_of_samples' distinct random values from
 * [0, number_of_elements[.
 */
inline void random_permutation(
        std::vector<int>& array,
        int number_of_elements,
        int number_of_samples)
{
    if (number_of_elements <= 0)
        return;
    std::vector<int> visited(number_of_elements, 0);
    for (int i = 0; i < number_of_samples; ++i) {
        int r = std::rand() % (number_of_elements - i);
        while (visited[r] == 1)
            r = (r + 1) % number_of_elements;
        array[i] = r;
        visited[r] = 1;
    }
}

/** Randomly shuffle the first 'number_of_elements' elements of 'array'. */
inline void random_shuffle(
        std::vector<int>& array,
        int number_of_elements)
{
    std::vector<int> shuffled_positions(number_of_elements);
    random_permutation(shuffled_positions, number_of_elements, number_of_elements);
    std::vector<int> original(array.begin(), array.begin() + number_of_elements);
    for (int i = 0; i < number_of_elements; ++i)
        array[i] = original[shuffled_positions[i]];
}

inline void selection_sort(
        std::vector<int>& values,
        int l,
        int r)
{
    for (int i = l; i < r; ++i) {
        int id = i;
        for (int j = i + 1; j <= r; ++j)
            if (values[j] < values[id])
                id = j;
        std::swap(values[i], values[id]);
    }
}

inline int quick_sort_partition(
        std::vector<int>& values,
        int l,
        int r)
{
    int id = l + std::rand() % (r - l + 1);
    std::swap(values[l], values[id]);
    id = l;
    for (int i = l + 1; i <= r; ++i)
        if (values[i] < values[l])
            std::swap(values[++id], values[i]);
    std::swap(values[l], values[id]);
    return id;
}

inline void quick_sort(
        std::vector<int>& values,
        int l,
        int r)
{
    if (l < r) {
        if (r - l < 20) { // utilizes selection sort for small batch of data
            selection_sort(values, l, r);
            return;
        }
        int mid = quick_sort_partition(values, l, r);
        quick_sort(values, l, mid - 1);
        quick_sort(values, mid + 1, r);
    }
}

/**
 * Fill the first 'number_of_indices' elements of 'sorted_indices' with the
 * indices, among the first 'number_of_values' elements of 'values', of the
 * 'number_of_indices' smallest ones, in increasing order of value.
 */
inline void sort_indices_ascending(
        const std::vector<int>& values,
        int number_of_values,
        std::vector<int>& sorted_indices,
        int number_of_indices)
{
    std::vector<int> checked(number_of_values, 0);
    for (int i = 0; i < number_of_indices; ++i) {
        int best_value = std::numeric_limits<int>::max();
        int best_index = 0;
        for (int j = 0; j < number_of_values; ++j) {
            if (values[j] < best_value && checked[j] == 0) {
                best_value = values[j];
                best_index = j;
            }
        }
        sorted_indices[i] = best_index;
        checked[best_index] = 1;
    }
}

/**
 * Fill the first 'number_of_indices' elements of 'sorted_indices' with the
 * indices, among the first 'number_of_values' elements of 'values', of the
 * 'number_of_indices' largest ones, in decreasing order of value.
 */
inline void sort_indices_descending(
        const std::vector<int>& values,
        int number_of_values,
        std::vector<int>& sorted_indices,
        int number_of_indices)
{
    std::vector<int> checked(number_of_values, 0);
    for (int i = 0; i < number_of_indices; ++i) {
        int best_value = std::numeric_limits<int>::min();
        int best_index = 0;
        for (int j = 0; j < number_of_values; ++j) {
            if (values[j] > best_value && checked[j] == 0) {
                best_value = values[j];
                best_index = j;
            }
        }
        sorted_indices[i] = best_index;
        checked[best_index] = 1;
    }
}

/** Sort the first 'number_of_values' elements of 'values' in increasing order. */
inline void sort_ascending(
        std::vector<int>& values,
        int number_of_values)
{
    quick_sort(values, 0, number_of_values - 1);
}
/**
 * All mutable working state for the EAX genetic algorithm, shared by every
 * free function below in place of the 'Evaluator'/'KOpt'/'Cross'/'Environment'
 * classes this was flattened from.
 *
 * A few fields that existed as same-named private members of two different
 * classes were deliberately renamed to avoid silently aliasing two distinct
 * concepts into one field once flattened (each original class kept its own
 * copy of these, so nothing here changes behavior, only the field name):
 * - 'KOpt::number_of_segments_' (current number of segments in the 2-opt
 *   tree representation) -> 'number_of_tree_segments', keeping
 *   'Cross::number_of_segments_' (path segments built while completing a
 *   child tour) as the plain 'number_of_segments'.
 * - 'Cross::max_stagnation_' ("Block2" eset local search stagnation limit,
 *   see 'search_eset()') -> 'eset_max_stagnation', keeping
 *   'Environment::max_stagnation_' (generation-level stage-1 -> stage-2
 *   stagnation limit) as the plain 'max_stagnation'.
 * Additionally, 'Cross::current_city_'/'previous_city_' (AB-cycle trace
 * state) were renamed to 'trace_current_city'/'trace_previous_city' to avoid
 * reading as a call to the free functions 'next_city()'/'previous_city()'
 * below (a collision only exposed once the class-private trailing
 * underscore was dropped).
 *
 * Fields that were genuinely the same value in two or more of the original
 * classes (all ultimately set from the same constructor argument) were
 * unified into one field: 'number_of_vertices' (was
 * 'Evaluator::number_of_vertices'/'KOpt::number_of_vertices_'/
 * 'Cross::number_of_vertices_'), 'population_size' (was
 * 'Cross::population_size_'/'Environment::population_size_'), and the
 * nearest-neighbor list size constant (was 'Evaluator::max_near_cities_'/
 * 'KOpt::max_near_cities_used_', both '50') into the file-scope constant
 * 'max_near_cities' below. 'Evaluator' itself no longer exists as a
 * sub-object: its fields ('distances', 'near_cities', 'number_of_vertices')
 * are now direct fields of 'LocalSearchData', accessed as 'data.distances'
 * instead of through 'evaluator_'.
 */
template <typename Distances>
struct LocalSearchData
{
    /** Constructor: binds 'distances'/'instance'/'parameters'/'output'/'algorithm_formatter' and sizes/computes everything that only depends on 'number_of_vertices'. */
    LocalSearchData(
            const Distances& distances,
            const Instance& instance,
            const LocalSearchParameters& parameters,
            const Output& output,
            AlgorithmFormatter& algorithm_formatter,
            VertexId number_of_vertices);

    ////////////////////////////////////////////////////////////////////////
    // Instance data (formerly 'Evaluator').
    ////////////////////////////////////////////////////////////////////////

    /** Distances between vertices. */
    const Distances& distances;

    /** Instance, used to build a 'Solution' every time a new best individual is found. */
    const Instance& instance;

    /** Algorithm parameters (used to check 'parameters.timer.needs_to_end()'). */
    const LocalSearchParameters& parameters;

    /** Used to cheaply check whether a candidate individual actually improves on the best solution found so far, before converting it to a 'Solution'. */
    const Output& output;

    /** Used to report every new best individual found, as soon as it is found. */
    AlgorithmFormatter& algorithm_formatter;

    /** Number of vertices. */
    VertexId number_of_vertices;

    /** 'near_cities[vertex_id][k]' is the k-th nearest vertex to 'vertex_id'. */
    std::vector<std::vector<VertexId>> near_cities;

    ////////////////////////////////////////////////////////////////////////
    // 2-opt local search state (formerly 'KOpt').
    ////////////////////////////////////////////////////////////////////////

    /** For each vertex, the vertices in whose nearest-neighbor list it appears. */
    std::vector<std::vector<VertexId>> inverse_near_list;

    /** Number of segments fixed at the start of the current 'optimize()' call. */
    int fixed_number_of_segments = 0;

    /** Current number of segments in the 2-opt tree representation. */
    int number_of_tree_segments = 0;

    /** Whether the move currently being applied reverses the tour direction. */
    int reversed = 0;

    /** Length of the tour, as tracked incrementally by the tree representation. */
    Distance tour_length = 0;

    /** 'tree_neighbors[vertex_id]' holds the two vertices adjacent to 'vertex_id' within its segment ('-1' at a segment end). */
    std::vector<std::array<VertexId, 2>> tree_neighbors;

    /** 'segment_neighbors[s]' holds the two segments adjacent to segment 's'. */
    std::vector<std::array<int, 2>> segment_neighbors;

    /** 'segment_endpoints[s]' holds the two endpoint vertices of segment 's'. */
    std::vector<std::array<VertexId, 2>> segment_endpoints;

    /** The four vertices ('t[1]'..'t[4]') involved in the move currently being considered. */
    std::array<VertexId, 5> t;

    /** 'city_segment[vertex_id]' is the segment 'vertex_id' belongs to. */
    std::vector<int> city_segment;

    /** Order of the vertices within their segment, used to tell direction/betweenness cheaply. */
    std::vector<int> city_order;

    /** Order of the segments along the tour. */
    std::vector<int> segment_order;

    /** Orientation of each segment (0 or 1). */
    std::vector<int> segment_orientation;

    /** Number of vertices in each segment. */
    std::vector<int> segment_size;

    /** Whether each vertex is still active (a candidate to start a new improving move from). */
    std::vector<int> active;

    /** Scratch array used to build 'individual_to_tree()'. */
    std::vector<VertexId> tree_array;

    /** Scratch array used by 'make_random_solution()'. */
    std::vector<VertexId> remaining;

    ////////////////////////////////////////////////////////////////////////
    // Crossover state (formerly 'Cross').
    ////////////////////////////////////////////////////////////////////////

    /** Population size. */
    int population_size;

    int random_pick;
    int start_new_trace;
    int cycle_complete;
    int traversal_type;
    int number_of_unbranched;
    int number_of_branching;

    /** Start vertex of the AB-cycle currently being traced by 'set_ab_cycle()'/'form_ab_cycle()'. */
    VertexId trace_start;

    /** Vertex currently being visited by 'set_ab_cycle()'/'form_ab_cycle()'. */
    VertexId trace_current_city;

    /** Vertex visited just before 'trace_current_city' by 'set_ab_cycle()'/'form_ab_cycle()'. */
    VertexId trace_previous_city;

    int start_appearance_count;
    int evaluation_type;
    int eset_strategy;
    int number_of_ab_cycles;
    int position_current;
    int max_number_of_ab_cycles;

    std::vector<VertexId> unbranched;
    std::vector<VertexId> branching;
    std::vector<int> unbranched_index;
    std::vector<int> branching_index;
    std::vector<int> first_visit_position;
    std::vector<VertexId> route;
    std::vector<int> permutation;
    std::vector<VertexId> cycle_buffer;

    std::vector<std::vector<VertexId>> near_data;
    std::vector<std::vector<int>> ab_cycle;

    // speeds up start
    int number_of_units;
    int number_of_segments;
    int number_of_segment_positions;
    int number_of_elements_in_center_unit;
    int number_of_segments_for_center;
    Distance modification_gain;
    int number_of_modified_edges;
    int number_of_best_modified_edges;
    int number_of_applied_cycles;
    int number_of_best_applied_cycles;

    std::vector<VertexId> order;
    std::vector<int> inverse_order;
    std::vector<int> segment_unit;
    std::vector<int> segment_position_list;
    std::vector<int> link_a_position;
    std::vector<int> position_segment;
    std::vector<int> number_of_elements_in_unit;
    std::vector<int> center_unit;
    std::vector<VertexId> list_of_center_unit;
    std::vector<int> segment_for_center;
    std::vector<Distance> gain_ab;
    std::vector<int> applied_cycle;
    std::vector<int> best_applied_cycle;

    std::vector<std::vector<int>> segment;
    std::vector<std::vector<int>> link_b_position;
    std::vector<std::vector<VertexId>> modified_edge;
    std::vector<std::vector<VertexId>> best_modified_edge;
    // speeds up end

    // block2
    int number_of_used_ab_cycles;
    int number_of_c_nodes;
    int number_of_e_edges;
    int t_max;

    /** "Block2" eset local search stagnation limit (see 'search_eset()'). */
    int eset_max_stagnation;

    int number_of_ab_cycles_in_eset;

    /**
     * NOTE: this is deliberately never written by 'set_parents()' -- see
     * 'set_parents()''s local 'distance_ab_local' and 'licenses/eax-ga/NOTICE.md'
     * for why this reproduces a pre-existing upstream bug byte-for-byte.
     */
    int distance_ab;

    int best_number_of_c_nodes;
    int best_number_of_e_edges;

    std::vector<int> number_of_elements_in_ab_cycle;
    std::vector<int> weight_sr;
    std::vector<int> weight_c;
    std::vector<int> used_ab_cycle;
    std::vector<int> moved_ab_cycle;
    std::vector<int> ab_cycle_in_eset;

    std::vector<std::vector<int>> in_effect_node;
    std::vector<std::vector<int>> weight_rr;

    /** Number of children generated by the last 'run_cross()' call. */
    int number_of_generated_children = 0;

    ////////////////////////////////////////////////////////////////////////
    // Genetic algorithm state (formerly 'Environment').
    ////////////////////////////////////////////////////////////////////////

    /** Number of children generated per generation. */
    int number_of_children;

    /** Current population. */
    std::vector<Individual> population;

    /** Best individual found so far. */
    Individual best_individual;

    /** Number of generations elapsed for the current population. */
    int current_number_of_generations = 0;

    /** Accumulated number of children generated so far. */
    long int accumulated_number_of_children = 0;

    /** Edge frequency across the population. */
    std::vector<std::vector<int>> edge_frequency;

    /** Average tour length in the population. */
    double average_value = 0.0;

    /** Tour length of the best individual in the population. */
    Distance best_value = 0;

    /** Index of the best individual in the population. */
    int best_index = 0;

    /** Random pairing of the population for the current generation's crossovers. */
    std::vector<int> index_for_mating;

    /** Number of consecutive generations without an improvement to 'best_individual'. */
    int stagnation_count = 0;

    /** EAX method and eset-selection strategy (see 'run_cross()'): '[0]' 1:Greedy 3:Distance 4:Entropy, '[1]' 1:Single-AB 2:Block2. */
    int flags[10];

    /** Current stage of the genetic algorithm (1: Single-AB, 2: Block2). */
    int stage = 1;

    /** 'stagnation_count' threshold at which the algorithm moves to the next stage. */
    int max_stagnation = 0;

    /** Number of generations elapsed when stage 1 completed. */
    int number_of_generations_stage_1 = 0;
};

/** Number of nearest neighbors stored per vertex in 'near_cities'/'inverse_near_list'. */
static constexpr int max_near_cities = 50;

/** The orientation opposite to 'orientation' (i.e. '1 - orientation'). */
inline int opposite_orientation(
        int orientation)
{
    return 1 - orientation;
}

/** Compute 'near_cities' from 'data.distances'. */
template <typename Distances>
void compute_near_cities(
        LocalSearchData<Distances>& data)
{
    std::vector<int> checked(data.number_of_vertices);
    for (VertexId vertex_id = 0; vertex_id < data.number_of_vertices; ++vertex_id) {
        std::fill(checked.begin(), checked.end(), 0);
        checked[vertex_id] = 1;
        data.near_cities[vertex_id][0] = vertex_id;
        for (int k = 1; k <= max_near_cities; ++k) {
            VertexId closest_vertex_id = -1;
            Distance min_distance = std::numeric_limits<Distance>::max();
            for (VertexId other_vertex_id = 0; other_vertex_id < data.number_of_vertices; ++other_vertex_id) {
                if (checked[other_vertex_id] == 0
                        && data.distances.distance(vertex_id, other_vertex_id) <= min_distance) {
                    closest_vertex_id = other_vertex_id;
                    min_distance = data.distances.distance(vertex_id, other_vertex_id);
                }
            }
            data.near_cities[vertex_id][k] = closest_vertex_id;
            checked[closest_vertex_id] = 1;
        }
    }
}

template <typename Distances>
LocalSearchData<Distances>::LocalSearchData(
        const Distances& distances,
        const Instance& instance,
        const LocalSearchParameters& parameters,
        const Output& output,
        AlgorithmFormatter& algorithm_formatter,
        VertexId number_of_vertices):
    distances(distances),
    instance(instance),
    parameters(parameters),
    output(output),
    algorithm_formatter(algorithm_formatter),
    number_of_vertices(number_of_vertices),
    near_cities(number_of_vertices, std::vector<VertexId>(max_near_cities + 1)),
    inverse_near_list(number_of_vertices),
    tree_neighbors(number_of_vertices),
    segment_neighbors(number_of_vertices),
    segment_endpoints(number_of_vertices),
    t{},
    city_segment(number_of_vertices),
    city_order(number_of_vertices),
    segment_order(number_of_vertices),
    segment_orientation(number_of_vertices),
    segment_size(number_of_vertices),
    active(number_of_vertices),
    tree_array(number_of_vertices + 2),
    remaining(number_of_vertices)
{
    compute_near_cities(*this);

    for (VertexId vertex_id = 0; vertex_id < number_of_vertices; ++vertex_id) {
        for (int k = 0; k < max_near_cities; ++k) {
            VertexId near_vertex_id = near_cities[vertex_id][k];
            inverse_near_list[near_vertex_id].push_back(vertex_id);
        }
    }
}

/** Size/fill every field of 'data' that depends on 'population_size'/'number_of_children'. */
template <typename Distances>
void init(
        LocalSearchData<Distances>& data,
        int population_size,
        int number_of_children)
{
    VertexId number_of_vertices = data.number_of_vertices;

    // Crossover state (formerly Cross's constructor).
    data.population_size = population_size;
    data.max_number_of_ab_cycles = 2000; // sets the maximum number of ab cycle

    data.near_data.clear();
    for (int i = 0; i < number_of_vertices; i++) {
        std::vector<VertexId> row(5);
        data.near_data.push_back(row);
    }

    data.ab_cycle.clear();
    for (int i = 0; i < data.max_number_of_ab_cycles; i++) {
        std::vector<int> row(2 * number_of_vertices + 4);
        data.ab_cycle.push_back(row);
    }

    data.unbranched.resize(number_of_vertices);
    data.branching.resize(number_of_vertices);
    data.unbranched_index.resize(number_of_vertices);
    data.branching_index.resize(number_of_vertices);
    data.first_visit_position.resize(number_of_vertices);
    data.route.resize(2 * number_of_vertices + 1);
    data.permutation.resize(data.max_number_of_ab_cycles);

    data.cycle_buffer.resize(2 * number_of_vertices + 4);

    // speeds up start
    data.order.resize(number_of_vertices);
    data.inverse_order.resize(number_of_vertices);

    // 'number_of_segment_positions'/'number_of_segments' are accumulated across every AB-cycle applied
    // within one crossover call (see 'change_sol()'/'make_unit()'), which in
    // the "Block2" eset mode (multiple AB-cycles per call) can exceed 'number_of_vertices'
    // even though each individual position is a valid index into a tour of
    // 'number_of_vertices' cities; size these generously (matching 'route'/'cycle_buffer' elsewhere
    // in this file) rather than assuming the count is bounded by 'number_of_vertices'.
    data.segment.clear();
    for (int i = 0; i < 2 * number_of_vertices; i++) {
        std::vector<int> row(2);
        data.segment.push_back(row);
    }

    data.segment_unit.resize(2 * number_of_vertices);
    data.segment_position_list.resize(2 * number_of_vertices);
    data.link_a_position.resize(number_of_vertices);

    data.link_b_position.clear();
    for (int i = 0; i < number_of_vertices; i++) {
        std::vector<int> row(2);
        data.link_b_position.push_back(row);
    }

    data.position_segment.resize(number_of_vertices);
    data.number_of_elements_in_unit.resize(number_of_vertices);

    data.center_unit.resize(number_of_vertices);
    for (int i = 0; i < number_of_vertices; i++) {
        data.center_unit[i] = 0;
    }

    data.list_of_center_unit.resize(number_of_vertices + 2);
    data.segment_for_center.resize(number_of_vertices);
    data.gain_ab.resize(number_of_vertices);

    data.modified_edge.clear();
    for (int i = 0; i < number_of_vertices; i++) {
        std::vector<VertexId> row(4);
        data.modified_edge.push_back(row);
    }

    data.best_modified_edge.clear();
    for (int i = 0; i < number_of_vertices; i++) {
        std::vector<VertexId> row(4);
        data.best_modified_edge.push_back(row);
    }

    data.applied_cycle.resize(number_of_vertices);
    data.best_applied_cycle.resize(number_of_vertices);
    // Speed Up End

    // block2
    data.number_of_elements_in_ab_cycle.resize(data.max_number_of_ab_cycles);

    data.in_effect_node.clear();
    for (int i = 0; i < number_of_vertices; i++) {
        std::vector<int> row(2);
        data.in_effect_node.push_back(row);
    }

    data.weight_rr.clear();
    for (int i = 0; i < data.max_number_of_ab_cycles; i++) {
        std::vector<int> row(data.max_number_of_ab_cycles);
        data.weight_rr.push_back(row);
    }

    data.weight_sr.resize(data.max_number_of_ab_cycles);
    data.weight_c.resize(data.max_number_of_ab_cycles);
    data.used_ab_cycle.resize(number_of_vertices);
    data.moved_ab_cycle.resize(number_of_vertices);
    data.ab_cycle_in_eset.resize(data.max_number_of_ab_cycles);

    // Genetic algorithm state (formerly Environment's constructor).
    data.number_of_children = number_of_children;
    data.population = std::vector<Individual>(population_size, Individual(number_of_vertices));
    data.best_individual = Individual(number_of_vertices);
    data.edge_frequency = std::vector<std::vector<int>>(number_of_vertices, std::vector<int>(number_of_vertices));
    data.index_for_mating.resize(population_size + 1);
}

/** Compute and store the length of an individual's tour. */
template <typename Distances>
void evaluate(
        const LocalSearchData<Distances>& data,
        Individual& individual)
{
    Distance d = 0;
    for (VertexId vertex_id = 0; vertex_id < data.number_of_vertices; ++vertex_id)
        d += data.distances.distance(vertex_id, individual.neighbors[vertex_id][0]) + data.distances.distance(vertex_id, individual.neighbors[vertex_id][1]);
    individual.length = d / 2;
}

/** Build a 'Solution' from the individual's tour. */
template <typename Distances>
Solution to_solution(
        const LocalSearchData<Distances>& data,
        const Individual& individual)
{
    // 'Solution' already starts at vertex 0, so it must not be added again
    // here.
    Solution solution(data.instance);
    VertexId current_vertex_id = 0;
    VertexId start_vertex_id = 0;
    VertexId previous_vertex_id = -1;
    for (int count = 0; count < data.number_of_vertices; ++count) {
        VertexId next_vertex_id = (individual.neighbors[current_vertex_id][0] == previous_vertex_id)?
            individual.neighbors[current_vertex_id][1]:
            individual.neighbors[current_vertex_id][0];
        previous_vertex_id = current_vertex_id;
        current_vertex_id = next_vertex_id;
        if (current_vertex_id == start_vertex_id)
            break;
        solution.add_vertex(data.distances, current_vertex_id);
    }
    return solution;
}

/** Build the tree representation of 'individual'. */
template <typename Distances>
void individual_to_tree(
        LocalSearchData<Distances>& data,
        const Individual& individual)
{
    data.tree_array[1] = 0;
    for (int i = 2; i <= data.number_of_vertices; ++i)
        data.tree_array[i] = individual.neighbors[data.tree_array[i - 1]][1];
    data.tree_array[0] = data.tree_array[data.number_of_vertices];
    data.tree_array[data.number_of_vertices + 1] = data.tree_array[1];

    int num = 1;
    data.number_of_tree_segments = 0;
    while (true) {
        int orientation = 1;
        int size = 0;
        data.segment_orientation[data.number_of_tree_segments] = orientation;
        data.segment_order[data.number_of_tree_segments] = data.number_of_tree_segments;

        data.tree_neighbors[data.tree_array[num]][0] = -1;
        data.tree_neighbors[data.tree_array[num]][1] = data.tree_array[num + 1];
        data.city_order[data.tree_array[num]] = size;
        data.city_segment[data.tree_array[num]] = data.number_of_tree_segments;
        data.segment_endpoints[data.number_of_tree_segments][opposite_orientation(orientation)] = data.tree_array[num];
        ++num;
        ++size;
        for (int i = 0; i < (int)std::sqrt(data.number_of_vertices * 1.0) - 1; ++i) {
            if (num == data.number_of_vertices)
                break;
            data.tree_neighbors[data.tree_array[num]][0] = data.tree_array[num - 1];
            data.tree_neighbors[data.tree_array[num]][1] = data.tree_array[num + 1];
            data.city_order[data.tree_array[num]] = size;
            data.city_segment[data.tree_array[num]] = data.number_of_tree_segments;
            ++num;
            ++size;
        }
        if (num == data.number_of_vertices - 1) {
            data.tree_neighbors[data.tree_array[num]][0] = data.tree_array[num - 1];
            data.tree_neighbors[data.tree_array[num]][1] = data.tree_array[num + 1];
            data.city_order[data.tree_array[num]] = size;
            data.city_segment[data.tree_array[num]] = data.number_of_tree_segments;
            ++num;
            ++size;
        }
        data.tree_neighbors[data.tree_array[num]][0] = data.tree_array[num - 1];
        data.tree_neighbors[data.tree_array[num]][1] = -1;
        data.city_order[data.tree_array[num]] = size;
        data.city_segment[data.tree_array[num]] = data.number_of_tree_segments;
        data.segment_endpoints[data.number_of_tree_segments][orientation] = data.tree_array[num];
        ++num;
        ++size;
        data.segment_size[data.number_of_tree_segments] = size;
        ++data.number_of_tree_segments;
        if (num == data.number_of_vertices + 1)
            break;
    }
    for (int s = 1; s < data.number_of_tree_segments - 1; ++s) {
        data.segment_neighbors[s][0] = s - 1;
        data.segment_neighbors[s][1] = s + 1;
    }
    data.segment_neighbors[0][0] = data.number_of_tree_segments - 1;
    data.segment_neighbors[0][1] = 1;
    data.segment_neighbors[data.number_of_tree_segments - 1][0] = data.number_of_tree_segments - 2;
    data.segment_neighbors[data.number_of_tree_segments - 1][1] = 0;
    data.tour_length = individual.length;
    data.fixed_number_of_segments = data.number_of_tree_segments;
}

/** Vertex immediately after 'vertex_id' in the tour. */
template <typename Distances>
VertexId next_city(
        const LocalSearchData<Distances>& data,
        VertexId vertex_id)
{
    int seg = data.city_segment[vertex_id];
    int orientation = data.segment_orientation[seg];
    VertexId next_vertex_id = data.tree_neighbors[vertex_id][orientation];
    if (next_vertex_id == -1) {
        seg = data.segment_neighbors[seg][orientation];
        orientation = opposite_orientation(data.segment_orientation[seg]);
        next_vertex_id = data.segment_endpoints[seg][orientation];
    }
    return next_vertex_id;
}

/** Vertex immediately before 'vertex_id' in the tour. */
template <typename Distances>
VertexId previous_city(
        const LocalSearchData<Distances>& data,
        VertexId vertex_id)
{
    int seg = data.city_segment[vertex_id];
    int orientation = data.segment_orientation[seg];
    VertexId previous_vertex_id = data.tree_neighbors[vertex_id][opposite_orientation(orientation)];
    if (previous_vertex_id == -1) {
        seg = data.segment_neighbors[seg][opposite_orientation(orientation)];
        orientation = data.segment_orientation[seg];
        previous_vertex_id = data.segment_endpoints[seg][orientation];
    }
    return previous_vertex_id;
}

/** Rebuild 'individual' from the tree representation. */
template <typename Distances>
void tree_to_individual(
        const LocalSearchData<Distances>& data,
        Individual& individual)
{
    for (VertexId vertex_id = 0; vertex_id < data.number_of_vertices; ++vertex_id) {
        individual.neighbors[vertex_id][0] = previous_city(data, vertex_id);
        individual.neighbors[vertex_id][1] = next_city(data, vertex_id);
    }
    evaluate(data, individual);
}

/** Merge segment 'segment_2' into segment 'segment_1'. */
template <typename Distances>
void merge_segments(
        LocalSearchData<Distances>& data,
        int segment_1,
        int segment_2)
{
    VertexId t_s = 0, t_e = 0;
    int direction = 0, ord = 0, increment = 0;

    if (data.segment_neighbors[segment_1][data.segment_orientation[segment_1]] == segment_2) {
        data.tree_neighbors[data.segment_endpoints[segment_1][data.segment_orientation[segment_1]]][data.segment_orientation[segment_1]] =
            data.segment_endpoints[segment_2][opposite_orientation(data.segment_orientation[segment_2])];
        data.tree_neighbors[data.segment_endpoints[segment_2][opposite_orientation(data.segment_orientation[segment_2])]][opposite_orientation(data.segment_orientation[segment_2])] =
            data.segment_endpoints[segment_1][data.segment_orientation[segment_1]];
        ord = data.city_order[data.segment_endpoints[segment_1][data.segment_orientation[segment_1]]];

        data.segment_endpoints[segment_1][data.segment_orientation[segment_1]] = data.segment_endpoints[segment_2][data.segment_orientation[segment_2]];
        data.segment_neighbors[segment_1][data.segment_orientation[segment_1]] = data.segment_neighbors[segment_2][data.segment_orientation[segment_2]];
        int seg = data.segment_neighbors[segment_2][data.segment_orientation[segment_2]];
        data.segment_neighbors[seg][opposite_orientation(data.segment_orientation[seg])] = segment_1;

        t_s = data.segment_endpoints[segment_2][opposite_orientation(data.segment_orientation[segment_2])];
        t_e = data.segment_endpoints[segment_2][data.segment_orientation[segment_2]];
        direction = data.segment_orientation[segment_2];

        increment = (data.segment_orientation[segment_1] == 1) ? 1 : -1;
    } else if (data.segment_neighbors[segment_1][opposite_orientation(data.segment_orientation[segment_1])] == segment_2) {
        data.tree_neighbors[data.segment_endpoints[segment_1][opposite_orientation(data.segment_orientation[segment_1])]][opposite_orientation(data.segment_orientation[segment_1])] =
            data.segment_endpoints[segment_2][data.segment_orientation[segment_2]];
        data.tree_neighbors[data.segment_endpoints[segment_2][data.segment_orientation[segment_2]]][data.segment_orientation[segment_2]] =
            data.segment_endpoints[segment_1][opposite_orientation(data.segment_orientation[segment_1])];
        ord = data.city_order[data.segment_endpoints[segment_1][opposite_orientation(data.segment_orientation[segment_1])]];

        data.segment_endpoints[segment_1][opposite_orientation(data.segment_orientation[segment_1])] = data.segment_endpoints[segment_2][opposite_orientation(data.segment_orientation[segment_2])];
        data.segment_neighbors[segment_1][opposite_orientation(data.segment_orientation[segment_1])] = data.segment_neighbors[segment_2][opposite_orientation(data.segment_orientation[segment_2])];
        int seg = data.segment_neighbors[segment_2][opposite_orientation(data.segment_orientation[segment_2])];
        data.segment_neighbors[seg][data.segment_orientation[seg]] = segment_1;

        t_s = data.segment_endpoints[segment_2][data.segment_orientation[segment_2]];
        t_e = data.segment_endpoints[segment_2][opposite_orientation(data.segment_orientation[segment_2])];
        direction = opposite_orientation(data.segment_orientation[segment_2]);

        increment = (data.segment_orientation[segment_1] == 1) ? -1 : 1;
    }
    VertexId curr = t_s;
    ord = ord + increment;
    while (true) {
        data.city_segment[curr] = segment_1;
        data.city_order[curr] = ord;

        VertexId next_vertex_id = data.tree_neighbors[curr][direction];
        if (data.segment_orientation[segment_1] != data.segment_orientation[segment_2])
            std::swap(data.tree_neighbors[curr][0], data.tree_neighbors[curr][1]);

        if (curr == t_e)
            break;
        curr = next_vertex_id;
        ord += increment;
    }
    data.segment_size[segment_1] += data.segment_size[segment_2];
    --data.number_of_tree_segments;
}

/** Apply the 2-opt move found by 'optimize()' to the tree representation. */
template <typename Distances>
void apply_move(
        LocalSearchData<Distances>& data)
{
    VertexId t1_s, t1_e, t2_s, t2_e;

    if (data.reversed == 0) {
        t1_s = data.t[1]; t1_e = data.t[3]; t2_s = data.t[4]; t2_e = data.t[2];
    } else {
        t1_s = data.t[2]; t1_e = data.t[4]; t2_s = data.t[3]; t2_e = data.t[1];
    }

    int seg_t1_s = data.city_segment[t1_s];
    int ordSeg_t1_s = data.segment_order[seg_t1_s];
    int orient_t1_s = data.segment_orientation[seg_t1_s];
    int seg_t1_e = data.city_segment[t1_e];
    int ordSeg_t1_e = data.segment_order[seg_t1_e];
    int orient_t1_e = data.segment_orientation[seg_t1_e];
    int seg_t2_s = data.city_segment[t2_s];
    int ordSeg_t2_s = data.segment_order[seg_t2_s];
    int orient_t2_s = data.segment_orientation[seg_t2_s];
    int seg_t2_e = data.city_segment[t2_e];
    int ordSeg_t2_e = data.segment_order[seg_t2_e];
    int orient_t2_e = data.segment_orientation[seg_t2_e];

    //////////////////// Type1 ////////////////////////
    if ((seg_t1_s == seg_t1_e) && (seg_t1_s == seg_t2_s) && (seg_t1_s == seg_t2_e)) {
        if ((data.segment_orientation[seg_t1_s] == 1 && (data.city_order[t1_s] > data.city_order[t1_e])) ||
            (data.segment_orientation[seg_t1_s] == 0 && (data.city_order[t1_s] < data.city_order[t1_e]))) {
            std::swap(t1_s, t2_s);
            std::swap(t1_e, t2_e);
            std::swap(seg_t1_s, seg_t2_s);
            std::swap(seg_t1_e, seg_t2_e);
            std::swap(ordSeg_t1_s, ordSeg_t2_s);
            std::swap(ordSeg_t1_e, ordSeg_t2_e);
            std::swap(orient_t1_s, orient_t2_s);
            std::swap(orient_t1_e, orient_t2_e);
        }
        VertexId curr = t1_s;
        int ord = data.city_order[t1_e];
        while (true) {
            std::swap(data.tree_neighbors[curr][0], data.tree_neighbors[curr][1]);
            data.city_order[curr] = ord;
            if (curr == t1_e)
                break;
            curr = data.tree_neighbors[curr][opposite_orientation(orient_t1_s)];
            if (orient_t1_s == 0)
                ++ord;
            else
                --ord;
        }

        data.tree_neighbors[t2_e][orient_t1_s] = t1_e;
        data.tree_neighbors[t2_s][opposite_orientation(orient_t1_s)] = t1_s;
        data.tree_neighbors[t1_s][orient_t1_s] = t2_s;
        data.tree_neighbors[t1_e][opposite_orientation(orient_t1_s)] = t2_e;

        return;
    }
    //////////////////// Type1 ///////////////////////

    int numOfSeg1 = (ordSeg_t1_e >= ordSeg_t1_s) ?
        ordSeg_t1_e - ordSeg_t1_s + 1 :
        ordSeg_t1_e - ordSeg_t1_s + 1 + data.number_of_tree_segments;
    int numOfSeg2 = (ordSeg_t2_e >= ordSeg_t2_s) ?
        ordSeg_t2_e - ordSeg_t2_s + 1 :
        ordSeg_t2_e - ordSeg_t2_s + 1 + data.number_of_tree_segments;

    if (numOfSeg1 > numOfSeg2) {
        std::swap(numOfSeg1, numOfSeg2);
        std::swap(t1_s, t2_s);
        std::swap(t1_e, t2_e);
        std::swap(seg_t1_s, seg_t2_s);
        std::swap(seg_t1_e, seg_t2_e);
        std::swap(ordSeg_t1_s, ordSeg_t2_s);
        std::swap(ordSeg_t1_e, ordSeg_t2_e);
        std::swap(orient_t1_s, orient_t2_s);
        std::swap(orient_t1_e, orient_t2_e);
    }
    int flag_t2e_t1s = (data.tree_neighbors[t2_e][orient_t2_e] == -1) ? 1 : 0;
    int flag_t2s_t1e = (data.tree_neighbors[t2_s][opposite_orientation(orient_t2_s)] == -1) ? 1 : 0;

    int length_t1s_seg = std::abs(data.city_order[t2_e] - data.city_order[data.segment_endpoints[seg_t2_e][orient_t2_e]]);
    int length_t1e_seg = std::abs(data.city_order[t2_s] - data.city_order[data.segment_endpoints[seg_t2_s][opposite_orientation(orient_t2_s)]]);

    ///////////////////// Type2 /////////////////
    if (seg_t1_s == seg_t1_e) {
        if (flag_t2e_t1s == 1 && flag_t2s_t1e == 1) {
            orient_t1_s = opposite_orientation(data.segment_orientation[seg_t1_s]);
            data.segment_orientation[seg_t1_s] = orient_t1_s;
            data.segment_endpoints[seg_t1_s][orient_t1_s] = t1_s;
            data.segment_endpoints[seg_t1_s][opposite_orientation(orient_t1_s)] = t1_e;
            data.segment_neighbors[seg_t1_s][orient_t1_s] = seg_t2_s;
            data.segment_neighbors[seg_t1_s][opposite_orientation(orient_t1_s)] = seg_t2_e;
            return;
        }
        if (flag_t2e_t1s == 0 && flag_t2s_t1e == 1) {
            VertexId curr = t1_e;
            int ord = data.city_order[t1_s];
            while (true) {
                std::swap(data.tree_neighbors[curr][0], data.tree_neighbors[curr][1]);
                data.city_order[curr] = ord;
                if (curr == t1_s)
                    break;
                curr = data.tree_neighbors[curr][orient_t2_e];
                if (orient_t2_e == 0)
                    --ord;
                else
                    ++ord;
            }
            data.tree_neighbors[t2_e][orient_t2_e] = t1_e;
            data.tree_neighbors[t1_s][orient_t2_e] = -1;
            data.tree_neighbors[t1_e][opposite_orientation(orient_t2_e)] = t2_e;
            data.segment_endpoints[seg_t2_e][orient_t2_e] = t1_s;
            return;
        }
        if (flag_t2e_t1s == 1 && flag_t2s_t1e == 0) {
            VertexId curr = t1_s;
            int ord = data.city_order[t1_e];
            while (true) {
                std::swap(data.tree_neighbors[curr][0], data.tree_neighbors[curr][1]);
                data.city_order[curr] = ord;
                if (curr == t1_e)
                    break;
                curr = data.tree_neighbors[curr][opposite_orientation(orient_t2_s)];
                if (orient_t2_s == 0)
                    ++ord;
                else
                    --ord;
            }
            data.tree_neighbors[t2_s][opposite_orientation(orient_t2_s)] = t1_s;
            data.tree_neighbors[t1_e][opposite_orientation(orient_t2_s)] = -1;
            data.tree_neighbors[t1_s][orient_t2_s] = t2_s;
            data.segment_endpoints[seg_t2_s][opposite_orientation(orient_t2_s)] = t1_e;
            return;
        }
    }

    ///////////////////// Type3 /////////////////

    if (flag_t2e_t1s == 1) {
        data.segment_neighbors[seg_t1_s][opposite_orientation(orient_t1_s)] = seg_t2_s;
    } else {
        seg_t1_s = data.number_of_tree_segments++;
        orient_t1_s = orient_t2_e;
        data.tree_neighbors[t1_s][opposite_orientation(orient_t1_s)] = -1;
        data.tree_neighbors[data.segment_endpoints[seg_t2_e][orient_t2_e]][orient_t1_s] = -1;
        data.segment_orientation[seg_t1_s] = orient_t1_s;
        data.segment_size[seg_t1_s] = length_t1s_seg;
        data.segment_endpoints[seg_t1_s][opposite_orientation(orient_t1_s)] = t1_s;
        data.segment_endpoints[seg_t1_s][orient_t1_s] = data.segment_endpoints[seg_t2_e][orient_t2_e];
        data.segment_neighbors[seg_t1_s][opposite_orientation(orient_t1_s)] = seg_t2_s;
        data.segment_neighbors[seg_t1_s][orient_t1_s] = data.segment_neighbors[seg_t2_e][orient_t2_e];
        int seg = data.segment_neighbors[seg_t2_e][orient_t2_e];
        data.segment_neighbors[seg][opposite_orientation(data.segment_orientation[seg])] = seg_t1_s;
    }

    if (flag_t2s_t1e == 1) {
        data.segment_neighbors[seg_t1_e][orient_t1_e] = seg_t2_e;
    } else {
        seg_t1_e = data.number_of_tree_segments++;
        orient_t1_e = orient_t2_s;
        data.tree_neighbors[t1_e][orient_t1_e] = -1;
        data.tree_neighbors[data.segment_endpoints[seg_t2_s][opposite_orientation(orient_t2_s)]][opposite_orientation(orient_t1_e)] = -1;
        data.segment_orientation[seg_t1_e] = orient_t1_e;
        data.segment_size[seg_t1_e] = length_t1e_seg;
        data.segment_endpoints[seg_t1_e][orient_t1_e] = t1_e;
        data.segment_endpoints[seg_t1_e][opposite_orientation(orient_t1_e)] = data.segment_endpoints[seg_t2_s][opposite_orientation(orient_t2_s)];
        data.segment_neighbors[seg_t1_e][orient_t1_e] = seg_t2_e;
        data.segment_neighbors[seg_t1_e][opposite_orientation(orient_t1_e)] = data.segment_neighbors[seg_t2_s][opposite_orientation(orient_t2_s)];
        int seg = data.segment_neighbors[seg_t2_s][opposite_orientation(orient_t2_s)];
        data.segment_neighbors[seg][data.segment_orientation[seg]] = seg_t1_e;
    }

    data.tree_neighbors[t2_e][orient_t2_e] = -1;
    data.segment_size[seg_t2_e] -= length_t1s_seg;
    data.segment_endpoints[seg_t2_e][orient_t2_e] = t2_e;
    data.segment_neighbors[seg_t2_e][orient_t2_e] = seg_t1_e;
    data.tree_neighbors[t2_s][opposite_orientation(orient_t2_s)] = -1;
    data.segment_size[seg_t2_s] -= length_t1e_seg;
    data.segment_endpoints[seg_t2_s][opposite_orientation(orient_t2_s)] = t2_s;
    data.segment_neighbors[seg_t2_s][opposite_orientation(orient_t2_s)] = seg_t1_s;

    {
        int seg = seg_t1_e;
        while (true) {
            data.segment_orientation[seg] = opposite_orientation(data.segment_orientation[seg]);
            if (seg == seg_t1_s)
                break;
            seg = data.segment_neighbors[seg][data.segment_orientation[seg]];
        }
    }

    if (data.segment_size[seg_t2_e] < length_t1s_seg) {
        int seg = data.segment_neighbors[seg_t2_e][opposite_orientation(data.segment_orientation[seg_t2_e])];
        data.segment_neighbors[seg][data.segment_orientation[seg]] = seg_t1_s;
        seg = data.segment_neighbors[seg_t2_e][data.segment_orientation[seg_t2_e]];
        data.segment_neighbors[seg][opposite_orientation(data.segment_orientation[seg])] = seg_t1_s;
        seg = data.segment_neighbors[seg_t1_s][opposite_orientation(data.segment_orientation[seg_t1_s])];
        data.segment_neighbors[seg][data.segment_orientation[seg]] = seg_t2_e;
        seg = data.segment_neighbors[seg_t1_s][data.segment_orientation[seg_t1_s]];
        data.segment_neighbors[seg][opposite_orientation(data.segment_orientation[seg])] = seg_t2_e;

        std::swap(data.segment_orientation[seg_t2_e], data.segment_orientation[seg_t1_s]);
        std::swap(data.segment_size[seg_t2_e], data.segment_size[seg_t1_s]);
        std::swap(data.segment_endpoints[seg_t2_e][0], data.segment_endpoints[seg_t1_s][0]);
        std::swap(data.segment_endpoints[seg_t2_e][1], data.segment_endpoints[seg_t1_s][1]);
        std::swap(data.segment_neighbors[seg_t2_e][0], data.segment_neighbors[seg_t1_s][0]);
        std::swap(data.segment_neighbors[seg_t2_e][1], data.segment_neighbors[seg_t1_s][1]);
        std::swap(seg_t2_e, seg_t1_s);
    }

    if (data.segment_size[seg_t2_s] < length_t1e_seg) {
        int seg = data.segment_neighbors[seg_t2_s][opposite_orientation(data.segment_orientation[seg_t2_s])];
        data.segment_neighbors[seg][data.segment_orientation[seg]] = seg_t1_e;
        seg = data.segment_neighbors[seg_t2_s][data.segment_orientation[seg_t2_s]];
        data.segment_neighbors[seg][opposite_orientation(data.segment_orientation[seg])] = seg_t1_e;
        seg = data.segment_neighbors[seg_t1_e][opposite_orientation(data.segment_orientation[seg_t1_e])];
        data.segment_neighbors[seg][data.segment_orientation[seg]] = seg_t2_s;
        seg = data.segment_neighbors[seg_t1_e][data.segment_orientation[seg_t1_e]];
        data.segment_neighbors[seg][opposite_orientation(data.segment_orientation[seg])] = seg_t2_s;

        std::swap(data.segment_orientation[seg_t2_s], data.segment_orientation[seg_t1_e]);
        std::swap(data.segment_size[seg_t2_s], data.segment_size[seg_t1_e]);
        std::swap(data.segment_endpoints[seg_t2_s][0], data.segment_endpoints[seg_t1_e][0]);
        std::swap(data.segment_endpoints[seg_t2_s][1], data.segment_endpoints[seg_t1_e][1]);
        std::swap(data.segment_neighbors[seg_t2_s][0], data.segment_neighbors[seg_t1_e][0]);
        std::swap(data.segment_neighbors[seg_t2_s][1], data.segment_neighbors[seg_t1_e][1]);
        std::swap(seg_t2_s, seg_t1_e);
    }

    while (data.number_of_tree_segments > data.fixed_number_of_segments) {
        if (data.segment_size[data.segment_neighbors[data.number_of_tree_segments - 1][0]] <
                data.segment_size[data.segment_neighbors[data.number_of_tree_segments - 1][1]])
            merge_segments(data, data.segment_neighbors[data.number_of_tree_segments - 1][0], data.number_of_tree_segments - 1);
        else
            merge_segments(data, data.segment_neighbors[data.number_of_tree_segments - 1][1], data.number_of_tree_segments - 1);
    }
    int ordSeg = 0;
    int seg = 0;
    while (true) {
        data.segment_order[seg] = ordSeg;
        ++ordSeg;
        seg = data.segment_neighbors[seg][data.segment_orientation[seg]];
        if (seg == 0)
            break;
    }
}

/** Repeatedly apply improving 2-opt moves until none remain. */
template <typename Distances>
void optimize(
        LocalSearchData<Distances>& data)
{
    std::fill(data.active.begin(), data.active.end(), 1);
BEGIN:
    {
        VertexId t1_start = random_integer(0, data.number_of_vertices - 1);
        data.t[1] = t1_start;
        while (true) {
            data.t[1] = next_city(data, data.t[1]);
            if (data.active[data.t[1]] == 0)
                goto RETURN;
            data.reversed = 0;
            data.t[2] = previous_city(data, data.t[1]);
            for (int num1 = 1; num1 < max_near_cities; ++num1) {
                data.t[4] = data.near_cities[data.t[1]][num1];
                data.t[3] = previous_city(data, data.t[4]);
                Distance dis1 = data.distances.distance(data.t[1], data.t[2]) - data.distances.distance(data.t[1], data.t[4]);
                if (dis1 > 0) {
                    Distance dis2 = dis1 + data.distances.distance(data.t[3], data.t[4]) - data.distances.distance(data.t[3], data.t[2]);
                    if (dis2 > 0) {
                        apply_move(data);
                        for (int a = 1; a <= 4; ++a)
                            for (VertexId near_vertex : data.inverse_near_list[data.t[a]])
                                data.active[near_vertex] = 1;
                        goto BEGIN;
                    }
                } else {
                    break;
                }
            }
            data.reversed = 1;
            data.t[2] = next_city(data, data.t[1]);
            for (int num1 = 1; num1 < max_near_cities; ++num1) {
                data.t[4] = data.near_cities[data.t[1]][num1];
                data.t[3] = next_city(data, data.t[4]);
                Distance dis1 = data.distances.distance(data.t[1], data.t[2]) - data.distances.distance(data.t[1], data.t[4]);
                if (dis1 > 0) {
                    Distance dis2 = dis1 + data.distances.distance(data.t[3], data.t[4]) - data.distances.distance(data.t[3], data.t[2]);
                    if (dis2 > 0) {
                        apply_move(data);
                        for (int a = 1; a <= 4; ++a)
                            for (VertexId near_vertex : data.inverse_near_list[data.t[a]])
                                data.active[near_vertex] = 1;
                        goto BEGIN;
                    }
                } else {
                    break;
                }
            }
            data.active[data.t[1]] = 0;
RETURN:
            if (data.t[1] == t1_start)
                break;
        }
    }
}

/** Run the local search on 'individual'. */
template <typename Distances>
void run_kopt(
        LocalSearchData<Distances>& data,
        Individual& individual)
{
    individual_to_tree(data, individual);
    optimize(data);
    tree_to_individual(data, individual);
}

/** Set 'individual' to a random tour and run the local search on it. */
template <typename Distances>
void make_random_solution(
        LocalSearchData<Distances>& data,
        Individual& individual)
{
    for (VertexId vertex_id = 0; vertex_id < data.number_of_vertices; ++vertex_id)
        data.remaining[vertex_id] = vertex_id;
    std::vector<VertexId> gene(data.number_of_vertices);
    for (int i = 0; i < data.number_of_vertices; ++i) {
        int r = random_integer(0, data.number_of_vertices - i - 1);
        gene[i] = data.remaining[r];
        data.remaining[r] = data.remaining[data.number_of_vertices - i - 1];
    }

    for (int j = 1; j < data.number_of_vertices - 1; ++j) {
        individual.neighbors[gene[j]][0] = gene[j - 1];
        individual.neighbors[gene[j]][1] = gene[j + 1];
    }
    individual.neighbors[gene[0]][0] = gene[data.number_of_vertices - 1];
    individual.neighbors[gene[0]][1] = gene[1];
    individual.neighbors[gene[data.number_of_vertices - 1]][0] = gene[data.number_of_vertices - 2];
    individual.neighbors[gene[data.number_of_vertices - 1]][1] = gene[0];

    evaluate(data, individual);
}

/** Record the AB-cycle currently being traced by 'set_ab_cycle()'. */
template <typename Distances>
void form_ab_cycle(
        LocalSearchData<Distances>& data)
{
    VertexId cycle_start;
    VertexId visiting_city;
    VertexId stock;
    int start_count;
    int edge_type;
    int cycle_length;
    Distance diff;

    edge_type = (data.position_current % 2 == 0)? 1: 2;
    cycle_start = data.route[data.position_current];
    cycle_length = 0;
    data.cycle_buffer[cycle_length] = cycle_start;

    start_count = 0;
    while (true) {
        ++cycle_length;
        --data.position_current;
        visiting_city = data.route[data.position_current];
        if (data.near_data[visiting_city][0] == 2) {
            data.unbranched[data.unbranched_index[visiting_city]] = data.unbranched[data.number_of_unbranched - 1];
            data.unbranched_index[data.unbranched[data.number_of_unbranched - 1]] = data.unbranched_index[visiting_city];
            --data.number_of_unbranched;
            data.branching[data.number_of_branching] = visiting_city;
            data.branching_index[visiting_city] = data.number_of_branching;
            ++data.number_of_branching;
        } else if (data.near_data[visiting_city][0] == 1) {
            data.branching[data.branching_index[visiting_city]] = data.branching[data.number_of_branching - 1];
            data.branching_index[data.branching[data.number_of_branching - 1]] = data.branching_index[visiting_city];
            --data.number_of_branching;
        }

        --data.near_data[visiting_city][0];
        if (visiting_city == cycle_start)
            ++start_count;
        if (start_count == data.start_appearance_count)
            break;
        data.cycle_buffer[cycle_length] = visiting_city;
    }

    if (cycle_length == 2)
        return;

    data.ab_cycle[data.number_of_ab_cycles][0] = cycle_length;

    if (edge_type == 2) {
        stock = data.cycle_buffer[0];
        for (int j = 0; j < cycle_length - 1; ++j)
            data.cycle_buffer[j] = data.cycle_buffer[j + 1];
        data.cycle_buffer[cycle_length - 1] = stock;
    }

    for (int j = 0; j < cycle_length; ++j)
        data.ab_cycle[data.number_of_ab_cycles][j + 2] = data.cycle_buffer[j];

    data.ab_cycle[data.number_of_ab_cycles][1] = data.cycle_buffer[cycle_length - 1];
    data.ab_cycle[data.number_of_ab_cycles][cycle_length + 2] = data.cycle_buffer[0];
    data.ab_cycle[data.number_of_ab_cycles][cycle_length + 3] = data.cycle_buffer[1];

    data.cycle_buffer[cycle_length] = data.cycle_buffer[0];
    data.cycle_buffer[cycle_length + 1] = data.cycle_buffer[1];
    diff = 0;
    for (int j = 0; j < cycle_length / 2; ++j)
        diff = diff + data.distances.distance(data.cycle_buffer[2 * j], data.cycle_buffer[1 + 2 * j]) - data.distances.distance(data.cycle_buffer[1 + 2 * j], data.cycle_buffer[2 + 2 * j]);

    data.gain_ab[data.number_of_ab_cycles] = diff;
    ++data.number_of_ab_cycles;
}

/** Find the AB-cycles of a given pair of parents. */
template <typename Distances>
void set_ab_cycle(
        LocalSearchData<Distances>& data,
        const Individual& parent1,
        const Individual& parent2,
        int number_of_kids)
{
    data.number_of_branching = 0;
    data.number_of_unbranched = 0;
    for (VertexId vertex_id = 0; vertex_id < data.number_of_vertices; ++vertex_id) {
        data.near_data[vertex_id][1] = parent1.neighbors[vertex_id][0];
        data.near_data[vertex_id][3] = parent1.neighbors[vertex_id][1];
        data.near_data[vertex_id][0] = 2;

        data.unbranched[data.number_of_unbranched] = vertex_id;
        data.number_of_unbranched++;

        data.near_data[vertex_id][2] = parent2.neighbors[vertex_id][0];
        data.near_data[vertex_id][4] = parent2.neighbors[vertex_id][1];
    }
    for (int i = 0; i < data.number_of_vertices; ++i) {
        data.first_visit_position[i] = -1;
        data.unbranched_index[data.unbranched[i]] = i;
    }
    data.number_of_ab_cycles = 0;
    data.start_new_trace = 1;
    while (data.number_of_unbranched != 0) {
        if (data.start_new_trace == 1) {
            data.position_current = 0;
            data.random_pick = rand() % data.number_of_unbranched;
            data.trace_start = data.unbranched[data.random_pick];
            data.first_visit_position[data.trace_start] = data.position_current;
            data.route[data.position_current] = data.trace_start;
            data.trace_current_city = data.trace_start;
            data.traversal_type = 2;
        } else if (data.start_new_trace == 0) {
            data.trace_current_city = data.route[data.position_current];
        }

        data.cycle_complete = 0;
        while (data.cycle_complete == 0) {
            data.position_current++;
            data.trace_previous_city = data.trace_current_city;
            switch (data.traversal_type) {
            case 1:
                data.trace_current_city = data.near_data[data.trace_previous_city][data.position_current % 2 + 1];
                break;
            case 2:
                data.random_pick = rand() % 2;
                data.trace_current_city = data.near_data[data.trace_previous_city][data.position_current % 2 + 1 + 2 * data.random_pick];
                if (data.random_pick == 0)
                    std::swap(data.near_data[data.trace_previous_city][data.position_current % 2 + 1], data.near_data[data.trace_previous_city][data.position_current % 2 + 3]);
                break;
            case 3:
                data.trace_current_city = data.near_data[data.trace_previous_city][data.position_current % 2 + 3];
            }
            data.route[data.position_current] = data.trace_current_city;
            if (data.near_data[data.trace_current_city][0] == 2) {
                if (data.trace_current_city == data.trace_start) {
                    if (data.first_visit_position[data.trace_start] == 0) {
                        if ((data.position_current - data.first_visit_position[data.trace_start]) % 2 == 0) {
                            if (data.near_data[data.trace_start][data.position_current % 2 + 1] == data.trace_previous_city)
                                std::swap(data.near_data[data.trace_current_city][data.position_current % 2 + 1], data.near_data[data.trace_current_city][data.position_current % 2 + 3]);

                            data.start_appearance_count = 1;
                            form_ab_cycle(data);
                            if (data.flags[1] == 1 && data.number_of_ab_cycles == number_of_kids)
                                goto RETURN;
                            if (data.number_of_ab_cycles == data.max_number_of_ab_cycles)
                                goto RETURN;

                            data.start_new_trace = 0;
                            data.cycle_complete = 1;
                            data.traversal_type = 1;
                        } else {
                            std::swap(data.near_data[data.trace_current_city][data.position_current % 2 + 1], data.near_data[data.trace_current_city][data.position_current % 2 + 3]);
                            data.traversal_type = 2;
                        }
                        data.first_visit_position[data.trace_start] = data.position_current;
                    } else {
                        data.start_appearance_count = 2;
                        form_ab_cycle(data);
                        if (data.flags[1] == 1 && data.number_of_ab_cycles == number_of_kids)
                            goto RETURN;
                        if (data.number_of_ab_cycles == data.max_number_of_ab_cycles)
                            goto RETURN;

                        data.start_new_trace = 1;
                        data.cycle_complete = 1;
                    }
                } else if (data.first_visit_position[data.trace_current_city] == -1) {
                    data.first_visit_position[data.trace_current_city] = data.position_current;
                    if (data.near_data[data.trace_current_city][data.position_current % 2 + 1] == data.trace_previous_city)
                        std::swap(data.near_data[data.trace_current_city][data.position_current % 2 + 1], data.near_data[data.trace_current_city][data.position_current % 2 + 3]);
                    data.traversal_type = 2;
                } else if (data.first_visit_position[data.trace_current_city] > 0) {
                    std::swap(data.near_data[data.trace_current_city][data.position_current % 2 + 1], data.near_data[data.trace_current_city][data.position_current % 2 + 3]);
                    if ((data.position_current - data.first_visit_position[data.trace_current_city]) % 2 == 0) {
                        data.start_appearance_count = 1;
                        form_ab_cycle(data);
                        if (data.flags[1] == 1 && data.number_of_ab_cycles == number_of_kids)
                            goto RETURN;
                        if (data.number_of_ab_cycles == data.max_number_of_ab_cycles)
                            goto RETURN;

                        data.start_new_trace = 0;
                        data.cycle_complete = 1;
                        data.traversal_type = 1;
                    } else {
                        std::swap(data.near_data[data.trace_current_city][(data.position_current + 1) % 2 + 1], data.near_data[data.trace_current_city][(data.position_current + 1) % 2 + 3]);
                        data.traversal_type = 3;
                    }
                }
            } else if (data.near_data[data.trace_current_city][0] == 1) {
                if (data.trace_current_city == data.trace_start) {
                    data.start_appearance_count = 1;
                    form_ab_cycle(data);
                    if (data.flags[1] == 1 && data.number_of_ab_cycles == number_of_kids)
                        goto RETURN;
                    if (data.number_of_ab_cycles == data.max_number_of_ab_cycles)
                        goto RETURN;
                    data.start_new_trace = 1;
                    data.cycle_complete = 1;
                } else {
                    data.traversal_type = 1;
                }
            }
        }
    }
    while (data.number_of_branching != 0) {
        data.position_current = 0;
        data.random_pick = rand() % data.number_of_branching;
        data.trace_start = data.branching[data.random_pick];
        data.route[data.position_current] = data.trace_start;
        data.trace_current_city = data.trace_start;

        data.cycle_complete = 0;
        while (data.cycle_complete == 0) {
            data.trace_previous_city = data.trace_current_city;
            data.position_current++;
            data.trace_current_city = data.near_data[data.trace_previous_city][data.position_current % 2 + 1];
            data.route[data.position_current] = data.trace_current_city;
            if (data.trace_current_city == data.trace_start) {
                data.start_appearance_count = 1;
                form_ab_cycle(data);
                if (data.flags[1] == 1 && data.number_of_ab_cycles == number_of_kids)
                    goto RETURN;
                if (data.number_of_ab_cycles == data.max_number_of_ab_cycles)
                    goto RETURN;

                data.cycle_complete = 1;
            }
        }
    }
RETURN:
    if (data.number_of_ab_cycles == data.max_number_of_ab_cycles) {
        printf("data.max_number_of_ab_cycles(%d) must be increased\n", data.max_number_of_ab_cycles);
        exit(1);
    }
}

/** Block2 eset selection: weight each AB-cycle by its interaction with the others. */
template <typename Distances>
void set_weight(
        LocalSearchData<Distances>& data,
        const Individual& parent1,
        const Individual& parent2)
{
    int cycle_length;
    VertexId red_vertex_id_1, red_vertex_id_2, current_vertex_id, next_vertex_id, previous_vertex_id;
    int ab_number;

    for (VertexId vertex_id = 0; vertex_id < data.number_of_vertices; ++vertex_id) {
        data.in_effect_node[vertex_id][0] = -1;
        data.in_effect_node[vertex_id][1] = -1;
    }

    // Step 1:
    for (int s = 0; s < data.number_of_ab_cycles; ++s) {
        cycle_length = data.ab_cycle[s][0];
        for (int j = 0; j < cycle_length / 2; ++j) {
            red_vertex_id_1 = data.ab_cycle[s][2 * j + 2]; // red edge
            red_vertex_id_2 = data.ab_cycle[s][2 * j + 3];

            if (data.in_effect_node[red_vertex_id_1][0] == -1) {
                data.in_effect_node[red_vertex_id_1][0] = s;
            } else if (data.in_effect_node[red_vertex_id_1][1] == -1) {
                data.in_effect_node[red_vertex_id_1][1] = s;
            }

            if (data.in_effect_node[red_vertex_id_2][0] == -1) {
                data.in_effect_node[red_vertex_id_2][0] = s;
            } else if (data.in_effect_node[red_vertex_id_2][1] == -1) {
                data.in_effect_node[red_vertex_id_2][1] = s;
            }
        }
    }

    // Step 2:
    for (VertexId vertex_id = 0; vertex_id < data.number_of_vertices; ++vertex_id) {
        if (data.in_effect_node[vertex_id][0] != -1 && data.in_effect_node[vertex_id][1] == -1) {
            ab_number = data.in_effect_node[vertex_id][0];
            current_vertex_id = vertex_id;

            if (parent1.neighbors[current_vertex_id][0] != parent2.neighbors[current_vertex_id][0] && parent1.neighbors[current_vertex_id][0] != parent2.neighbors[current_vertex_id][1]) {
                previous_vertex_id = parent1.neighbors[current_vertex_id][0];
            } else if (parent1.neighbors[current_vertex_id][1] != parent2.neighbors[current_vertex_id][0] && parent1.neighbors[current_vertex_id][1] != parent2.neighbors[current_vertex_id][1]) {
                previous_vertex_id = parent1.neighbors[current_vertex_id][1];
            }

            while (true) {
                data.in_effect_node[current_vertex_id][1] = ab_number;

                if (parent1.neighbors[current_vertex_id][0] != previous_vertex_id) {
                    next_vertex_id = parent1.neighbors[current_vertex_id][0];
                } else if (parent1.neighbors[current_vertex_id][1] != previous_vertex_id) {
                    next_vertex_id = parent1.neighbors[current_vertex_id][1];
                }

                if (data.in_effect_node[next_vertex_id][0] == -1) {
                    data.in_effect_node[next_vertex_id][0] = ab_number;
                } else if (data.in_effect_node[next_vertex_id][1] == -1) {
                    data.in_effect_node[next_vertex_id][1] = ab_number;
                }

                if (data.in_effect_node[next_vertex_id][1] != -1)
                    break;
                previous_vertex_id = current_vertex_id;
                current_vertex_id = next_vertex_id;
            }
        }
    }

    // Step 3:

    for (int s1 = 0; s1 < data.number_of_ab_cycles; ++s1) {
        data.weight_c[s1] = 0;
        for (int s2 = 0; s2 < data.number_of_ab_cycles; ++s2)
            data.weight_rr[s1][s2] = 0;
    }

    for (VertexId vertex_id = 0; vertex_id < data.number_of_vertices; ++vertex_id) {
        if (data.in_effect_node[vertex_id][0] != -1 && data.in_effect_node[vertex_id][1] != -1) {
            ++data.weight_rr[data.in_effect_node[vertex_id][0]][data.in_effect_node[vertex_id][1]];
            ++data.weight_rr[data.in_effect_node[vertex_id][1]][data.in_effect_node[vertex_id][0]];
        }
        if (data.in_effect_node[vertex_id][0] != data.in_effect_node[vertex_id][1]) {
            ++data.weight_c[data.in_effect_node[vertex_id][0]];
            ++data.weight_c[data.in_effect_node[vertex_id][1]];
        }
    }
    for (int s1 = 0; s1 < data.number_of_ab_cycles; ++s1)
        data.weight_rr[s1][s1] = 0;
}

/** Prepare AB-cycles for a given pair of parents. */
template <typename Distances>
void set_parents(
        LocalSearchData<Distances>& data,
        const Individual& parent1,
        const Individual& parent2,
        int number_of_kids)
{
    set_ab_cycle(data, parent1, parent2, number_of_kids);

    // NOTE: this local is deliberately never written back to 'data.distance_ab'
    // -- this reproduces a pre-existing upstream quirk byte-for-byte: the
    // field 'data.distance_ab', read later in 'run_cross()' (guarding
    // '2 * data.best_number_of_e_edges < data.distance_ab'), is never
    // actually updated by this function, so 'run_cross()' always sees a
    // stale value from a previous call (or 0, on the very first call). See
    // NOTICE.md.
    int distance_ab_local = 0;
    VertexId start_vertex_id = 0;
    VertexId current_vertex_id = -1;
    VertexId next_vertex_id = start_vertex_id;
    VertexId previous_vertex_id;
    for (int i = 0; i < data.number_of_vertices; ++i) {
        previous_vertex_id = current_vertex_id;
        current_vertex_id = next_vertex_id;
        if (parent1.neighbors[current_vertex_id][0] != previous_vertex_id) {
            next_vertex_id = parent1.neighbors[current_vertex_id][0];
        } else {
            next_vertex_id = parent1.neighbors[current_vertex_id][1];
        }
        if (parent2.neighbors[current_vertex_id][0] != next_vertex_id && parent2.neighbors[current_vertex_id][1] != next_vertex_id)
            ++distance_ab_local;
        data.order[i] = current_vertex_id;
        data.inverse_order[current_vertex_id] = i;
    }

    if (data.flags[1] == 2) {
        data.t_max = 10;
        data.eset_max_stagnation = 20; // 1:Greedy LS, 20:Tabu Search
        set_weight(data, parent1, parent2);
    }
}

/** Apply (type == 1) or roll back (type == 2) an AB-cycle on 'child'. */
template <typename Distances>
void change_sol(
        LocalSearchData<Distances>& data,
        Individual& child,
        int ab_number,
        int type)
{
    int j;
    int cycle_length;
    VertexId red_vertex_id_1, red_vertex_id_2, blue_vertex_id_1, blue_vertex_id_2;
    int red_position_1, red_position_2, blue_position_1, blue_position_2;

    cycle_length = data.ab_cycle[ab_number][0];
    data.cycle_buffer[0] = data.ab_cycle[ab_number][0];

    if (type == 2) {
        for (j = 0; j < cycle_length + 3; ++j)
            data.cycle_buffer[cycle_length + 3 - j] = data.ab_cycle[ab_number][j + 1];
    } else {
        for (j = 1; j <= cycle_length + 3; ++j)
            data.cycle_buffer[j] = data.ab_cycle[ab_number][j];
    }

    for (j = 0; j < cycle_length / 2; ++j) {
        red_vertex_id_1 = data.cycle_buffer[2 + 2 * j];
        red_vertex_id_2 = data.cycle_buffer[3 + 2 * j];
        blue_vertex_id_1 = data.cycle_buffer[1 + 2 * j];
        blue_vertex_id_2 = data.cycle_buffer[4 + 2 * j];

        if (child.neighbors[red_vertex_id_1][0] == red_vertex_id_2) {
            child.neighbors[red_vertex_id_1][0] = blue_vertex_id_1;
        } else {
            child.neighbors[red_vertex_id_1][1] = blue_vertex_id_1;
        }
        if (child.neighbors[red_vertex_id_2][0] == red_vertex_id_1) {
            child.neighbors[red_vertex_id_2][0] = blue_vertex_id_2;
        } else {
            child.neighbors[red_vertex_id_2][1] = blue_vertex_id_2;
        }

        red_position_1 = data.inverse_order[red_vertex_id_1];
        red_position_2 = data.inverse_order[red_vertex_id_2];
        blue_position_1 = data.inverse_order[blue_vertex_id_1];
        blue_position_2 = data.inverse_order[blue_vertex_id_2];

        if (red_position_1 == 0 && red_position_2 == data.number_of_vertices - 1) {
            data.segment_position_list[data.number_of_segment_positions++] = red_position_1;
        } else if (red_position_1 == data.number_of_vertices - 1 && red_position_2 == 0) {
            data.segment_position_list[data.number_of_segment_positions++] = red_position_2;
        } else if (red_position_1 < red_position_2) {
            data.segment_position_list[data.number_of_segment_positions++] = red_position_2;
        } else if (red_position_2 < red_position_1) {
            data.segment_position_list[data.number_of_segment_positions++] = red_position_1;
        }

        data.link_b_position[red_position_1][1] = data.link_b_position[red_position_1][0];
        data.link_b_position[red_position_2][1] = data.link_b_position[red_position_2][0];
        data.link_b_position[red_position_1][0] = blue_position_1;
        data.link_b_position[red_position_2][0] = blue_position_2;
    }
}

/** The 5th step of EAX: complete the tour from the remaining path segments. */
template <typename Distances>
void make_complete_sol(
        LocalSearchData<Distances>& data,
        Individual& child)
{
    int j1, j2;
    VertexId unit_start, previous_vertex_id, current_vertex_id, next_vertex_id;
    VertexId vertex_id_1, vertex_id_2, vertex_id_3, vertex_id_4;
    VertexId best_vertex_id_1, best_vertex_id_2, best_vertex_id_3, best_vertex_id_4;
    int min_unit_city;
    int center_unit_index, selected_unit_index;
    Distance diff, max_diff;
    int near_num, near_search_limit;

    data.modification_gain = 0;
    while (data.number_of_units != 1) {
        min_unit_city = data.number_of_vertices + 12345;
        for (int u = 0; u < data.number_of_units; ++u)
            if (data.number_of_elements_in_unit[u] < min_unit_city) {
                center_unit_index = u;
                min_unit_city = data.number_of_elements_in_unit[u];
            }

        unit_start = -1;
        data.number_of_segments_for_center = 0;
        for (int s = 0; s < data.number_of_segments; ++s)
            if (data.segment_unit[s] == center_unit_index) {
                int posi = data.segment[s][0];
                unit_start = data.order[posi];
                data.segment_for_center[data.number_of_segments_for_center++] = s;
            }
        current_vertex_id = -1;
        next_vertex_id = unit_start;
        data.number_of_elements_in_center_unit = 0;
        while (true) {
            previous_vertex_id = current_vertex_id;
            current_vertex_id = next_vertex_id;
            data.center_unit[current_vertex_id] = 1;
            data.list_of_center_unit[data.number_of_elements_in_center_unit] = current_vertex_id;
            ++data.number_of_elements_in_center_unit;
            if (child.neighbors[current_vertex_id][0] != previous_vertex_id) {
                next_vertex_id = child.neighbors[current_vertex_id][0];
            } else {
                next_vertex_id = child.neighbors[current_vertex_id][1];
            }
            if (next_vertex_id == unit_start)
                break;
        }
        data.list_of_center_unit[data.number_of_elements_in_center_unit] = data.list_of_center_unit[0];
        data.list_of_center_unit[data.number_of_elements_in_center_unit + 1] = data.list_of_center_unit[1];

        max_diff = std::numeric_limits<Distance>::min();
        best_vertex_id_3 = -1;
        best_vertex_id_4 = -1;
        near_search_limit = 10;   // N_near
        // near_search_limit <= max_near_cities

    RESTART:
        for (int s = 1; s <= data.number_of_elements_in_center_unit; ++s) {
            vertex_id_1 = data.list_of_center_unit[s];

            for (near_num = 1; near_num <= near_search_limit; ++near_num) {
                vertex_id_3 = data.near_cities[vertex_id_1][near_num];
                if (data.center_unit[vertex_id_3] == 0) {
                    for (j1 = 0; j1 < 2; ++j1) {
                        vertex_id_2 = data.list_of_center_unit[s - 1 + 2 * j1];
                        for (j2 = 0; j2 < 2; ++j2) {
                            vertex_id_4 = child.neighbors[vertex_id_3][j2];
                            diff = data.distances.distance(vertex_id_1, vertex_id_2) + data.distances.distance(vertex_id_3, vertex_id_4) - data.distances.distance(vertex_id_1, vertex_id_3) - data.distances.distance(vertex_id_2, vertex_id_4);
                            if (diff > max_diff) {
                                best_vertex_id_1 = vertex_id_1;
                                best_vertex_id_2 = vertex_id_2;
                                best_vertex_id_3 = vertex_id_3;
                                best_vertex_id_4 = vertex_id_4;
                                max_diff = diff;
                            }
                            diff = data.distances.distance(vertex_id_1, vertex_id_2) + data.distances.distance(vertex_id_4, vertex_id_3) -
                                data.distances.distance(vertex_id_1, vertex_id_4) - data.distances.distance(vertex_id_2, vertex_id_3);
                            if (diff > max_diff) {
                                best_vertex_id_1 = vertex_id_1;
                                best_vertex_id_2 = vertex_id_2;
                                best_vertex_id_3 = vertex_id_4;
                                best_vertex_id_4 = vertex_id_3;
                                max_diff = diff;
                            }
                        }
                    }
                }
            }
        }

        if (best_vertex_id_3 == -1 && near_search_limit == 10) {
            near_search_limit = 50;
            goto RESTART;
        } else if (best_vertex_id_3 == -1 && near_search_limit == 50) {
            int random_center_index = rand() % (data.number_of_elements_in_center_unit - 1);
            vertex_id_1 = data.list_of_center_unit[random_center_index];
            vertex_id_2 = data.list_of_center_unit[random_center_index + 1];
            for (VertexId vertex_id = 0; vertex_id < data.number_of_vertices; ++vertex_id) {
                if (data.center_unit[vertex_id] == 0) {
                    best_vertex_id_1 = vertex_id_1;
                    best_vertex_id_2 = vertex_id_2;
                    best_vertex_id_3 = vertex_id;
                    best_vertex_id_4 = child.neighbors[vertex_id][0];
                    break;
                }
            }
            max_diff = data.distances.distance(best_vertex_id_1, best_vertex_id_2) + data.distances.distance(best_vertex_id_3, best_vertex_id_4) - data.distances.distance(vertex_id_1, best_vertex_id_3) - data.distances.distance(vertex_id_2, best_vertex_id_4);
        }

        if (child.neighbors[best_vertex_id_1][0] == best_vertex_id_2) {
            child.neighbors[best_vertex_id_1][0] = best_vertex_id_3;
        } else {
            child.neighbors[best_vertex_id_1][1] = best_vertex_id_3;
        }
        if (child.neighbors[best_vertex_id_2][0] == best_vertex_id_1) {
            child.neighbors[best_vertex_id_2][0] = best_vertex_id_4;
        } else {
            child.neighbors[best_vertex_id_2][1] = best_vertex_id_4;
        }
        if (child.neighbors[best_vertex_id_3][0] == best_vertex_id_4) {
            child.neighbors[best_vertex_id_3][0] = best_vertex_id_1;
        } else {
            child.neighbors[best_vertex_id_3][1] = best_vertex_id_1;
        }
        if (child.neighbors[best_vertex_id_4][0] == best_vertex_id_3) {
            child.neighbors[best_vertex_id_4][0] = best_vertex_id_2;
        } else {
            child.neighbors[best_vertex_id_4][1] = best_vertex_id_2;
        }

        data.modified_edge[data.number_of_modified_edges][0] = best_vertex_id_1;
        data.modified_edge[data.number_of_modified_edges][1] = best_vertex_id_2;
        data.modified_edge[data.number_of_modified_edges][2] = best_vertex_id_3;
        data.modified_edge[data.number_of_modified_edges][3] = best_vertex_id_4;
        ++data.number_of_modified_edges;

        data.modification_gain += max_diff;

        int best_position_3 = data.inverse_order[best_vertex_id_3];
        selected_unit_index = -1;
        for (int s = 0; s < data.number_of_segments; ++s)
            if (data.segment[s][0] <= best_position_3 && best_position_3 <= data.segment[s][1]) {
                selected_unit_index = data.segment_unit[s];
                break;
            }

        for (int s = 0; s < data.number_of_segments; ++s)
            if (data.segment_unit[s] == selected_unit_index)
                data.segment_unit[s] = center_unit_index;

        data.number_of_elements_in_unit[center_unit_index] += data.number_of_elements_in_unit[selected_unit_index];

        for (int s = 0; s < data.number_of_segments; ++s)
            if (data.segment_unit[s] == data.number_of_units - 1)
                data.segment_unit[s] = selected_unit_index;

        data.number_of_elements_in_unit[selected_unit_index] = data.number_of_elements_in_unit[data.number_of_units - 1];
        --data.number_of_units;

        for (int s = 0; s < data.number_of_elements_in_center_unit; ++s) {
            VertexId vertex_id = data.list_of_center_unit[s];
            data.center_unit[vertex_id] = 0;
        }
    }
}

/** The 5-1th step of EAX: group path segments into units. */
template <typename Distances>
void make_unit(
        LocalSearchData<Distances>& data)
{
    int flag = 1;
    for (int s = 0; s < data.number_of_segment_positions; ++s) {
        if (data.segment_position_list[s] == 0) {
            flag = 0;
            break;
        }
    }
    if (flag == 1) {
        data.segment_position_list[data.number_of_segment_positions++] = 0;
        data.link_b_position[data.number_of_vertices - 1][1]  = data.link_b_position[data.number_of_vertices - 1][0];
        data.link_b_position[0][1] = data.link_b_position[0][0];
        data.link_b_position[data.number_of_vertices - 1][0] = 0;
        data.link_b_position[0][0] = data.number_of_vertices - 1;
    }

    sort_ascending(data.segment_position_list, data.number_of_segment_positions);
    data.number_of_segments = data.number_of_segment_positions;
    for (int s = 0; s < data.number_of_segments - 1; ++s) {
        data.segment[s][0] = data.segment_position_list[s];
        data.segment[s][1] = data.segment_position_list[s + 1] - 1;
    }

    data.segment[data.number_of_segments - 1][0] = data.segment_position_list[data.number_of_segments - 1];
    data.segment[data.number_of_segments - 1][1] = data.number_of_vertices - 1;

    for (int s = 0; s < data.number_of_segments; ++s) {
        data.link_a_position[data.segment[s][0]] = data.segment[s][1];
        data.link_a_position[data.segment[s][1]] = data.segment[s][0];
        data.position_segment[data.segment[s][0]] = s;
        data.position_segment[data.segment[s][1]] = s;
    }

    for (int s = 0; s < data.number_of_segments; ++s)
        data.segment_unit[s] = -1;
    data.number_of_units = 0;

    int start_position, position1, position2, next_position, previous_position;
    int segment_number;
    while (1) {
        flag = 0;
        for (int s = 0; s < data.number_of_segments; ++s) {
            if (data.segment_unit[s] == -1) {
                start_position = data.segment[s][0];
                previous_position = -1;
                position1 = start_position;
                flag = 1;
                break;
            }
        }
        if (flag == 0)
            break;

        while (1) {
            segment_number = data.position_segment[position1];
            data.segment_unit[segment_number] = data.number_of_units;

            position2 = data.link_a_position[position1];
            next_position = data.link_b_position[position2][0];
            if (position1 == position2)
                if (next_position == previous_position)
                    next_position = data.link_b_position[position2][1];

            if (next_position == start_position) {
                ++data.number_of_units;
                break;
            }

            previous_position = position2;
            position1 = next_position;
        }
    }

    for (int s = 0; s < data.number_of_units; ++s)
        data.number_of_elements_in_unit[s] = 0;

    int unit_number = -1;
    int merged_segment_count = -1;
    for (int s = 0; s < data.number_of_segments; ++s) {
        if (data.segment_unit[s] != unit_number) {
            ++merged_segment_count;
            data.segment[merged_segment_count][0] = data.segment[s][0];
            data.segment[merged_segment_count][1] = data.segment[s][1];
            unit_number = data.segment_unit[s];
            data.segment_unit[merged_segment_count] = unit_number;
            data.number_of_elements_in_unit[unit_number] += data.segment[s][1] - data.segment[s][0] + 1;
        } else {
            data.segment[merged_segment_count][1] = data.segment[s][1];
            data.number_of_elements_in_unit[unit_number] += data.segment[s][1] - data.segment[s][0] + 1;
        }
    }
    data.number_of_segments = merged_segment_count + 1;
}

/** Roll back 'child' to parent 1. */
template <typename Distances>
void back_to_pa1(
        LocalSearchData<Distances>& data,
        Individual& child)
{
    VertexId vertex_id_1, vertex_id_2, vertex_id_3, vertex_id_4;
    int jnum;

    for (int s = data.number_of_modified_edges - 1; s >= 0; --s) {
        vertex_id_1 = data.modified_edge[s][0];
        vertex_id_3 = data.modified_edge[s][1];
        vertex_id_2 = data.modified_edge[s][2];
        vertex_id_4 = data.modified_edge[s][3];

        if (child.neighbors[vertex_id_1][0] == vertex_id_2) {
            child.neighbors[vertex_id_1][0] = vertex_id_3;
        } else {
            child.neighbors[vertex_id_1][1] = vertex_id_3;
        }
        if (child.neighbors[vertex_id_4][0] == vertex_id_3) {
            child.neighbors[vertex_id_4][0] = vertex_id_2;
        } else {
            child.neighbors[vertex_id_4][1] = vertex_id_2;
        }
        if (child.neighbors[vertex_id_2][0] == vertex_id_1) {
            child.neighbors[vertex_id_2][0] = vertex_id_4;
        } else {
            child.neighbors[vertex_id_2][1] = vertex_id_4;
        }
        if (child.neighbors[vertex_id_3][0] == vertex_id_4) {
            child.neighbors[vertex_id_3][0] = vertex_id_1;
        } else {
            child.neighbors[vertex_id_3][1] = vertex_id_1;
        }
    }

    for (int s = 0; s < data.number_of_applied_cycles; ++s) {
        jnum = data.applied_cycle[s];
        change_sol(data, child, jnum, 2);
    }
}

/** Set 'child' to the best child found by the last 'run_cross()' call. */
template <typename Distances>
void go_to_best(
        LocalSearchData<Distances>& data,
        Individual& child)
{
    VertexId vertex_id_1, vertex_id_2, vertex_id_3, vertex_id_4;
    int jnum;

    for (int s = 0; s < data.number_of_best_applied_cycles; ++s) {
        jnum = data.best_applied_cycle[s];
        change_sol(data, child, jnum, 1);
    }

    for (int s = 0; s < data.number_of_best_modified_edges; ++s) {
        vertex_id_1 = data.best_modified_edge[s][0];
        vertex_id_2 = data.best_modified_edge[s][1];
        vertex_id_3 = data.best_modified_edge[s][2];
        vertex_id_4 = data.best_modified_edge[s][3];

        if (child.neighbors[vertex_id_1][0] == vertex_id_2) {
            child.neighbors[vertex_id_1][0] = vertex_id_3;
        } else {
            child.neighbors[vertex_id_1][1] = vertex_id_3;
        }
        if (child.neighbors[vertex_id_2][0] == vertex_id_1) {
            child.neighbors[vertex_id_2][0] = vertex_id_4;
        } else {
            child.neighbors[vertex_id_2][1] = vertex_id_4;
        }
        if (child.neighbors[vertex_id_3][0] == vertex_id_4) {
            child.neighbors[vertex_id_3][0] = vertex_id_1;
        } else {
            child.neighbors[vertex_id_3][1] = vertex_id_1;
        }
        if (child.neighbors[vertex_id_4][0] == vertex_id_3) {
            child.neighbors[vertex_id_4][0] = vertex_id_2;
        } else {
            child.neighbors[vertex_id_4][1] = vertex_id_2;
        }
    }
}

/** Update 'data.edge_frequency' for the best child found by the last 'run_cross()' call. */
template <typename Distances>
void increment_edge_freq(
        LocalSearchData<Distances>& data)
{
    int j, jnum, cycle_length;
    VertexId red_vertex_id_1, red_vertex_id_2, blue_vertex_id_1, blue_vertex_id_2;
    VertexId vertex_id_1, vertex_id_2, vertex_id_3, vertex_id_4;

    for (int s = 0; s < data.number_of_best_applied_cycles; ++s) {
        jnum = data.best_applied_cycle[s];

        cycle_length = data.ab_cycle[jnum][0];
        data.cycle_buffer[0] = data.ab_cycle[jnum][0];

        for (j = 1; j <= cycle_length + 3; ++j)
            data.cycle_buffer[j] = data.ab_cycle[jnum][j];

        for (j = 0; j < cycle_length / 2; ++j) {
            red_vertex_id_1 = data.cycle_buffer[2 + 2 * j];
            red_vertex_id_2 = data.cycle_buffer[3 + 2 * j];
            blue_vertex_id_1 = data.cycle_buffer[1 + 2 * j];
            blue_vertex_id_2 = data.cycle_buffer[4 + 2 * j];

            ++data.edge_frequency[red_vertex_id_1][blue_vertex_id_1];
            --data.edge_frequency[red_vertex_id_1][red_vertex_id_2];
            --data.edge_frequency[red_vertex_id_2][red_vertex_id_1];
            ++data.edge_frequency[red_vertex_id_2][blue_vertex_id_2];
        }
    }
    for (int s = 0; s < data.number_of_best_modified_edges; ++s) {
        vertex_id_1 = data.best_modified_edge[s][0];
        vertex_id_2 = data.best_modified_edge[s][1];
        vertex_id_3 = data.best_modified_edge[s][2];
        vertex_id_4 = data.best_modified_edge[s][3];

        --data.edge_frequency[vertex_id_1][vertex_id_2];
        --data.edge_frequency[vertex_id_3][vertex_id_4];
        ++data.edge_frequency[vertex_id_1][vertex_id_3];
        ++data.edge_frequency[vertex_id_2][vertex_id_4];
        --data.edge_frequency[vertex_id_2][vertex_id_1];
        --data.edge_frequency[vertex_id_4][vertex_id_3];
        ++data.edge_frequency[vertex_id_3][vertex_id_1];
        ++data.edge_frequency[vertex_id_4][vertex_id_2];
    }
}

/** Change in the average tour length implied by preserving 'data.edge_frequency'. */
template <typename Distances>
int calc_adaptive_loss(
        LocalSearchData<Distances>& data)
{
    int j, jnum, cycle_length;
    VertexId red_vertex_id_1, red_vertex_id_2, blue_vertex_id_1, blue_vertex_id_2;
    VertexId vertex_id_1, vertex_id_2, vertex_id_3, vertex_id_4;
    double loss;

    loss = 0;
    for (int s = 0; s < data.number_of_applied_cycles; ++s) {
        jnum = data.applied_cycle[s];

        cycle_length = data.ab_cycle[jnum][0];
        data.cycle_buffer[0] = data.ab_cycle[jnum][0];

        for (j = 1; j <= cycle_length + 3; ++j)
            data.cycle_buffer[j] = data.ab_cycle[jnum][j];

        for (j = 0; j < cycle_length / 2; ++j) {
            red_vertex_id_1 = data.cycle_buffer[2 + 2 * j];
            red_vertex_id_2 = data.cycle_buffer[3 + 2 * j];
            blue_vertex_id_1 = data.cycle_buffer[1 + 2 * j];
            blue_vertex_id_2 = data.cycle_buffer[4 + 2 * j];

            loss -= (data.edge_frequency[red_vertex_id_1][red_vertex_id_2] - 1);
            loss -= (data.edge_frequency[red_vertex_id_2][red_vertex_id_1] - 1);
            loss += data.edge_frequency[red_vertex_id_2][blue_vertex_id_2];
            loss += data.edge_frequency[blue_vertex_id_2][red_vertex_id_2];

            --data.edge_frequency[red_vertex_id_1][red_vertex_id_2];
            --data.edge_frequency[red_vertex_id_2][red_vertex_id_1];
            ++data.edge_frequency[red_vertex_id_2][blue_vertex_id_2];
            ++data.edge_frequency[blue_vertex_id_2][red_vertex_id_2];
        }
    }
    for (int s = 0; s < data.number_of_modified_edges; ++s) {
        vertex_id_1 = data.modified_edge[s][0];
        vertex_id_2 = data.modified_edge[s][1];
        vertex_id_3 = data.modified_edge[s][2];
        vertex_id_4 = data.modified_edge[s][3];

        loss -= (data.edge_frequency[vertex_id_1][vertex_id_2] - 1);
        loss -= (data.edge_frequency[vertex_id_2][vertex_id_1] - 1);
        loss -= (data.edge_frequency[vertex_id_3][vertex_id_4] - 1);
        loss -= (data.edge_frequency[vertex_id_4][vertex_id_3] - 1);

        loss += data.edge_frequency[vertex_id_1][vertex_id_3];
        loss += data.edge_frequency[vertex_id_3][vertex_id_1];
        loss += data.edge_frequency[vertex_id_2][vertex_id_4];
        loss += data.edge_frequency[vertex_id_4][vertex_id_2];

        --data.edge_frequency[vertex_id_1][vertex_id_2];
        --data.edge_frequency[vertex_id_2][vertex_id_1];
        --data.edge_frequency[vertex_id_3][vertex_id_4];
        --data.edge_frequency[vertex_id_4][vertex_id_3];

        ++data.edge_frequency[vertex_id_1][vertex_id_3];
        ++data.edge_frequency[vertex_id_3][vertex_id_1];
        ++data.edge_frequency[vertex_id_2][vertex_id_4];
        ++data.edge_frequency[vertex_id_4][vertex_id_2];
    }
    for (int s = 0; s < data.number_of_applied_cycles; ++s) {
        jnum = data.applied_cycle[s];
        cycle_length = data.ab_cycle[jnum][0];
        data.cycle_buffer[0] = data.ab_cycle[jnum][0];
        for (j = 1; j <= cycle_length + 3; ++j)
            data.cycle_buffer[j] = data.ab_cycle[jnum][j];

        for (j = 0; j < cycle_length / 2; ++j) {
            red_vertex_id_1 = data.cycle_buffer[2 + 2 * j];
            red_vertex_id_2 = data.cycle_buffer[3 + 2 * j];
            blue_vertex_id_1 = data.cycle_buffer[1 + 2 * j];
            blue_vertex_id_2 = data.cycle_buffer[4 + 2 * j];

            ++data.edge_frequency[red_vertex_id_1][red_vertex_id_2];
            ++data.edge_frequency[red_vertex_id_2][red_vertex_id_1];
            --data.edge_frequency[red_vertex_id_2][blue_vertex_id_2];
            --data.edge_frequency[blue_vertex_id_2][red_vertex_id_2];
        }
    }
    for (int s = 0; s < data.number_of_modified_edges; ++s) {
        vertex_id_1 = data.modified_edge[s][0];
        vertex_id_2 = data.modified_edge[s][1];
        vertex_id_3 = data.modified_edge[s][2];
        vertex_id_4 = data.modified_edge[s][3];

        ++data.edge_frequency[vertex_id_1][vertex_id_2];
        ++data.edge_frequency[vertex_id_2][vertex_id_1];
        ++data.edge_frequency[vertex_id_3][vertex_id_4];
        ++data.edge_frequency[vertex_id_4][vertex_id_3];

        --data.edge_frequency[vertex_id_1][vertex_id_3];
        --data.edge_frequency[vertex_id_3][vertex_id_1];
        --data.edge_frequency[vertex_id_2][vertex_id_4];
        --data.edge_frequency[vertex_id_4][vertex_id_2];
    }
    return int(loss / 2);
}

/** Change in edge-frequency entropy implied by 'data.edge_frequency'. */
template <typename Distances>
double calc_entropy_loss(
        LocalSearchData<Distances>& data)
{
    int j, jnum, cycle_length;
    VertexId red_vertex_id_1, red_vertex_id_2, blue_vertex_id_1, blue_vertex_id_2;
    VertexId vertex_id_1, vertex_id_2, vertex_id_3, vertex_id_4;
    double loss;
    double h1, h2;

    loss = 0;  // AB-cycle
    for (int s = 0; s < data.number_of_applied_cycles; ++s) {
        jnum = data.applied_cycle[s];
        cycle_length = data.ab_cycle[jnum][0];
        data.cycle_buffer[0] = data.ab_cycle[jnum][0];

        for (j = 1; j <= cycle_length + 3; ++j)
            data.cycle_buffer[j] = data.ab_cycle[jnum][j];

        for (j = 0; j < cycle_length / 2; ++j) {
            red_vertex_id_1 = data.cycle_buffer[2 + 2 * j];
            red_vertex_id_2 = data.cycle_buffer[3 + 2 * j];
            blue_vertex_id_1 = data.cycle_buffer[1 + 2 * j];
            blue_vertex_id_2 = data.cycle_buffer[4 + 2 * j];

            h1 = (double)(data.edge_frequency[red_vertex_id_1][red_vertex_id_2] - 1) / (double)data.population_size;
            h2 = (double)(data.edge_frequency[red_vertex_id_1][red_vertex_id_2]) / (double)data.population_size;
            if (data.edge_frequency[red_vertex_id_1][red_vertex_id_2] - 1 != 0)
                loss -= h1 * log(h1);
            loss += h2 * log(h2);
            --data.edge_frequency[red_vertex_id_1][red_vertex_id_2];
            --data.edge_frequency[red_vertex_id_2][red_vertex_id_1];

            h1 = (double)(data.edge_frequency[red_vertex_id_2][blue_vertex_id_2] + 1) / (double)data.population_size;
            h2 = (double)(data.edge_frequency[red_vertex_id_2][blue_vertex_id_2]) / (double)data.population_size;
            loss -= h1 * log(h1);
            if (data.edge_frequency[red_vertex_id_2][blue_vertex_id_2] != 0)
                loss += h2 * log(h2);
            ++data.edge_frequency[red_vertex_id_2][blue_vertex_id_2];
            ++data.edge_frequency[blue_vertex_id_2][red_vertex_id_2];
        }
    }

    for (int s = 0; s < data.number_of_modified_edges; ++s) {
        vertex_id_1 = data.modified_edge[s][0];
        vertex_id_2 = data.modified_edge[s][1];
        vertex_id_3 = data.modified_edge[s][2];
        vertex_id_4 = data.modified_edge[s][3];

        h1 = (double)(data.edge_frequency[vertex_id_1][vertex_id_2] - 1) / (double)data.population_size;
        h2 = (double)(data.edge_frequency[vertex_id_1][vertex_id_2]) / (double)data.population_size;
        if (data.edge_frequency[vertex_id_1][vertex_id_2] - 1 != 0)
            loss -= h1 * log(h1);
        loss += h2 * log(h2);
        --data.edge_frequency[vertex_id_1][vertex_id_2];
        --data.edge_frequency[vertex_id_2][vertex_id_1];

        h1 = (double)(data.edge_frequency[vertex_id_3][vertex_id_4] - 1) / (double)data.population_size;
        h2 = (double)(data.edge_frequency[vertex_id_3][vertex_id_4]) / (double)data.population_size;
        if (data.edge_frequency[vertex_id_3][vertex_id_4] - 1 != 0)
            loss -= h1 * log(h1);
        loss += h2 * log(h2);
        --data.edge_frequency[vertex_id_3][vertex_id_4];
        --data.edge_frequency[vertex_id_4][vertex_id_3];

        h1 = (double)(data.edge_frequency[vertex_id_1][vertex_id_3] + 1) / (double)data.population_size;
        h2 = (double)(data.edge_frequency[vertex_id_1][vertex_id_3]) / (double)data.population_size;
        loss -= h1 * log(h1);
        if (data.edge_frequency[vertex_id_1][vertex_id_3] != 0)
            loss += h2 * log(h2);
        ++data.edge_frequency[vertex_id_1][vertex_id_3];
        ++data.edge_frequency[vertex_id_3][vertex_id_1];

        h1 = (double)(data.edge_frequency[vertex_id_2][vertex_id_4] + 1) / (double)data.population_size;
        h2 = (double)(data.edge_frequency[vertex_id_2][vertex_id_4]) / (double)data.population_size;
        loss -= h1 * log(h1);
        if (data.edge_frequency[vertex_id_2][vertex_id_4] != 0)
            loss += h2 * log(h2);
        ++data.edge_frequency[vertex_id_2][vertex_id_4];
        ++data.edge_frequency[vertex_id_4][vertex_id_2];
    }
    loss = -loss;

    // restores data.edge_frequency
    for (int s = 0; s < data.number_of_applied_cycles; ++s) {
        jnum = data.applied_cycle[s];

        cycle_length = data.ab_cycle[jnum][0];
        data.cycle_buffer[0] = data.ab_cycle[jnum][0];

        for (j = 1; j <= cycle_length + 3; ++j)
            data.cycle_buffer[j] = data.ab_cycle[jnum][j];

        for (j = 0; j < cycle_length / 2; ++j) {
            red_vertex_id_1 = data.cycle_buffer[2 + 2 * j];
            red_vertex_id_2 = data.cycle_buffer[3 + 2 * j];
            blue_vertex_id_1 = data.cycle_buffer[1 + 2 * j];
            blue_vertex_id_2 = data.cycle_buffer[4 + 2 * j];

            ++data.edge_frequency[red_vertex_id_1][red_vertex_id_2];
            ++data.edge_frequency[red_vertex_id_2][red_vertex_id_1];
            --data.edge_frequency[red_vertex_id_2][blue_vertex_id_2];
            --data.edge_frequency[blue_vertex_id_2][red_vertex_id_2];
        }
    }
    for (int s = 0; s < data.number_of_modified_edges; ++s) {
        vertex_id_1 = data.modified_edge[s][0];
        vertex_id_2 = data.modified_edge[s][1];
        vertex_id_3 = data.modified_edge[s][2];
        vertex_id_4 = data.modified_edge[s][3];

        ++data.edge_frequency[vertex_id_1][vertex_id_2];
        ++data.edge_frequency[vertex_id_2][vertex_id_1];
        ++data.edge_frequency[vertex_id_3][vertex_id_4];
        ++data.edge_frequency[vertex_id_4][vertex_id_3];

        --data.edge_frequency[vertex_id_1][vertex_id_3];
        --data.edge_frequency[vertex_id_3][vertex_id_1];
        --data.edge_frequency[vertex_id_2][vertex_id_4];
        --data.edge_frequency[vertex_id_4][vertex_id_2];
    }
    return loss;
}

/** Add AB-cycle 'num' to the eset being built by 'search_eset()'. */
template <typename Distances>
void add_ab(
        LocalSearchData<Distances>& data,
        int num)
{
    data.number_of_c_nodes += data.weight_c[num] - 2 * data.weight_sr[num];
    data.number_of_e_edges += data.ab_cycle[num][0] / 2;

    data.used_ab_cycle[num] = 1;
    ++data.number_of_used_ab_cycles;
    for (int s1 = 0; s1 < data.number_of_ab_cycles; ++s1)
        data.weight_sr[s1] += data.weight_rr[s1][num];
}

/** Remove AB-cycle 'num' from the eset being built by 'search_eset()'. */
template <typename Distances>
void delete_ab(
        LocalSearchData<Distances>& data,
        int num)
{
    data.number_of_c_nodes -= data.weight_c[num] - 2 * data.weight_sr[num];
    data.number_of_e_edges -= data.ab_cycle[num][0] / 2;

    data.used_ab_cycle[num] = 0;
    --data.number_of_used_ab_cycles;
    for (int s1 = 0; s1 < data.number_of_ab_cycles; ++s1)
        data.weight_sr[s1] -= data.weight_rr[s1][num];
}

/** Block2 eset selection: local search over which AB-cycles to include around 'num'. */
template <typename Distances>
void search_eset(
        LocalSearchData<Distances>& data,
        int center_ab)
{
    int iteration, stagnation;
    int delta_weight, min_delta_weight_non_tabu;
    int improving_change, non_tabu_change;
    int selected_ab_cycle, selected_ab_cycle_non_tabu;
    int jnum;

    data.number_of_c_nodes = 0; // Number of C nodes in E-set
    data.number_of_e_edges = 0; // Number of Edges in E-set

    data.number_of_used_ab_cycles = 0;
    for (int s1 = 0; s1 < data.number_of_ab_cycles; ++s1) {
        data.used_ab_cycle[s1] = 0;
        data.weight_sr[s1] = 0;
        data.moved_ab_cycle[s1] = 0;
    }

    for (int s = 0; s < data.number_of_ab_cycles_in_eset; ++s) {
        jnum = data.ab_cycle_in_eset[s];
        add_ab(data, jnum);
    }
    data.best_number_of_c_nodes = data.number_of_c_nodes;
    data.best_number_of_e_edges = data.number_of_e_edges;

    stagnation = 0;
    iteration = 0;
    while (1) {
        ++iteration;
        min_delta_weight_non_tabu = 99999999;
        improving_change = 0;
        non_tabu_change = 0;
        for (int s1 = 0; s1 < data.number_of_ab_cycles; ++s1) {
            if (data.used_ab_cycle[s1] == 0 && data.weight_sr[s1] > 0) {
                delta_weight = data.weight_c[s1] - 2 * data.weight_sr[s1];
                if (data.number_of_c_nodes + delta_weight < data.best_number_of_c_nodes) {
                    selected_ab_cycle = s1;
                    improving_change = 1;
                    data.best_number_of_c_nodes = data.number_of_c_nodes + delta_weight;
                }
                if (delta_weight < min_delta_weight_non_tabu && iteration > data.moved_ab_cycle[s1]) {
                    selected_ab_cycle_non_tabu = s1;
                    non_tabu_change = 1;
                    min_delta_weight_non_tabu = delta_weight;
                }
            } else if (data.used_ab_cycle[s1] == 1 && s1 != center_ab) {
                delta_weight = - data.weight_c[s1] + 2 * data.weight_sr[s1];
                if (data.number_of_c_nodes + delta_weight < data.best_number_of_c_nodes) {
                    selected_ab_cycle = s1;
                    improving_change = -1;
                    data.best_number_of_c_nodes = data.number_of_c_nodes + delta_weight;
                }
                if (delta_weight < min_delta_weight_non_tabu && iteration > data.moved_ab_cycle[s1]) {
                    selected_ab_cycle_non_tabu = s1;
                    non_tabu_change = -1;
                    min_delta_weight_non_tabu = delta_weight;
                }
            }
        }

        if (improving_change != 0) {
            if (improving_change == 1) {
                add_ab(data, selected_ab_cycle);
            } else if (improving_change == -1) {
                delete_ab(data, selected_ab_cycle);
            }

            data.moved_ab_cycle[selected_ab_cycle] = iteration + random_integer(1, data.t_max);

            data.best_number_of_e_edges = data.number_of_e_edges;

            data.number_of_ab_cycles_in_eset = 0;
            for (int s1 = 0; s1 < data.number_of_ab_cycles; ++s1)
                if (data.used_ab_cycle[s1] == 1)
                    data.ab_cycle_in_eset[data.number_of_ab_cycles_in_eset++] = s1;

            stagnation = 0;
        } else if (non_tabu_change != 0) {
            if (non_tabu_change == 1) {
                add_ab(data, selected_ab_cycle_non_tabu);
            } else if (non_tabu_change == -1) {
                delete_ab(data, selected_ab_cycle_non_tabu);
            }
            data.moved_ab_cycle[selected_ab_cycle_non_tabu] = iteration + random_integer(1, data.t_max);
        }
        if (improving_change == 0)
            ++stagnation;
        if (stagnation == data.eset_max_stagnation)
            break;
    }
}

/** Generate 'number_of_kids' children from 'child'/'parent2' and set 'child' to the best one found. */
template <typename Distances>
void run_cross(
        LocalSearchData<Distances>& data,
        Individual& child,
        Individual& parent2,
        int number_of_kids,
        int flag_p)
{
    int number_of_candidates;
    int jnum, center_ab;
    Distance gain;
    Distance best_gain;
    double best_point, point;
    double loss;

    data.evaluation_type = data.flags[0]; // 1:Greedy, 2:---, 3:Distance, 4:Entropy
    data.eset_strategy = data.flags[1]; // 1:Single-AB, 2:Block2

    if (number_of_kids <= data.number_of_ab_cycles) {
        number_of_candidates = number_of_kids;
    } else {
        number_of_candidates = data.number_of_ab_cycles;
    }

    if (data.eset_strategy == 1) { // Single-AB
        random_permutation(data.permutation, data.number_of_ab_cycles, data.number_of_ab_cycles);
    } else if (data.eset_strategy == 2) { // Block2
        for (int k = 0; k < data.number_of_ab_cycles; ++k)
            data.number_of_elements_in_ab_cycle[k] = data.ab_cycle[k][0];
        sort_indices_descending(data.number_of_elements_in_ab_cycle, data.number_of_ab_cycles, data.permutation, data.number_of_ab_cycles);
    }
    data.number_of_generated_children = 0;
    best_point = 0.0;
    best_gain = 0;
    bool improved = false;
    for (int j = 0; j < number_of_candidates; ++j) {
        data.number_of_ab_cycles_in_eset = 0;
        if (data.eset_strategy == 1) { //Single-AB
            jnum = data.permutation[j];
            data.ab_cycle_in_eset[data.number_of_ab_cycles_in_eset++] = jnum;
        } else if (data.eset_strategy == 2) { //Block2
            jnum = data.permutation[j];
            center_ab = jnum;
            for (int s = 0; s < data.number_of_ab_cycles; ++s) {
                if (s == center_ab) {
                    data.ab_cycle_in_eset[data.number_of_ab_cycles_in_eset++] = s;
                } else {
                    if (data.weight_rr[center_ab][s] > 0 && data.ab_cycle[s][0] < data.ab_cycle[center_ab][0]) {
                        if (rand() % 2 == 0)
                            data.ab_cycle_in_eset[data.number_of_ab_cycles_in_eset++] = s;
                    }
                }
            }
            search_eset(data, center_ab);
        }
        data.number_of_segment_positions = 0;
        gain = 0;
        data.number_of_applied_cycles = 0;
        data.number_of_modified_edges = 0;

        data.number_of_applied_cycles = data.number_of_ab_cycles_in_eset;
        for (int k = 0; k < data.number_of_applied_cycles; ++k) {
            data.applied_cycle[k] = data.ab_cycle_in_eset[k];
            jnum = data.applied_cycle[k];
            change_sol(data, child, jnum, flag_p);
            gain += data.gain_ab[jnum];
        }

        make_unit(data);
        make_complete_sol(data, child);
        gain += data.modification_gain;

        ++data.number_of_generated_children;

        if (data.evaluation_type == 1) { // Greedy
            loss = 1.0;
        } else if (data.evaluation_type == 3) { // Distance preservation
            loss = calc_adaptive_loss(data);
        } else if (data.evaluation_type == 4) { // Entropy preservation
            loss = calc_entropy_loss(data);
        }

        if (loss <= 0.0)
            loss = 0.00000001;

        point = (double)gain / loss;
        child.length = child.length - gain;

        if (best_point < point && (2 * data.best_number_of_e_edges < data.distance_ab || child.length != parent2.length)) {
            best_point = point;
            best_gain = gain;
            improved = true;

            data.number_of_best_applied_cycles = data.number_of_applied_cycles;
            for (int s = 0; s < data.number_of_best_applied_cycles; ++s)
                data.best_applied_cycle[s] = data.applied_cycle[s];

            data.number_of_best_modified_edges = data.number_of_modified_edges;
            for (int s = 0; s < data.number_of_best_modified_edges; ++s) {
                data.best_modified_edge[s][0] = data.modified_edge[s][0];
                data.best_modified_edge[s][1] = data.modified_edge[s][1];
                data.best_modified_edge[s][2] = data.modified_edge[s][2];
                data.best_modified_edge[s][3] = data.modified_edge[s][3];
            }

        }
        back_to_pa1(data, child);
        child.length = child.length + gain;
    }
    if (improved) {
        go_to_best(data, child);
        child.length = child.length - best_gain;
        increment_edge_freq(data);
    }
}

/** Number of "C nodes" (unmatched half-edges) if only naively combining the used AB-cycles. */
template <typename Distances>
int calc_c_naive(
        LocalSearchData<Distances>& data)
{
    int count_c_nodes;
    int tally;

    count_c_nodes = 0;

    for (VertexId vertex_id = 0; vertex_id < data.number_of_vertices; ++vertex_id) {
        if (data.in_effect_node[vertex_id][0] != -1 && data.in_effect_node[vertex_id][1] != -1) {
            tally = 0;
            if (data.used_ab_cycle[data.in_effect_node[vertex_id][0]] == 1)
                ++tally;
            if (data.used_ab_cycle[data.in_effect_node[vertex_id][1]] == 1)
                ++tally;
            if (tally == 1)
                ++count_c_nodes;
        }
    }
    return count_c_nodes;
}

/** Reset the generation counters and eset-selection strategy at the start of a run. */
template <typename Distances>
void reset_state(
        LocalSearchData<Distances>& data)
{
    data.accumulated_number_of_children = 0;
    data.current_number_of_generations = 0;
    data.stagnation_count = 0;
    data.max_stagnation = 0;
    data.stage = 1; // sets stage to 1
    data.flags[0] = 4; // maintains population diversity  1:Greedy, 2:---, 3:Distance, 4:Entropy
    data.flags[1] = 1; // the type of Eset: 1:Single-AB, 2:Block2
}

/** Whether the algorithm should stop (time limit, population convergence, or stagnation). */
template <typename Distances>
bool termination_condition(
        LocalSearchData<Distances>& data)
{
    if (data.parameters.timer.needs_to_end())
        return true;
    if (data.average_value - data.best_value < 0.001)
        return true;
    if (data.stage == 1) {
        if (data.stagnation_count == int(1500 / data.number_of_children) && data.max_stagnation == 0) { // 1500/Nch
            data.max_stagnation = int(data.current_number_of_generations / 10); // data.max_stagnation = G/10
        } else if (data.max_stagnation != 0 && data.max_stagnation <= data.stagnation_count) {
            data.stagnation_count = 0;
            data.max_stagnation = 0;
            data.number_of_generations_stage_1 = data.current_number_of_generations;
            data.flags[1] = 2;
            data.stage = 2;
        }
        return false;
    }
    if (data.stage == 2) {
        if (data.stagnation_count == int(1500 / data.number_of_children) && data.max_stagnation == 0) { // 1500/Nch
            data.max_stagnation = int((data.current_number_of_generations - data.number_of_generations_stage_1) / 10); // data.max_stagnation = G/10
        } else if (data.max_stagnation != 0 && data.max_stagnation <= data.stagnation_count) {
            return true;
        }
        return false;
    }

    return true;
}

/** Update 'average_value'/'best_value'/'best_individual' from the current population. */
template <typename Distances>
void set_average_best(
        LocalSearchData<Distances>& data)
{
    Distance stock_best = data.best_individual.length;
    data.average_value = 0.0;
    data.best_index = 0;
    data.best_value = data.population[0].length;
    for (int i = 0; i < data.population_size; ++i) {
        data.average_value += data.population[i].length;
        if (data.population[i].length < data.best_value) {
            data.best_index = i;
            data.best_value = data.population[i].length;
        }
    }
    data.best_individual = data.population[data.best_index];
    data.average_value /= (double)data.population_size;
    if (data.best_individual.length < stock_best) {
        data.stagnation_count = 0;
    } else {
        ++data.stagnation_count;
    }
}

/** Set the population to random tours, then locally optimize each of them. */
template <typename Distances>
void init_population(
        LocalSearchData<Distances>& data)
{
    for (int i = 0; i < data.population_size; ++i) {
        // Always build at least the first individual regardless of the time
        // budget, so there is always a solution to report.
        if (i > 0 && data.parameters.timer.needs_to_end())
            break;

        make_random_solution(data, data.population[i]); // randomly sets a data.route
        run_kopt(data, data.population[i]); // local search (2-opt neighborhood)

        if (!data.output.solution.feasible()
                || data.population[i].length < data.output.solution.objective_value()) {
            data.algorithm_formatter.update_solution(
                    to_solution(data, data.population[i]),
                    "individual " + std::to_string(i));
        }
    }
}

/** Draw a random pairing of the population for crossover. */
template <typename Distances>
void select_for_mating(
        LocalSearchData<Distances>& data)
{
    random_permutation(data.index_for_mating, data.population_size, data.population_size);
    data.index_for_mating[data.population_size] = data.index_for_mating[0];
}

/** Cross the pair of individuals at 'index_for_mating[s]'/'[s + 1]'. */
template <typename Distances>
void generate_kids(
        LocalSearchData<Distances>& data,
        int s)
{
    // data.population[data.index_for_mating[s]] gets replaced by the best solution found by 'run_cross()'
    // 'data.edge_frequency' gets updated at the same time
    set_parents(data, data.population[data.index_for_mating[s]], data.population[data.index_for_mating[s + 1]], data.number_of_children);
    run_cross(data, data.population[data.index_for_mating[s]], data.population[data.index_for_mating[s + 1]], data.number_of_children, 1);
    data.accumulated_number_of_children += data.number_of_generated_children;
}

/** Recompute 'edge_frequency' from the current population. */
template <typename Distances>
void compute_edge_frequencies(
        LocalSearchData<Distances>& data)
{
    for (VertexId vertex_id_1 = 0; vertex_id_1 < data.number_of_vertices; ++vertex_id_1)
        for (VertexId vertex_id_2 = 0; vertex_id_2 < data.number_of_vertices; ++vertex_id_2)
            data.edge_frequency[vertex_id_1][vertex_id_2] = 0;

    for (int i = 0; i < data.population_size; ++i)
        for (VertexId vertex_id = 0; vertex_id < data.number_of_vertices; ++vertex_id) {
            VertexId neighbor_vertex_id_1 = data.population[i].neighbors[vertex_id][0];
            VertexId neighbor_vertex_id_2 = data.population[i].neighbors[vertex_id][1];
            ++data.edge_frequency[vertex_id][neighbor_vertex_id_1];
            ++data.edge_frequency[vertex_id][neighbor_vertex_id_2];
        }
}
}

template <typename Distances>
const Output local_search(
        const Distances& distances,
        const Instance& instance,
        const LocalSearchParameters& parameters)
{
    Output output(instance);
    AlgorithmFormatter algorithm_formatter(parameters, output);
    algorithm_formatter.start("Local Search");
    algorithm_formatter.print_header();

    VertexId number_of_vertices = instance.number_of_vertices();

    // The population/local-search RNG is seeded here (matching what
    // 'InitURandom(seed)' did in the original vendored code); the sort/RNG
    // helper objects the original code lazily allocated as process-wide
    // globals ('tRand'/'tSort') held no state of their own, so this
    // integration replaced them with plain functions ('random_*'/'sort_*'
    // above) instead of carrying that non-reentrancy hazard forward.
    seed_random(parameters.seed);

    LocalSearchData<Distances> data(distances, instance, parameters, output, algorithm_formatter, number_of_vertices);
    init(data, parameters.population_size, parameters.number_of_children);

    init_population(data);

    if (!data.parameters.timer.needs_to_end()) {
        reset_state(data);

        // The genetic algorithm's main loop: maintains a population of
        // individuals, combining pairs of them with 'run_cross()' and locally
        // optimizing the result with 'run_kopt()' each generation.
        compute_edge_frequencies(data);
        while (true) {
            set_average_best(data);
            if (!data.output.solution.feasible()
                    || data.best_individual.length < data.output.solution.objective_value()) {
                data.algorithm_formatter.update_solution(
                        to_solution(data, data.best_individual),
                        "generation " + std::to_string(data.current_number_of_generations));
            }
            if (termination_condition(data))
                break;

            select_for_mating(data);
            for (int s = 0; s < data.population_size; ++s)
                generate_kids(data, s);

            ++data.current_number_of_generations;
        }
    }

    algorithm_formatter.end();
    return output;
}

}
