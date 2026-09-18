#pragma once

#include "travelingsalesmansolver/solution.hpp"
#include "travelingsalesmansolver/algorithm_formatter.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <ctime>
#include <limits>
#include <vector>

namespace travelingsalesmansolver
{

struct EaxParameters: Parameters
{
    /** Population size. */
    int population_size = 100;

    /** Number of children generated per generation. */
    int number_of_children = 30;

    /** Seed of the random number generator. */
    int seed = 0;
};

const Output eax(
        const Instance& instance,
        const EaxParameters& parameters = {});

template <typename Distances>
const Output eax(
        const Distances& distances,
        const Instance& instance,
        const EaxParameters& parameters = {});

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/**
 * Implementation of the Edge Assembly Crossover (EAX) genetic algorithm for
 * the TSP, adapted from Shujia Liu's C++ implementation
 * (https://github.com/Sugia/GA-for-TSP); see 'licenses/eax-ga/NOTICE.md' for
 * the full list of changes made to the original source.
 */
namespace
{

/** A TSP tour, represented as a doubly-linked list over its vertices. */
class Individual
{

public:

    /** Create an individual with no vertices. */
    Individual() { }

    /** Create an individual for a given number of vertices. */
    Individual(VertexId number_of_vertices):
        neighbors(number_of_vertices) { }

    /** Checks if two individuals represent the same tour. */
    bool operator==(
            const Individual& individual) const;

    /** 'neighbors[vertex_id]' holds the two vertices adjacent to 'vertex_id'. */
    std::vector<std::array<VertexId, 2>> neighbors;

    /** Length of the tour. */
    Distance length = 0;

};

/** Seed the random number generator used throughout the algorithm. */
void seed_random(int seed);

/** Random integer in ['min', 'max']. */
int random_integer(int min, int max);

/** Random number drawn from the normal distribution of mean 'mu' and standard deviation 'sigma'. */
double random_normal(double mu, double sigma);

/**
 * Fill 'array' with 'number_of_samples' distinct random values from
 * [0, number_of_elements[.
 */
void random_permutation(
        std::vector<int>& array,
        int number_of_elements,
        int number_of_samples);

/** Randomly shuffle the first 'number_of_elements' elements of 'array'. */
void random_shuffle(
        std::vector<int>& array,
        int number_of_elements);

/**
 * Fill the first 'number_of_indices' elements of 'sorted_indices' with the
 * indices, among the first 'number_of_values' elements of 'values', of the
 * 'number_of_indices' smallest ones, in increasing order of value.
 */
void sort_indices_ascending(
        const std::vector<int>& values,
        int number_of_values,
        std::vector<int>& sorted_indices,
        int number_of_indices);

/**
 * Fill the first 'number_of_indices' elements of 'sorted_indices' with the
 * indices, among the first 'number_of_values' elements of 'values', of the
 * 'number_of_indices' largest ones, in decreasing order of value.
 */
void sort_indices_descending(
        const std::vector<int>& values,
        int number_of_values,
        std::vector<int>& sorted_indices,
        int number_of_indices);

/** Sort the first 'number_of_values' elements of 'values' in increasing order. */
void sort_ascending(
        std::vector<int>& values,
        int number_of_values);

/**
 * Precomputed nearest-neighbor lists for a set of vertices.
 *
 * Distances are not stored/copied here: 'distance(i, j)' is looked up
 * directly on the original (external) 'Distances' object, dispatched once
 * per 'eax()' call via the 'Distances' template parameter, the same way
 * every other algorithm in this library consumes distances.
 */
template <typename Distances>
class Evaluator
{

public:

    /** Constructor. */
    Evaluator(
            const Distances& distances,
            VertexId number_of_vertices);

    /** Compute and store the length of an individual's tour. */
    void evaluate(
            Individual& individual) const;

    /** Return the individual's tour as a 0-indexed list of vertices. */
    std::vector<VertexId> get_tour(
            const Individual& individual) const;

    /** Distance between two vertices. */
    inline Distance distance(
            VertexId vertex_id_1,
            VertexId vertex_id_2) const
    {
        return distances_.distance(vertex_id_1, vertex_id_2);
    }

    /** Number of vertices. */
    VertexId number_of_vertices;

    /** 'near_cities[vertex_id][k]' is the k-th nearest vertex to 'vertex_id'. */
    std::vector<std::vector<VertexId>> near_cities;

private:

    /** Compute 'near_cities' from 'distance()'. */
    void compute_near_cities();

    /** Distances between vertices. */
    const Distances& distances_;

    /** Number of nearest neighbors stored per vertex in 'near_cities'. */
    static constexpr int max_near_cities_ = 50;

};

template <typename Distances>
Evaluator<Distances>::Evaluator(
        const Distances& distances,
        VertexId number_of_vertices):
    number_of_vertices(number_of_vertices),
    near_cities(number_of_vertices, std::vector<VertexId>(max_near_cities_ + 1)),
    distances_(distances)
{
    compute_near_cities();
}

template <typename Distances>
void Evaluator<Distances>::compute_near_cities()
{
    std::vector<int> checked(number_of_vertices);
    for (VertexId vertex_id = 0; vertex_id < number_of_vertices; ++vertex_id) {
        std::fill(checked.begin(), checked.end(), 0);
        checked[vertex_id] = 1;
        near_cities[vertex_id][0] = vertex_id;
        for (int k = 1; k <= max_near_cities_; ++k) {
            VertexId closest_vertex_id = -1;
            Distance min_distance = std::numeric_limits<Distance>::max();
            for (VertexId other_vertex_id = 0; other_vertex_id < number_of_vertices; ++other_vertex_id) {
                if (checked[other_vertex_id] == 0
                        && distance(vertex_id, other_vertex_id) <= min_distance) {
                    closest_vertex_id = other_vertex_id;
                    min_distance = distance(vertex_id, other_vertex_id);
                }
            }
            near_cities[vertex_id][k] = closest_vertex_id;
            checked[closest_vertex_id] = 1;
        }
    }
}

template <typename Distances>
void Evaluator<Distances>::evaluate(
        Individual& individual) const
{
    Distance d = 0;
    for (VertexId vertex_id = 0; vertex_id < number_of_vertices; ++vertex_id)
        d += distance(vertex_id, individual.neighbors[vertex_id][0]) + distance(vertex_id, individual.neighbors[vertex_id][1]);
    individual.length = d / 2;
}

template <typename Distances>
std::vector<VertexId> Evaluator<Distances>::get_tour(
        const Individual& individual) const
{
    std::vector<VertexId> tour(number_of_vertices);
    VertexId current_vertex_id = 0;
    VertexId start_vertex_id = 0;
    VertexId previous_vertex_id = -1;
    for (int count = 0; count < number_of_vertices; ++count) {
        tour[count] = current_vertex_id;
        VertexId next_vertex_id = (individual.neighbors[current_vertex_id][0] == previous_vertex_id)?
            individual.neighbors[current_vertex_id][1]:
            individual.neighbors[current_vertex_id][0];
        previous_vertex_id = current_vertex_id;
        current_vertex_id = next_vertex_id;
        if (current_vertex_id == start_vertex_id)
            break;
    }
    return tour;
}

/**
 * Local search based on the 2-opt neighborhood, operating on a segment-tree
 * representation of the tour to keep 'next'/'previous'/'between' queries and
 * moves cheap even on large instances.
 */
template <typename Distances>
class KOpt
{

public:

    /** Constructor. */
    KOpt(
            Evaluator<Distances>& evaluator,
            VertexId number_of_vertices);

    /** Run the local search on 'individual'. */
    void run(
            Individual& individual);

    /** Set 'individual' to a random tour and run the local search on it. */
    void make_random_solution(
            Individual& individual);

private:

    /** Build the tree representation of 'individual'. */
    void individual_to_tree(
            const Individual& individual);

    /** Rebuild 'individual' from the tree representation. */
    void tree_to_individual(
            Individual& individual) const;

    /** Repeatedly apply improving 2-opt moves until none remain. */
    void optimize();

    /** Vertex immediately after 'vertex_id' in the tour. */
    VertexId next_city(
            VertexId vertex_id) const;

    /** Vertex immediately before 'vertex_id' in the tour. */
    VertexId previous_city(
            VertexId vertex_id) const;

    /** The orientation opposite to 'orientation' (i.e. '1 - orientation'). */
    static int opposite_orientation(
            int orientation)
    {
        return 1 - orientation;
    }

    /** Apply the 2-opt move found by 'optimize()' to the tree representation. */
    void apply_move();

    /** Merge segment 'segment_2' into segment 'segment_1'. */
    void merge_segments(
            int segment_1,
            int segment_2);

    /** Number of vertices. */
    VertexId number_of_vertices_;

    /** Evaluator, used to look up distances and nearest-neighbor lists. */
    Evaluator<Distances>& evaluator_;

    /** For each vertex, the vertices in whose nearest-neighbor list it appears. */
    std::vector<std::vector<VertexId>> inverse_near_list_;

    /** Number of segments fixed at the start of the current 'optimize()' call. */
    int fixed_number_of_segments_ = 0;

    /** Current number of segments. */
    int number_of_segments_ = 0;

    /** Whether the move currently being applied reverses the tour direction. */
    int reversed_ = 0;

    /** Length of the tour, as tracked incrementally by the tree representation. */
    Distance tour_length_ = 0;

    /** 'neighbors_[vertex_id]' holds the two vertices adjacent to 'vertex_id' within its segment ('-1' at a segment end). */
    std::vector<std::array<VertexId, 2>> neighbors_;

    /** 'segment_neighbors_[s]' holds the two segments adjacent to segment 's'. */
    std::vector<std::array<int, 2>> segment_neighbors_;

    /** 'segment_endpoints_[s]' holds the two endpoint vertices of segment 's'. */
    std::vector<std::array<VertexId, 2>> segment_endpoints_;

    /** The four vertices ('t_[1]'..'t_[4]') involved in the move currently being considered. */
    std::array<VertexId, 5> t_;

    /** 'city_segment_[vertex_id]' is the segment 'vertex_id' belongs to. */
    std::vector<int> city_segment_;

    /** Order of the vertices within their segment, used to tell direction/betweenness cheaply. */
    std::vector<int> city_order_;

    /** Order of the segments along the tour. */
    std::vector<int> segment_order_;

    /** Orientation of each segment (0 or 1). */
    std::vector<int> segment_orientation_;

    /** Number of vertices in each segment. */
    std::vector<int> segment_size_;

    /** Whether each vertex is still active (a candidate to start a new improving move from). */
    std::vector<int> active_;

    /** Scratch array used to build 'individual_to_tree()'. */
    std::vector<VertexId> array_;

    /** Scratch array used by 'make_random_solution()'. */
    std::vector<VertexId> remaining_;

    /** Number of nearest neighbors considered per vertex (matches 'Evaluator::max_near_cities_'). */
    static constexpr int max_near_cities_used_ = 50;

};

template <typename Distances>
KOpt<Distances>::KOpt(
        Evaluator<Distances>& evaluator,
        VertexId number_of_vertices):
    number_of_vertices_(number_of_vertices),
    evaluator_(evaluator),
    inverse_near_list_(number_of_vertices),
    neighbors_(number_of_vertices),
    segment_neighbors_(number_of_vertices),
    segment_endpoints_(number_of_vertices),
    t_{},
    city_segment_(number_of_vertices),
    city_order_(number_of_vertices),
    segment_order_(number_of_vertices),
    segment_orientation_(number_of_vertices),
    segment_size_(number_of_vertices),
    active_(number_of_vertices),
    array_(number_of_vertices + 2),
    remaining_(number_of_vertices)
{
    for (VertexId vertex_id = 0; vertex_id < number_of_vertices_; ++vertex_id) {
        for (int k = 0; k < max_near_cities_used_; ++k) {
            VertexId near_vertex_id = evaluator_.near_cities[vertex_id][k];
            inverse_near_list_[near_vertex_id].push_back(vertex_id);
        }
    }
}

template <typename Distances>
void KOpt<Distances>::individual_to_tree(
        const Individual& individual)
{
    array_[1] = 0;
    for (int i = 2; i <= number_of_vertices_; ++i)
        array_[i] = individual.neighbors[array_[i - 1]][1];
    array_[0] = array_[number_of_vertices_];
    array_[number_of_vertices_ + 1] = array_[1];

    int num = 1;
    number_of_segments_ = 0;
    while (true) {
        int orientation = 1;
        int size = 0;
        segment_orientation_[number_of_segments_] = orientation;
        segment_order_[number_of_segments_] = number_of_segments_;

        neighbors_[array_[num]][0] = -1;
        neighbors_[array_[num]][1] = array_[num + 1];
        city_order_[array_[num]] = size;
        city_segment_[array_[num]] = number_of_segments_;
        segment_endpoints_[number_of_segments_][opposite_orientation(orientation)] = array_[num];
        ++num;
        ++size;
        for (int i = 0; i < (int)std::sqrt(number_of_vertices_ * 1.0) - 1; ++i) {
            if (num == number_of_vertices_)
                break;
            neighbors_[array_[num]][0] = array_[num - 1];
            neighbors_[array_[num]][1] = array_[num + 1];
            city_order_[array_[num]] = size;
            city_segment_[array_[num]] = number_of_segments_;
            ++num;
            ++size;
        }
        if (num == number_of_vertices_ - 1) {
            neighbors_[array_[num]][0] = array_[num - 1];
            neighbors_[array_[num]][1] = array_[num + 1];
            city_order_[array_[num]] = size;
            city_segment_[array_[num]] = number_of_segments_;
            ++num;
            ++size;
        }
        neighbors_[array_[num]][0] = array_[num - 1];
        neighbors_[array_[num]][1] = -1;
        city_order_[array_[num]] = size;
        city_segment_[array_[num]] = number_of_segments_;
        segment_endpoints_[number_of_segments_][orientation] = array_[num];
        ++num;
        ++size;
        segment_size_[number_of_segments_] = size;
        ++number_of_segments_;
        if (num == number_of_vertices_ + 1)
            break;
    }
    for (int s = 1; s < number_of_segments_ - 1; ++s) {
        segment_neighbors_[s][0] = s - 1;
        segment_neighbors_[s][1] = s + 1;
    }
    segment_neighbors_[0][0] = number_of_segments_ - 1;
    segment_neighbors_[0][1] = 1;
    segment_neighbors_[number_of_segments_ - 1][0] = number_of_segments_ - 2;
    segment_neighbors_[number_of_segments_ - 1][1] = 0;
    tour_length_ = individual.length;
    fixed_number_of_segments_ = number_of_segments_;
}

template <typename Distances>
void KOpt<Distances>::tree_to_individual(
        Individual& individual) const
{
    for (VertexId vertex_id = 0; vertex_id < number_of_vertices_; ++vertex_id) {
        individual.neighbors[vertex_id][0] = previous_city(vertex_id);
        individual.neighbors[vertex_id][1] = next_city(vertex_id);
    }
    evaluator_.evaluate(individual);
}

template <typename Distances>
void KOpt<Distances>::run(
        Individual& individual)
{
    individual_to_tree(individual);
    optimize();
    tree_to_individual(individual);
}

template <typename Distances>
void KOpt<Distances>::optimize()
{
    std::fill(active_.begin(), active_.end(), 1);
BEGIN:
    {
        VertexId t1_start = random_integer(0, number_of_vertices_ - 1);
        t_[1] = t1_start;
        while (true) {
            t_[1] = next_city(t_[1]);
            if (active_[t_[1]] == 0)
                goto RETURN;
            reversed_ = 0;
            t_[2] = previous_city(t_[1]);
            for (int num1 = 1; num1 < max_near_cities_used_; ++num1) {
                t_[4] = evaluator_.near_cities[t_[1]][num1];
                t_[3] = previous_city(t_[4]);
                Distance dis1 = evaluator_.distance(t_[1], t_[2]) - evaluator_.distance(t_[1], t_[4]);
                if (dis1 > 0) {
                    Distance dis2 = dis1 + evaluator_.distance(t_[3], t_[4]) - evaluator_.distance(t_[3], t_[2]);
                    if (dis2 > 0) {
                        apply_move();
                        for (int a = 1; a <= 4; ++a)
                            for (VertexId near_vertex : inverse_near_list_[t_[a]])
                                active_[near_vertex] = 1;
                        goto BEGIN;
                    }
                } else {
                    break;
                }
            }
            reversed_ = 1;
            t_[2] = next_city(t_[1]);
            for (int num1 = 1; num1 < max_near_cities_used_; ++num1) {
                t_[4] = evaluator_.near_cities[t_[1]][num1];
                t_[3] = next_city(t_[4]);
                Distance dis1 = evaluator_.distance(t_[1], t_[2]) - evaluator_.distance(t_[1], t_[4]);
                if (dis1 > 0) {
                    Distance dis2 = dis1 + evaluator_.distance(t_[3], t_[4]) - evaluator_.distance(t_[3], t_[2]);
                    if (dis2 > 0) {
                        apply_move();
                        for (int a = 1; a <= 4; ++a)
                            for (VertexId near_vertex : inverse_near_list_[t_[a]])
                                active_[near_vertex] = 1;
                        goto BEGIN;
                    }
                } else {
                    break;
                }
            }
            active_[t_[1]] = 0;
RETURN:
            if (t_[1] == t1_start)
                break;
        }
    }
}

template <typename Distances>
VertexId KOpt<Distances>::next_city(
        VertexId vertex_id) const
{
    int seg = city_segment_[vertex_id];
    int orientation = segment_orientation_[seg];
    VertexId next_vertex_id = neighbors_[vertex_id][orientation];
    if (next_vertex_id == -1) {
        seg = segment_neighbors_[seg][orientation];
        orientation = opposite_orientation(segment_orientation_[seg]);
        next_vertex_id = segment_endpoints_[seg][orientation];
    }
    return next_vertex_id;
}

template <typename Distances>
VertexId KOpt<Distances>::previous_city(
        VertexId vertex_id) const
{
    int seg = city_segment_[vertex_id];
    int orientation = segment_orientation_[seg];
    VertexId previous_vertex_id = neighbors_[vertex_id][opposite_orientation(orientation)];
    if (previous_vertex_id == -1) {
        seg = segment_neighbors_[seg][opposite_orientation(orientation)];
        orientation = segment_orientation_[seg];
        previous_vertex_id = segment_endpoints_[seg][orientation];
    }
    return previous_vertex_id;
}

template <typename Distances>
void KOpt<Distances>::apply_move()
{
    VertexId t1_s, t1_e, t2_s, t2_e;

    if (reversed_ == 0) {
        t1_s = t_[1]; t1_e = t_[3]; t2_s = t_[4]; t2_e = t_[2];
    } else {
        t1_s = t_[2]; t1_e = t_[4]; t2_s = t_[3]; t2_e = t_[1];
    }

    int seg_t1_s = city_segment_[t1_s];
    int ordSeg_t1_s = segment_order_[seg_t1_s];
    int orient_t1_s = segment_orientation_[seg_t1_s];
    int seg_t1_e = city_segment_[t1_e];
    int ordSeg_t1_e = segment_order_[seg_t1_e];
    int orient_t1_e = segment_orientation_[seg_t1_e];
    int seg_t2_s = city_segment_[t2_s];
    int ordSeg_t2_s = segment_order_[seg_t2_s];
    int orient_t2_s = segment_orientation_[seg_t2_s];
    int seg_t2_e = city_segment_[t2_e];
    int ordSeg_t2_e = segment_order_[seg_t2_e];
    int orient_t2_e = segment_orientation_[seg_t2_e];

    //////////////////// Type1 ////////////////////////
    if ((seg_t1_s == seg_t1_e) && (seg_t1_s == seg_t2_s) && (seg_t1_s == seg_t2_e)) {
        if ((segment_orientation_[seg_t1_s] == 1 && (city_order_[t1_s] > city_order_[t1_e])) ||
            (segment_orientation_[seg_t1_s] == 0 && (city_order_[t1_s] < city_order_[t1_e]))) {
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
        int ord = city_order_[t1_e];
        while (true) {
            std::swap(neighbors_[curr][0], neighbors_[curr][1]);
            city_order_[curr] = ord;
            if (curr == t1_e)
                break;
            curr = neighbors_[curr][opposite_orientation(orient_t1_s)];
            if (orient_t1_s == 0)
                ++ord;
            else
                --ord;
        }

        neighbors_[t2_e][orient_t1_s] = t1_e;
        neighbors_[t2_s][opposite_orientation(orient_t1_s)] = t1_s;
        neighbors_[t1_s][orient_t1_s] = t2_s;
        neighbors_[t1_e][opposite_orientation(orient_t1_s)] = t2_e;

        return;
    }
    //////////////////// Type1 ///////////////////////

    int numOfSeg1 = (ordSeg_t1_e >= ordSeg_t1_s) ?
        ordSeg_t1_e - ordSeg_t1_s + 1 :
        ordSeg_t1_e - ordSeg_t1_s + 1 + number_of_segments_;
    int numOfSeg2 = (ordSeg_t2_e >= ordSeg_t2_s) ?
        ordSeg_t2_e - ordSeg_t2_s + 1 :
        ordSeg_t2_e - ordSeg_t2_s + 1 + number_of_segments_;

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
    int flag_t2e_t1s = (neighbors_[t2_e][orient_t2_e] == -1) ? 1 : 0;
    int flag_t2s_t1e = (neighbors_[t2_s][opposite_orientation(orient_t2_s)] == -1) ? 1 : 0;

    int length_t1s_seg = std::abs(city_order_[t2_e] - city_order_[segment_endpoints_[seg_t2_e][orient_t2_e]]);
    int length_t1e_seg = std::abs(city_order_[t2_s] - city_order_[segment_endpoints_[seg_t2_s][opposite_orientation(orient_t2_s)]]);

    ///////////////////// Type2 /////////////////
    if (seg_t1_s == seg_t1_e) {
        if (flag_t2e_t1s == 1 && flag_t2s_t1e == 1) {
            orient_t1_s = opposite_orientation(segment_orientation_[seg_t1_s]);
            segment_orientation_[seg_t1_s] = orient_t1_s;
            segment_endpoints_[seg_t1_s][orient_t1_s] = t1_s;
            segment_endpoints_[seg_t1_s][opposite_orientation(orient_t1_s)] = t1_e;
            segment_neighbors_[seg_t1_s][orient_t1_s] = seg_t2_s;
            segment_neighbors_[seg_t1_s][opposite_orientation(orient_t1_s)] = seg_t2_e;
            return;
        }
        if (flag_t2e_t1s == 0 && flag_t2s_t1e == 1) {
            VertexId curr = t1_e;
            int ord = city_order_[t1_s];
            while (true) {
                std::swap(neighbors_[curr][0], neighbors_[curr][1]);
                city_order_[curr] = ord;
                if (curr == t1_s)
                    break;
                curr = neighbors_[curr][orient_t2_e];
                if (orient_t2_e == 0)
                    --ord;
                else
                    ++ord;
            }
            neighbors_[t2_e][orient_t2_e] = t1_e;
            neighbors_[t1_s][orient_t2_e] = -1;
            neighbors_[t1_e][opposite_orientation(orient_t2_e)] = t2_e;
            segment_endpoints_[seg_t2_e][orient_t2_e] = t1_s;
            return;
        }
        if (flag_t2e_t1s == 1 && flag_t2s_t1e == 0) {
            VertexId curr = t1_s;
            int ord = city_order_[t1_e];
            while (true) {
                std::swap(neighbors_[curr][0], neighbors_[curr][1]);
                city_order_[curr] = ord;
                if (curr == t1_e)
                    break;
                curr = neighbors_[curr][opposite_orientation(orient_t2_s)];
                if (orient_t2_s == 0)
                    ++ord;
                else
                    --ord;
            }
            neighbors_[t2_s][opposite_orientation(orient_t2_s)] = t1_s;
            neighbors_[t1_e][opposite_orientation(orient_t2_s)] = -1;
            neighbors_[t1_s][orient_t2_s] = t2_s;
            segment_endpoints_[seg_t2_s][opposite_orientation(orient_t2_s)] = t1_e;
            return;
        }
    }

    ///////////////////// Type3 /////////////////

    if (flag_t2e_t1s == 1) {
        segment_neighbors_[seg_t1_s][opposite_orientation(orient_t1_s)] = seg_t2_s;
    } else {
        seg_t1_s = number_of_segments_++;
        orient_t1_s = orient_t2_e;
        neighbors_[t1_s][opposite_orientation(orient_t1_s)] = -1;
        neighbors_[segment_endpoints_[seg_t2_e][orient_t2_e]][orient_t1_s] = -1;
        segment_orientation_[seg_t1_s] = orient_t1_s;
        segment_size_[seg_t1_s] = length_t1s_seg;
        segment_endpoints_[seg_t1_s][opposite_orientation(orient_t1_s)] = t1_s;
        segment_endpoints_[seg_t1_s][orient_t1_s] = segment_endpoints_[seg_t2_e][orient_t2_e];
        segment_neighbors_[seg_t1_s][opposite_orientation(orient_t1_s)] = seg_t2_s;
        segment_neighbors_[seg_t1_s][orient_t1_s] = segment_neighbors_[seg_t2_e][orient_t2_e];
        int seg = segment_neighbors_[seg_t2_e][orient_t2_e];
        segment_neighbors_[seg][opposite_orientation(segment_orientation_[seg])] = seg_t1_s;
    }

    if (flag_t2s_t1e == 1) {
        segment_neighbors_[seg_t1_e][orient_t1_e] = seg_t2_e;
    } else {
        seg_t1_e = number_of_segments_++;
        orient_t1_e = orient_t2_s;
        neighbors_[t1_e][orient_t1_e] = -1;
        neighbors_[segment_endpoints_[seg_t2_s][opposite_orientation(orient_t2_s)]][opposite_orientation(orient_t1_e)] = -1;
        segment_orientation_[seg_t1_e] = orient_t1_e;
        segment_size_[seg_t1_e] = length_t1e_seg;
        segment_endpoints_[seg_t1_e][orient_t1_e] = t1_e;
        segment_endpoints_[seg_t1_e][opposite_orientation(orient_t1_e)] = segment_endpoints_[seg_t2_s][opposite_orientation(orient_t2_s)];
        segment_neighbors_[seg_t1_e][orient_t1_e] = seg_t2_e;
        segment_neighbors_[seg_t1_e][opposite_orientation(orient_t1_e)] = segment_neighbors_[seg_t2_s][opposite_orientation(orient_t2_s)];
        int seg = segment_neighbors_[seg_t2_s][opposite_orientation(orient_t2_s)];
        segment_neighbors_[seg][segment_orientation_[seg]] = seg_t1_e;
    }

    neighbors_[t2_e][orient_t2_e] = -1;
    segment_size_[seg_t2_e] -= length_t1s_seg;
    segment_endpoints_[seg_t2_e][orient_t2_e] = t2_e;
    segment_neighbors_[seg_t2_e][orient_t2_e] = seg_t1_e;
    neighbors_[t2_s][opposite_orientation(orient_t2_s)] = -1;
    segment_size_[seg_t2_s] -= length_t1e_seg;
    segment_endpoints_[seg_t2_s][opposite_orientation(orient_t2_s)] = t2_s;
    segment_neighbors_[seg_t2_s][opposite_orientation(orient_t2_s)] = seg_t1_s;

    {
        int seg = seg_t1_e;
        while (true) {
            segment_orientation_[seg] = opposite_orientation(segment_orientation_[seg]);
            if (seg == seg_t1_s)
                break;
            seg = segment_neighbors_[seg][segment_orientation_[seg]];
        }
    }

    if (segment_size_[seg_t2_e] < length_t1s_seg) {
        int seg = segment_neighbors_[seg_t2_e][opposite_orientation(segment_orientation_[seg_t2_e])];
        segment_neighbors_[seg][segment_orientation_[seg]] = seg_t1_s;
        seg = segment_neighbors_[seg_t2_e][segment_orientation_[seg_t2_e]];
        segment_neighbors_[seg][opposite_orientation(segment_orientation_[seg])] = seg_t1_s;
        seg = segment_neighbors_[seg_t1_s][opposite_orientation(segment_orientation_[seg_t1_s])];
        segment_neighbors_[seg][segment_orientation_[seg]] = seg_t2_e;
        seg = segment_neighbors_[seg_t1_s][segment_orientation_[seg_t1_s]];
        segment_neighbors_[seg][opposite_orientation(segment_orientation_[seg])] = seg_t2_e;

        std::swap(segment_orientation_[seg_t2_e], segment_orientation_[seg_t1_s]);
        std::swap(segment_size_[seg_t2_e], segment_size_[seg_t1_s]);
        std::swap(segment_endpoints_[seg_t2_e][0], segment_endpoints_[seg_t1_s][0]);
        std::swap(segment_endpoints_[seg_t2_e][1], segment_endpoints_[seg_t1_s][1]);
        std::swap(segment_neighbors_[seg_t2_e][0], segment_neighbors_[seg_t1_s][0]);
        std::swap(segment_neighbors_[seg_t2_e][1], segment_neighbors_[seg_t1_s][1]);
        std::swap(seg_t2_e, seg_t1_s);
    }

    if (segment_size_[seg_t2_s] < length_t1e_seg) {
        int seg = segment_neighbors_[seg_t2_s][opposite_orientation(segment_orientation_[seg_t2_s])];
        segment_neighbors_[seg][segment_orientation_[seg]] = seg_t1_e;
        seg = segment_neighbors_[seg_t2_s][segment_orientation_[seg_t2_s]];
        segment_neighbors_[seg][opposite_orientation(segment_orientation_[seg])] = seg_t1_e;
        seg = segment_neighbors_[seg_t1_e][opposite_orientation(segment_orientation_[seg_t1_e])];
        segment_neighbors_[seg][segment_orientation_[seg]] = seg_t2_s;
        seg = segment_neighbors_[seg_t1_e][segment_orientation_[seg_t1_e]];
        segment_neighbors_[seg][opposite_orientation(segment_orientation_[seg])] = seg_t2_s;

        std::swap(segment_orientation_[seg_t2_s], segment_orientation_[seg_t1_e]);
        std::swap(segment_size_[seg_t2_s], segment_size_[seg_t1_e]);
        std::swap(segment_endpoints_[seg_t2_s][0], segment_endpoints_[seg_t1_e][0]);
        std::swap(segment_endpoints_[seg_t2_s][1], segment_endpoints_[seg_t1_e][1]);
        std::swap(segment_neighbors_[seg_t2_s][0], segment_neighbors_[seg_t1_e][0]);
        std::swap(segment_neighbors_[seg_t2_s][1], segment_neighbors_[seg_t1_e][1]);
        std::swap(seg_t2_s, seg_t1_e);
    }

    while (number_of_segments_ > fixed_number_of_segments_) {
        if (segment_size_[segment_neighbors_[number_of_segments_ - 1][0]] <
                segment_size_[segment_neighbors_[number_of_segments_ - 1][1]])
            merge_segments(segment_neighbors_[number_of_segments_ - 1][0], number_of_segments_ - 1);
        else
            merge_segments(segment_neighbors_[number_of_segments_ - 1][1], number_of_segments_ - 1);
    }
    int ordSeg = 0;
    int seg = 0;
    while (true) {
        segment_order_[seg] = ordSeg;
        ++ordSeg;
        seg = segment_neighbors_[seg][segment_orientation_[seg]];
        if (seg == 0)
            break;
    }
}

template <typename Distances>
void KOpt<Distances>::merge_segments(
        int segment_1,
        int segment_2)
{
    VertexId t_s = 0, t_e = 0;
    int direction = 0, ord = 0, increment = 0;

    if (segment_neighbors_[segment_1][segment_orientation_[segment_1]] == segment_2) {
        neighbors_[segment_endpoints_[segment_1][segment_orientation_[segment_1]]][segment_orientation_[segment_1]] =
            segment_endpoints_[segment_2][opposite_orientation(segment_orientation_[segment_2])];
        neighbors_[segment_endpoints_[segment_2][opposite_orientation(segment_orientation_[segment_2])]][opposite_orientation(segment_orientation_[segment_2])] =
            segment_endpoints_[segment_1][segment_orientation_[segment_1]];
        ord = city_order_[segment_endpoints_[segment_1][segment_orientation_[segment_1]]];

        segment_endpoints_[segment_1][segment_orientation_[segment_1]] = segment_endpoints_[segment_2][segment_orientation_[segment_2]];
        segment_neighbors_[segment_1][segment_orientation_[segment_1]] = segment_neighbors_[segment_2][segment_orientation_[segment_2]];
        int seg = segment_neighbors_[segment_2][segment_orientation_[segment_2]];
        segment_neighbors_[seg][opposite_orientation(segment_orientation_[seg])] = segment_1;

        t_s = segment_endpoints_[segment_2][opposite_orientation(segment_orientation_[segment_2])];
        t_e = segment_endpoints_[segment_2][segment_orientation_[segment_2]];
        direction = segment_orientation_[segment_2];

        increment = (segment_orientation_[segment_1] == 1) ? 1 : -1;
    } else if (segment_neighbors_[segment_1][opposite_orientation(segment_orientation_[segment_1])] == segment_2) {
        neighbors_[segment_endpoints_[segment_1][opposite_orientation(segment_orientation_[segment_1])]][opposite_orientation(segment_orientation_[segment_1])] =
            segment_endpoints_[segment_2][segment_orientation_[segment_2]];
        neighbors_[segment_endpoints_[segment_2][segment_orientation_[segment_2]]][segment_orientation_[segment_2]] =
            segment_endpoints_[segment_1][opposite_orientation(segment_orientation_[segment_1])];
        ord = city_order_[segment_endpoints_[segment_1][opposite_orientation(segment_orientation_[segment_1])]];

        segment_endpoints_[segment_1][opposite_orientation(segment_orientation_[segment_1])] = segment_endpoints_[segment_2][opposite_orientation(segment_orientation_[segment_2])];
        segment_neighbors_[segment_1][opposite_orientation(segment_orientation_[segment_1])] = segment_neighbors_[segment_2][opposite_orientation(segment_orientation_[segment_2])];
        int seg = segment_neighbors_[segment_2][opposite_orientation(segment_orientation_[segment_2])];
        segment_neighbors_[seg][segment_orientation_[seg]] = segment_1;

        t_s = segment_endpoints_[segment_2][segment_orientation_[segment_2]];
        t_e = segment_endpoints_[segment_2][opposite_orientation(segment_orientation_[segment_2])];
        direction = opposite_orientation(segment_orientation_[segment_2]);

        increment = (segment_orientation_[segment_1] == 1) ? -1 : 1;
    }
    VertexId curr = t_s;
    ord = ord + increment;
    while (true) {
        city_segment_[curr] = segment_1;
        city_order_[curr] = ord;

        VertexId next_vertex_id = neighbors_[curr][direction];
        if (segment_orientation_[segment_1] != segment_orientation_[segment_2])
            std::swap(neighbors_[curr][0], neighbors_[curr][1]);

        if (curr == t_e)
            break;
        curr = next_vertex_id;
        ord += increment;
    }
    segment_size_[segment_1] += segment_size_[segment_2];
    --number_of_segments_;
}

template <typename Distances>
void KOpt<Distances>::make_random_solution(
        Individual& individual)
{
    for (VertexId vertex_id = 0; vertex_id < number_of_vertices_; ++vertex_id)
        remaining_[vertex_id] = vertex_id;
    std::vector<VertexId> gene(number_of_vertices_);
    for (int i = 0; i < number_of_vertices_; ++i) {
        int r = random_integer(0, number_of_vertices_ - i - 1);
        gene[i] = remaining_[r];
        remaining_[r] = remaining_[number_of_vertices_ - i - 1];
    }

    for (int j = 1; j < number_of_vertices_ - 1; ++j) {
        individual.neighbors[gene[j]][0] = gene[j - 1];
        individual.neighbors[gene[j]][1] = gene[j + 1];
    }
    individual.neighbors[gene[0]][0] = gene[number_of_vertices_ - 1];
    individual.neighbors[gene[0]][1] = gene[1];
    individual.neighbors[gene[number_of_vertices_ - 1]][0] = gene[number_of_vertices_ - 2];
    individual.neighbors[gene[number_of_vertices_ - 1]][1] = gene[0];

    evaluator_.evaluate(individual);
}

/**
 * The Edge Assembly Crossover (EAX) operator: combines two parent tours
 * into a set of candidate children and keeps the best one found.
 */
template <typename Distances>
class Cross
{

public:

    /** Constructor. */
    Cross(
            Evaluator<Distances>& evaluator,
            VertexId number_of_vertices,
            int population_size);

    /** Generate 'number_of_kids' children from 'child'/'parent2' and set 'child' to the best one found. */
    void run(
            Individual& child,
            Individual& parent2,
            int number_of_kids,
            int flag_p,
            int flags[10],
            std::vector<std::vector<int>>& edge_frequency);

    /** Prepare AB-cycles for a given pair of parents. */
    void set_parents(
            const Individual& parent1,
            const Individual& parent2,
            int flags[10],
            int number_of_kids);

    /** Number of children generated by the last 'run()' call. */
    int number_of_generated_children = 0;

private:

    /** Find the AB-cycles of a given pair of parents. */
    void set_ab_cycle(
            const Individual& parent1,
            const Individual& parent2,
            int flags[10],
            int number_of_kids);

    /** Record the AB-cycle currently being traced by 'set_ab_cycle()'. */
    void form_ab_cycle();

    /** Apply (type == 1) or roll back (type == 2) an AB-cycle on 'child'. */
    void change_sol(
            Individual& child,
            int ab_number,
            int type);

    /** The 5th step of EAX: complete the tour from the remaining path segments. */
    void make_complete_sol(
            Individual& child);

    /** The 5-1th step of EAX: group path segments into units. */
    void make_unit();

    /** Roll back 'child' to parent 1. */
    void back_to_pa1(
            Individual& child);

    /** Set 'child' to the best child found by the last 'run()' call. */
    void go_to_best(
            Individual& child);

    /** Update 'edge_frequency' for the best child found by the last 'run()' call. */
    void increment_edge_freq(
            std::vector<std::vector<int>>& edge_frequency);

    /** Change in the average tour length implied by preserving 'edge_frequency'. */
    int calc_adaptive_loss(
            std::vector<std::vector<int>>& edge_frequency);

    /** Change in edge-frequency entropy implied by 'edge_frequency'. */
    double calc_entropy_loss(
            std::vector<std::vector<int>>& edge_frequency);

    /** Block2 eset selection: weight each AB-cycle by its interaction with the others. */
    void set_weight(
            const Individual& parent1,
            const Individual& parent2);

    /** Number of "C nodes" (unmatched half-edges) if only naively combining the used AB-cycles. */
    int calc_c_naive();

    /** Block2 eset selection: local search over which AB-cycles to include around 'num'. */
    void search_eset(
            int num);

    /** Add AB-cycle 'num' to the eset being built by 'search_eset()'. */
    void add_ab(
            int num);

    /** Remove AB-cycle 'num' from the eset being built by 'search_eset()'. */
    void delete_ab(
            int num);

    /** Evaluator, used to look up distances. */
    Evaluator<Distances>& evaluator_;

    /** Population size. */
    int population_size_;

    /** Number of vertices. */
    VertexId number_of_vertices_;
    int random_pick_;
    int start_new_trace_;
    int cycle_complete_;
    int traversal_type_;
    int number_of_unbranched_;
    int number_of_branching_;
    VertexId trace_start_;
    VertexId current_city_;
    VertexId previous_city_;
    int start_appearance_count_;
    int evaluation_type_;
    int eset_strategy_;
    int number_of_ab_cycles_;
    int position_current_;
    int max_number_of_ab_cycles_;

    std::vector<VertexId> unbranched_;
    std::vector<VertexId> branching_;
    std::vector<int> unbranched_index_;
    std::vector<int> branching_index_;
    std::vector<int> first_visit_position_;
    std::vector<VertexId> route_;
    std::vector<int> permutation_;
    std::vector<VertexId> c_;

    std::vector<std::vector<VertexId>> near_data_;
    std::vector<std::vector<int>> ab_cycle_;

    // speeds up start
    int number_of_units_;
    int number_of_segments_;
    int number_of_segment_positions_;
    int number_of_elements_in_center_unit_;
    int number_of_segments_for_center_;
    Distance modification_gain_;
    int number_of_modified_edges_;
    int number_of_best_modified_edges_;
    int number_of_applied_cycles_;
    int number_of_best_applied_cycles_;

    std::vector<VertexId> order_;
    std::vector<int> inverse_order_;
    std::vector<int> segment_unit_;
    std::vector<int> segment_position_list_;
    std::vector<int> link_a_position_;
    std::vector<int> position_segment_;
    std::vector<int> number_of_elements_in_unit_;
    std::vector<int> center_unit_;
    std::vector<VertexId> list_of_center_unit_;
    std::vector<int> segment_for_center_;
    std::vector<Distance> gain_ab_;
    std::vector<int> applied_cycle_;
    std::vector<int> best_applied_cycle_;

    std::vector<std::vector<int>> segment_;
    std::vector<std::vector<int>> link_b_position_;
    std::vector<std::vector<VertexId>> modified_edge_;
    std::vector<std::vector<VertexId>> best_modified_edge_;
    // speeds up end

    // block2
    int number_of_used_ab_cycles_;
    int number_of_c_nodes_;
    int number_of_e_edges_;
    int t_max_;
    int max_stagnation_;
    int number_of_ab_cycles_in_eset_;
    int distance_ab_;
    int best_number_of_c_nodes_;
    int best_number_of_e_edges_;

    std::vector<int> number_of_elements_in_ab_cycle_;
    std::vector<int> weight_sr_;
    std::vector<int> weight_c_;
    std::vector<int> used_ab_cycle_;
    std::vector<int> moved_ab_cycle_;
    std::vector<int> ab_cycle_in_eset_;

    std::vector<std::vector<int>> in_effect_node_;
    std::vector<std::vector<int>> weight_rr_;
};

template <typename Distances>
Cross<Distances>::Cross(
        Evaluator<Distances>& evaluator,
        VertexId number_of_vertices,
        int population_size):
    evaluator_(evaluator),
    population_size_(population_size)
{
    max_number_of_ab_cycles_ = 2000; // sets the maximum number of ab cycle
    number_of_vertices_ = number_of_vertices;

    near_data_.clear();
    for (int i = 0; i < number_of_vertices_; i++) {
        std::vector<VertexId> row(5);
        near_data_.push_back(row);
    }

    ab_cycle_.clear();
    for (int i = 0; i < max_number_of_ab_cycles_; i++) {
        std::vector<int> row(2 * number_of_vertices_ + 4);
        ab_cycle_.push_back(row);
    }

    unbranched_.resize(number_of_vertices_);
    branching_.resize(number_of_vertices_);
    unbranched_index_.resize(number_of_vertices_);
    branching_index_.resize(number_of_vertices_);
    first_visit_position_.resize(number_of_vertices_);
    route_.resize(2 * number_of_vertices_ + 1);
    permutation_.resize(max_number_of_ab_cycles_);

    c_.resize(2 * number_of_vertices_ + 4);

    // speeds up start
    order_.resize(number_of_vertices_);
    inverse_order_.resize(number_of_vertices_);

    // 'number_of_segment_positions_'/'number_of_segments_' are accumulated across every AB-cycle applied
    // within one crossover call (see 'change_sol()'/'make_unit()'), which in
    // the "Block2" eset mode (multiple AB-cycles per call) can exceed 'number_of_vertices_'
    // even though each individual position is a valid index into a tour of
    // 'number_of_vertices_' cities; size these generously (matching 'route_'/'c_' elsewhere
    // in this file) rather than assuming the count is bounded by 'number_of_vertices_'.
    segment_.clear();
    for (int i = 0; i < 2 * number_of_vertices_; i++) {
        std::vector<int> row(2);
        segment_.push_back(row);
    }

    segment_unit_.resize(2 * number_of_vertices_);
    segment_position_list_.resize(2 * number_of_vertices_);
    link_a_position_.resize(number_of_vertices_);

    link_b_position_.clear();
    for (int i = 0; i < number_of_vertices_; i++) {
        std::vector<int> row(2);
        link_b_position_.push_back(row);
    }

    position_segment_.resize(number_of_vertices_);
    number_of_elements_in_unit_.resize(number_of_vertices_);

    center_unit_.resize(number_of_vertices_);
    for (int i = 0; i < number_of_vertices_; i++) {
        center_unit_[i] = 0;
    }

    list_of_center_unit_.resize(number_of_vertices_ + 2);
    segment_for_center_.resize(number_of_vertices_);
    gain_ab_.resize(number_of_vertices_);

    modified_edge_.clear();
    for (int i = 0; i < number_of_vertices_; i++) {
        std::vector<VertexId> row(4);
        modified_edge_.push_back(row);
    }

    best_modified_edge_.clear();
    for (int i = 0; i < number_of_vertices_; i++) {
        std::vector<VertexId> row(4);
        best_modified_edge_.push_back(row);
    }

    applied_cycle_.resize(number_of_vertices_);
    best_applied_cycle_.resize(number_of_vertices_);
    // Speed Up End

    // block2
    number_of_elements_in_ab_cycle_.resize(max_number_of_ab_cycles_);

    in_effect_node_.clear();
    for (int i = 0; i < number_of_vertices_; i++) {
        std::vector<int> row(2);
        in_effect_node_.push_back(row);
    }

    weight_rr_.clear();
    for (int i = 0; i < max_number_of_ab_cycles_; i++) {
        std::vector<int> row(max_number_of_ab_cycles_);
        weight_rr_.push_back(row);
    }

    weight_sr_.resize(max_number_of_ab_cycles_);
    weight_c_.resize(max_number_of_ab_cycles_);
    used_ab_cycle_.resize(number_of_vertices_);
    moved_ab_cycle_.resize(number_of_vertices_);
    ab_cycle_in_eset_.resize(max_number_of_ab_cycles_);
}

template <typename Distances>
void Cross<Distances>::set_parents(
        const Individual& parent1,
        const Individual& parent2,
        int flags[10],
        int number_of_kids)
{
    this->set_ab_cycle(parent1, parent2, flags, number_of_kids);

    // NOTE: this local deliberately shadows the member 'distance_ab_' under
    // a different name (rather than writing through to it) -- this
    // reproduces a pre-existing upstream quirk byte-for-byte: the member
    // 'distance_ab_', read later in 'run()' (guarding
    // '2 * best_number_of_e_edges_ < distance_ab_'), is never actually
    // updated by this method, so 'run()' always sees a stale value from a
    // previous call (or 0, on the very first call). See NOTICE.md.
    int distance_ab_local = 0;
    VertexId start_vertex_id = 0;
    VertexId current_vertex_id = -1;
    VertexId next_vertex_id = start_vertex_id;
    VertexId previous_vertex_id;
    for (int i = 0; i < number_of_vertices_; ++i) {
        previous_vertex_id = current_vertex_id;
        current_vertex_id = next_vertex_id;
        if (parent1.neighbors[current_vertex_id][0] != previous_vertex_id) {
            next_vertex_id = parent1.neighbors[current_vertex_id][0];
        } else {
            next_vertex_id = parent1.neighbors[current_vertex_id][1];
        }
        if (parent2.neighbors[current_vertex_id][0] != next_vertex_id && parent2.neighbors[current_vertex_id][1] != next_vertex_id)
            ++distance_ab_local;
        order_[i] = current_vertex_id;
        inverse_order_[current_vertex_id] = i;
    }

    if (flags[1] == 2) {
        t_max_ = 10;
        max_stagnation_ = 20; // 1:Greedy LS, 20:Tabu Search
        this->set_weight(parent1, parent2);
    }
}

template <typename Distances>
void Cross<Distances>::run(
        Individual& child,
        Individual& parent2,
        int number_of_kids,
        int flag_p,
        int flags[10],
        std::vector<std::vector<int>>& edge_frequency)
{
    int number_of_candidates;
    int jnum, center_ab;
    Distance gain;
    Distance best_gain;
    double best_point, point;
    double loss;

    evaluation_type_ = flags[0]; // 1:Greedy, 2:---, 3:Distance, 4:Entropy
    eset_strategy_ = flags[1]; // 1:Single-AB, 2:Block2

    if (number_of_kids <= number_of_ab_cycles_) {
        number_of_candidates = number_of_kids;
    } else {
        number_of_candidates = number_of_ab_cycles_;
    }

    if (eset_strategy_ == 1) { // Single-AB
        random_permutation(permutation_, number_of_ab_cycles_, number_of_ab_cycles_);
    } else if (eset_strategy_ == 2) { // Block2
        for (int k = 0; k < number_of_ab_cycles_; ++k)
            number_of_elements_in_ab_cycle_[k] = ab_cycle_[k][0];
        sort_indices_descending(number_of_elements_in_ab_cycle_, number_of_ab_cycles_, permutation_, number_of_ab_cycles_);
    }
    number_of_generated_children = 0;
    best_point = 0.0;
    best_gain = 0;
    bool improved = false;
    for (int j = 0; j < number_of_candidates; ++j) {
        number_of_ab_cycles_in_eset_ = 0;
        if (eset_strategy_ == 1) { //Single-AB
            jnum = permutation_[j];
            ab_cycle_in_eset_[number_of_ab_cycles_in_eset_++] = jnum;
        } else if (eset_strategy_ == 2) { //Block2
            jnum = permutation_[j];
            center_ab = jnum;
            for (int s = 0; s < number_of_ab_cycles_; ++s) {
                if (s == center_ab) {
                    ab_cycle_in_eset_[number_of_ab_cycles_in_eset_++] = s;
                } else {
                    if (weight_rr_[center_ab][s] > 0 && ab_cycle_[s][0] < ab_cycle_[center_ab][0]) {
                        if (rand() % 2 == 0)
                            ab_cycle_in_eset_[number_of_ab_cycles_in_eset_++] = s;
                    }
                }
            }
            this->search_eset(center_ab);
        }
        number_of_segment_positions_ = 0;
        gain = 0;
        number_of_applied_cycles_ = 0;
        number_of_modified_edges_ = 0;

        number_of_applied_cycles_ = number_of_ab_cycles_in_eset_;
        for (int k = 0; k < number_of_applied_cycles_; ++k) {
            applied_cycle_[k] = ab_cycle_in_eset_[k];
            jnum = applied_cycle_[k];
            this->change_sol(child, jnum, flag_p);
            gain += gain_ab_[jnum];
        }

        this->make_unit();
        this->make_complete_sol(child);
        gain += modification_gain_;

        ++number_of_generated_children;

        if (evaluation_type_ == 1) { // Greedy
            loss = 1.0;
        } else if (evaluation_type_ == 3) { // Distance preservation
            loss = this->calc_adaptive_loss(edge_frequency);
        } else if (evaluation_type_ == 4) { // Entropy preservation
            loss = this->calc_entropy_loss(edge_frequency);
        }

        if (loss <= 0.0)
            loss = 0.00000001;

        point = (double)gain / loss;
        child.length = child.length - gain;

        if (best_point < point && (2 * best_number_of_e_edges_ < distance_ab_ || child.length != parent2.length)) {
            best_point = point;
            best_gain = gain;
            improved = true;

            number_of_best_applied_cycles_ = number_of_applied_cycles_;
            for (int s = 0; s < number_of_best_applied_cycles_; ++s)
                best_applied_cycle_[s] = applied_cycle_[s];

            number_of_best_modified_edges_ = number_of_modified_edges_;
            for (int s = 0; s < number_of_best_modified_edges_; ++s) {
                best_modified_edge_[s][0] = modified_edge_[s][0];
                best_modified_edge_[s][1] = modified_edge_[s][1];
                best_modified_edge_[s][2] = modified_edge_[s][2];
                best_modified_edge_[s][3] = modified_edge_[s][3];
            }

        }
        this->back_to_pa1(child);
        child.length = child.length + gain;
    }
    if (improved) {
        this->go_to_best(child);
        child.length = child.length - best_gain;
        this->increment_edge_freq(edge_frequency);
    }
}

template <typename Distances>
void Cross<Distances>::set_ab_cycle(
        const Individual& parent1,
        const Individual& parent2,
        int flags[10],
        int number_of_kids)
{
    number_of_branching_ = 0;
    number_of_unbranched_ = 0;
    for (VertexId vertex_id = 0; vertex_id < number_of_vertices_; ++vertex_id) {
        near_data_[vertex_id][1] = parent1.neighbors[vertex_id][0];
        near_data_[vertex_id][3] = parent1.neighbors[vertex_id][1];
        near_data_[vertex_id][0] = 2;

        unbranched_[number_of_unbranched_] = vertex_id;
        number_of_unbranched_++;

        near_data_[vertex_id][2] = parent2.neighbors[vertex_id][0];
        near_data_[vertex_id][4] = parent2.neighbors[vertex_id][1];
    }
    for (int i = 0; i < number_of_vertices_; ++i) {
        first_visit_position_[i] = -1;
        unbranched_index_[unbranched_[i]] = i;
    }
    number_of_ab_cycles_ = 0;
    start_new_trace_ = 1;
    while (number_of_unbranched_ != 0) {
        if (start_new_trace_ == 1) {
            position_current_ = 0;
            random_pick_ = rand() % number_of_unbranched_;
            trace_start_ = unbranched_[random_pick_];
            first_visit_position_[trace_start_] = position_current_;
            route_[position_current_] = trace_start_;
            current_city_ = trace_start_;
            traversal_type_ = 2;
        } else if (start_new_trace_ == 0) {
            current_city_ = route_[position_current_];
        }

        cycle_complete_ = 0;
        while (cycle_complete_ == 0) {
            position_current_++;
            previous_city_ = current_city_;
            switch (traversal_type_) {
            case 1:
                current_city_ = near_data_[previous_city_][position_current_ % 2 + 1];
                break;
            case 2:
                random_pick_ = rand() % 2;
                current_city_ = near_data_[previous_city_][position_current_ % 2 + 1 + 2 * random_pick_];
                if (random_pick_ == 0)
                    std::swap(near_data_[previous_city_][position_current_ % 2 + 1], near_data_[previous_city_][position_current_ % 2 + 3]);
                break;
            case 3:
                current_city_ = near_data_[previous_city_][position_current_ % 2 + 3];
            }
            route_[position_current_] = current_city_;
            if (near_data_[current_city_][0] == 2) {
                if (current_city_ == trace_start_) {
                    if (first_visit_position_[trace_start_] == 0) {
                        if ((position_current_ - first_visit_position_[trace_start_]) % 2 == 0) {
                            if (near_data_[trace_start_][position_current_ % 2 + 1] == previous_city_)
                                std::swap(near_data_[current_city_][position_current_ % 2 + 1], near_data_[current_city_][position_current_ % 2 + 3]);

                            start_appearance_count_ = 1;
                            this->form_ab_cycle();
                            if (flags[1] == 1 && number_of_ab_cycles_ == number_of_kids)
                                goto RETURN;
                            if (number_of_ab_cycles_ == max_number_of_ab_cycles_)
                                goto RETURN;

                            start_new_trace_ = 0;
                            cycle_complete_ = 1;
                            traversal_type_ = 1;
                        } else {
                            std::swap(near_data_[current_city_][position_current_ % 2 + 1], near_data_[current_city_][position_current_ % 2 + 3]);
                            traversal_type_ = 2;
                        }
                        first_visit_position_[trace_start_] = position_current_;
                    } else {
                        start_appearance_count_ = 2;
                        this->form_ab_cycle();
                        if (flags[1] == 1 && number_of_ab_cycles_ == number_of_kids)
                            goto RETURN;
                        if (number_of_ab_cycles_ == max_number_of_ab_cycles_)
                            goto RETURN;

                        start_new_trace_ = 1;
                        cycle_complete_ = 1;
                    }
                } else if (first_visit_position_[current_city_] == -1) {
                    first_visit_position_[current_city_] = position_current_;
                    if (near_data_[current_city_][position_current_ % 2 + 1] == previous_city_)
                        std::swap(near_data_[current_city_][position_current_ % 2 + 1], near_data_[current_city_][position_current_ % 2 + 3]);
                    traversal_type_ = 2;
                } else if (first_visit_position_[current_city_] > 0) {
                    std::swap(near_data_[current_city_][position_current_ % 2 + 1], near_data_[current_city_][position_current_ % 2 + 3]);
                    if ((position_current_ - first_visit_position_[current_city_]) % 2 == 0) {
                        start_appearance_count_ = 1;
                        this->form_ab_cycle();
                        if (flags[1] == 1 && number_of_ab_cycles_ == number_of_kids)
                            goto RETURN;
                        if (number_of_ab_cycles_ == max_number_of_ab_cycles_)
                            goto RETURN;

                        start_new_trace_ = 0;
                        cycle_complete_ = 1;
                        traversal_type_ = 1;
                    } else {
                        std::swap(near_data_[current_city_][(position_current_ + 1) % 2 + 1], near_data_[current_city_][(position_current_ + 1) % 2 + 3]);
                        traversal_type_ = 3;
                    }
                }
            } else if (near_data_[current_city_][0] == 1) {
                if (current_city_ == trace_start_) {
                    start_appearance_count_ = 1;
                    this->form_ab_cycle();
                    if (flags[1] == 1 && number_of_ab_cycles_ == number_of_kids)
                        goto RETURN;
                    if (number_of_ab_cycles_ == max_number_of_ab_cycles_)
                        goto RETURN;
                    start_new_trace_ = 1;
                    cycle_complete_ = 1;
                } else {
                    traversal_type_ = 1;
                }
            }
        }
    }
    while (number_of_branching_ != 0) {
        position_current_ = 0;
        random_pick_ = rand() % number_of_branching_;
        trace_start_ = branching_[random_pick_];
        route_[position_current_] = trace_start_;
        current_city_ = trace_start_;

        cycle_complete_ = 0;
        while (cycle_complete_ == 0) {
            previous_city_ = current_city_;
            position_current_++;
            current_city_ = near_data_[previous_city_][position_current_ % 2 + 1];
            route_[position_current_] = current_city_;
            if (current_city_ == trace_start_) {
                start_appearance_count_ = 1;
                this->form_ab_cycle();
                if (flags[1] == 1 && number_of_ab_cycles_ == number_of_kids)
                    goto RETURN;
                if (number_of_ab_cycles_ == max_number_of_ab_cycles_)
                    goto RETURN;

                cycle_complete_ = 1;
            }
        }
    }
RETURN:
    if (number_of_ab_cycles_ == max_number_of_ab_cycles_) {
        printf("max_number_of_ab_cycles_(%d) must be increased\n", max_number_of_ab_cycles_);
        exit(1);
    }
}

template <typename Distances>
void Cross<Distances>::form_ab_cycle()
{
    VertexId cycle_start;
    VertexId visiting_city;
    VertexId stock;
    int start_count;
    int edge_type;
    int cycle_length;
    Distance diff;

    edge_type = (position_current_ % 2 == 0)? 1: 2;
    cycle_start = route_[position_current_];
    cycle_length = 0;
    c_[cycle_length] = cycle_start;

    start_count = 0;
    while (true) {
        ++cycle_length;
        --position_current_;
        visiting_city = route_[position_current_];
        if (near_data_[visiting_city][0] == 2) {
            unbranched_[unbranched_index_[visiting_city]] = unbranched_[number_of_unbranched_ - 1];
            unbranched_index_[unbranched_[number_of_unbranched_ - 1]] = unbranched_index_[visiting_city];
            --number_of_unbranched_;
            branching_[number_of_branching_] = visiting_city;
            branching_index_[visiting_city] = number_of_branching_;
            ++number_of_branching_;
        } else if (near_data_[visiting_city][0] == 1) {
            branching_[branching_index_[visiting_city]] = branching_[number_of_branching_ - 1];
            branching_index_[branching_[number_of_branching_ - 1]] = branching_index_[visiting_city];
            --number_of_branching_;
        }

        --near_data_[visiting_city][0];
        if (visiting_city == cycle_start)
            ++start_count;
        if (start_count == start_appearance_count_)
            break;
        c_[cycle_length] = visiting_city;
    }

    if (cycle_length == 2)
        return;

    ab_cycle_[number_of_ab_cycles_][0] = cycle_length;

    if (edge_type == 2) {
        stock = c_[0];
        for (int j = 0; j < cycle_length - 1; ++j)
            c_[j] = c_[j + 1];
        c_[cycle_length - 1] = stock;
    }

    for (int j = 0; j < cycle_length; ++j)
        ab_cycle_[number_of_ab_cycles_][j + 2] = c_[j];

    ab_cycle_[number_of_ab_cycles_][1] = c_[cycle_length - 1];
    ab_cycle_[number_of_ab_cycles_][cycle_length + 2] = c_[0];
    ab_cycle_[number_of_ab_cycles_][cycle_length + 3] = c_[1];

    c_[cycle_length] = c_[0];
    c_[cycle_length + 1] = c_[1];
    diff = 0;
    for (int j = 0; j < cycle_length / 2; ++j)
        diff = diff + evaluator_.distance(c_[2 * j], c_[1 + 2 * j]) - evaluator_.distance(c_[1 + 2 * j], c_[2 + 2 * j]);

    gain_ab_[number_of_ab_cycles_] = diff;
    ++number_of_ab_cycles_;
}

template <typename Distances>
void Cross<Distances>::change_sol(
        Individual& child,
        int ab_number,
        int type)
{
    int j;
    int cycle_length;
    VertexId red_vertex_id_1, red_vertex_id_2, blue_vertex_id_1, blue_vertex_id_2;
    int red_position_1, red_position_2, blue_position_1, blue_position_2;

    cycle_length = ab_cycle_[ab_number][0];
    c_[0] = ab_cycle_[ab_number][0];

    if (type == 2) {
        for (j = 0; j < cycle_length + 3; ++j)
            c_[cycle_length + 3 - j] = ab_cycle_[ab_number][j + 1];
    } else {
        for (j = 1; j <= cycle_length + 3; ++j)
            c_[j] = ab_cycle_[ab_number][j];
    }

    for (j = 0; j < cycle_length / 2; ++j) {
        red_vertex_id_1 = c_[2 + 2 * j];
        red_vertex_id_2 = c_[3 + 2 * j];
        blue_vertex_id_1 = c_[1 + 2 * j];
        blue_vertex_id_2 = c_[4 + 2 * j];

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

        red_position_1 = inverse_order_[red_vertex_id_1];
        red_position_2 = inverse_order_[red_vertex_id_2];
        blue_position_1 = inverse_order_[blue_vertex_id_1];
        blue_position_2 = inverse_order_[blue_vertex_id_2];

        if (red_position_1 == 0 && red_position_2 == number_of_vertices_ - 1) {
            segment_position_list_[number_of_segment_positions_++] = red_position_1;
        } else if (red_position_1 == number_of_vertices_ - 1 && red_position_2 == 0) {
            segment_position_list_[number_of_segment_positions_++] = red_position_2;
        } else if (red_position_1 < red_position_2) {
            segment_position_list_[number_of_segment_positions_++] = red_position_2;
        } else if (red_position_2 < red_position_1) {
            segment_position_list_[number_of_segment_positions_++] = red_position_1;
        }

        link_b_position_[red_position_1][1] = link_b_position_[red_position_1][0];
        link_b_position_[red_position_2][1] = link_b_position_[red_position_2][0];
        link_b_position_[red_position_1][0] = blue_position_1;
        link_b_position_[red_position_2][0] = blue_position_2;
    }
}

template <typename Distances>
void Cross<Distances>::make_complete_sol(Individual& child)
{
    int j1, j2;
    VertexId unit_start, previous_vertex_id, current_vertex_id, next_vertex_id;
    VertexId vertex_id_1, vertex_id_2, vertex_id_3, vertex_id_4;
    VertexId best_vertex_id_1, best_vertex_id_2, best_vertex_id_3, best_vertex_id_4;
    int min_unit_city;
    int center_unit_index, selected_unit_index;
    Distance diff, max_diff;
    int near_num, near_search_limit;

    modification_gain_ = 0;
    while (number_of_units_ != 1) {
        min_unit_city = number_of_vertices_ + 12345;
        for (int u = 0; u < number_of_units_; ++u)
            if (number_of_elements_in_unit_[u] < min_unit_city) {
                center_unit_index = u;
                min_unit_city = number_of_elements_in_unit_[u];
            }

        unit_start = -1;
        number_of_segments_for_center_ = 0;
        for (int s = 0; s < number_of_segments_; ++s)
            if (segment_unit_[s] == center_unit_index) {
                int posi = segment_[s][0];
                unit_start = order_[posi];
                segment_for_center_[number_of_segments_for_center_++] = s;
            }
        current_vertex_id = -1;
        next_vertex_id = unit_start;
        number_of_elements_in_center_unit_ = 0;
        while (true) {
            previous_vertex_id = current_vertex_id;
            current_vertex_id = next_vertex_id;
            center_unit_[current_vertex_id] = 1;
            list_of_center_unit_[number_of_elements_in_center_unit_] = current_vertex_id;
            ++number_of_elements_in_center_unit_;
            if (child.neighbors[current_vertex_id][0] != previous_vertex_id) {
                next_vertex_id = child.neighbors[current_vertex_id][0];
            } else {
                next_vertex_id = child.neighbors[current_vertex_id][1];
            }
            if (next_vertex_id == unit_start)
                break;
        }
        list_of_center_unit_[number_of_elements_in_center_unit_] = list_of_center_unit_[0];
        list_of_center_unit_[number_of_elements_in_center_unit_ + 1] = list_of_center_unit_[1];

        max_diff = std::numeric_limits<Distance>::min();
        best_vertex_id_3 = -1;
        best_vertex_id_4 = -1;
        near_search_limit = 10;   // N_near
        // near_search_limit <= eva->fNearNumMax (kopt.cpp)

    RESTART:
        for (int s = 1; s <= number_of_elements_in_center_unit_; ++s) {
            vertex_id_1 = list_of_center_unit_[s];

            for (near_num = 1; near_num <= near_search_limit; ++near_num) {
                vertex_id_3 = evaluator_.near_cities[vertex_id_1][near_num];
                if (center_unit_[vertex_id_3] == 0) {
                    for (j1 = 0; j1 < 2; ++j1) {
                        vertex_id_2 = list_of_center_unit_[s - 1 + 2 * j1];
                        for (j2 = 0; j2 < 2; ++j2) {
                            vertex_id_4 = child.neighbors[vertex_id_3][j2];
                            diff = evaluator_.distance(vertex_id_1, vertex_id_2) + evaluator_.distance(vertex_id_3, vertex_id_4) - evaluator_.distance(vertex_id_1, vertex_id_3) - evaluator_.distance(vertex_id_2, vertex_id_4);
                            if (diff > max_diff) {
                                best_vertex_id_1 = vertex_id_1;
                                best_vertex_id_2 = vertex_id_2;
                                best_vertex_id_3 = vertex_id_3;
                                best_vertex_id_4 = vertex_id_4;
                                max_diff = diff;
                            }
                            diff = evaluator_.distance(vertex_id_1, vertex_id_2) + evaluator_.distance(vertex_id_4, vertex_id_3) -
                                evaluator_.distance(vertex_id_1, vertex_id_4) - evaluator_.distance(vertex_id_2, vertex_id_3);
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
            int random_center_index = rand() % (number_of_elements_in_center_unit_ - 1);
            vertex_id_1 = list_of_center_unit_[random_center_index];
            vertex_id_2 = list_of_center_unit_[random_center_index + 1];
            for (VertexId vertex_id = 0; vertex_id < number_of_vertices_; ++vertex_id) {
                if (center_unit_[vertex_id] == 0) {
                    best_vertex_id_1 = vertex_id_1;
                    best_vertex_id_2 = vertex_id_2;
                    best_vertex_id_3 = vertex_id;
                    best_vertex_id_4 = child.neighbors[vertex_id][0];
                    break;
                }
            }
            max_diff = evaluator_.distance(best_vertex_id_1, best_vertex_id_2) + evaluator_.distance(best_vertex_id_3, best_vertex_id_4) - evaluator_.distance(vertex_id_1, best_vertex_id_3) - evaluator_.distance(vertex_id_2, best_vertex_id_4);
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

        modified_edge_[number_of_modified_edges_][0] = best_vertex_id_1;
        modified_edge_[number_of_modified_edges_][1] = best_vertex_id_2;
        modified_edge_[number_of_modified_edges_][2] = best_vertex_id_3;
        modified_edge_[number_of_modified_edges_][3] = best_vertex_id_4;
        ++number_of_modified_edges_;

        modification_gain_ += max_diff;

        int best_position_3 = inverse_order_[best_vertex_id_3];
        selected_unit_index = -1;
        for (int s = 0; s < number_of_segments_; ++s)
            if (segment_[s][0] <= best_position_3 && best_position_3 <= segment_[s][1]) {
                selected_unit_index = segment_unit_[s];
                break;
            }

        for (int s = 0; s < number_of_segments_; ++s)
            if (segment_unit_[s] == selected_unit_index)
                segment_unit_[s] = center_unit_index;

        number_of_elements_in_unit_[center_unit_index] += number_of_elements_in_unit_[selected_unit_index];

        for (int s = 0; s < number_of_segments_; ++s)
            if (segment_unit_[s] == number_of_units_ - 1)
                segment_unit_[s] = selected_unit_index;

        number_of_elements_in_unit_[selected_unit_index] = number_of_elements_in_unit_[number_of_units_ - 1];
        --number_of_units_;

        for (int s = 0; s < number_of_elements_in_center_unit_; ++s) {
            VertexId vertex_id = list_of_center_unit_[s];
            center_unit_[vertex_id] = 0;
        }
    }
}

template <typename Distances>
void Cross<Distances>::make_unit()
{
    int flag = 1;
    for (int s = 0; s < number_of_segment_positions_; ++s) {
        if (segment_position_list_[s] == 0) {
            flag = 0;
            break;
        }
    }
    if (flag == 1) {
        segment_position_list_[number_of_segment_positions_++] = 0;
        link_b_position_[number_of_vertices_ - 1][1]  = link_b_position_[number_of_vertices_ - 1][0];
        link_b_position_[0][1] = link_b_position_[0][0];
        link_b_position_[number_of_vertices_ - 1][0] = 0;
        link_b_position_[0][0] = number_of_vertices_ - 1;
    }

    sort_ascending(segment_position_list_, number_of_segment_positions_);
    number_of_segments_ = number_of_segment_positions_;
    for (int s = 0; s < number_of_segments_ - 1; ++s) {
        segment_[s][0] = segment_position_list_[s];
        segment_[s][1] = segment_position_list_[s + 1] - 1;
    }

    segment_[number_of_segments_ - 1][0] = segment_position_list_[number_of_segments_ - 1];
    segment_[number_of_segments_ - 1][1] = number_of_vertices_ - 1;

    for (int s = 0; s < number_of_segments_; ++s) {
        link_a_position_[segment_[s][0]] = segment_[s][1];
        link_a_position_[segment_[s][1]] = segment_[s][0];
        position_segment_[segment_[s][0]] = s;
        position_segment_[segment_[s][1]] = s;
    }

    for (int s = 0; s < number_of_segments_; ++s)
        segment_unit_[s] = -1;
    number_of_units_ = 0;

    int start_position, position1, position2, next_position, previous_position;
    int segment_number;
    while (1) {
        flag = 0;
        for (int s = 0; s < number_of_segments_; ++s) {
            if (segment_unit_[s] == -1) {
                start_position = segment_[s][0];
                previous_position = -1;
                position1 = start_position;
                flag = 1;
                break;
            }
        }
        if (flag == 0)
            break;

        while (1) {
            segment_number = position_segment_[position1];
            segment_unit_[segment_number] = number_of_units_;

            position2 = link_a_position_[position1];
            next_position = link_b_position_[position2][0];
            if (position1 == position2)
                if (next_position == previous_position)
                    next_position = link_b_position_[position2][1];

            if (next_position == start_position) {
                ++number_of_units_;
                break;
            }

            previous_position = position2;
            position1 = next_position;
        }
    }

    for (int s = 0; s < number_of_units_; ++s)
        number_of_elements_in_unit_[s] = 0;

    int unit_number = -1;
    int merged_segment_count = -1;
    for (int s = 0; s < number_of_segments_; ++s) {
        if (segment_unit_[s] != unit_number) {
            ++merged_segment_count;
            segment_[merged_segment_count][0] = segment_[s][0];
            segment_[merged_segment_count][1] = segment_[s][1];
            unit_number = segment_unit_[s];
            segment_unit_[merged_segment_count] = unit_number;
            number_of_elements_in_unit_[unit_number] += segment_[s][1] - segment_[s][0] + 1;
        } else {
            segment_[merged_segment_count][1] = segment_[s][1];
            number_of_elements_in_unit_[unit_number] += segment_[s][1] - segment_[s][0] + 1;
        }
    }
    number_of_segments_ = merged_segment_count + 1;
}

template <typename Distances>
void Cross<Distances>::back_to_pa1(Individual& child)
{
    VertexId vertex_id_1, vertex_id_2, vertex_id_3, vertex_id_4;
    int jnum;

    for (int s = number_of_modified_edges_ - 1; s >= 0; --s) {
        vertex_id_1 = modified_edge_[s][0];
        vertex_id_3 = modified_edge_[s][1];
        vertex_id_2 = modified_edge_[s][2];
        vertex_id_4 = modified_edge_[s][3];

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

    for (int s = 0; s < number_of_applied_cycles_; ++s) {
        jnum = applied_cycle_[s];
        this->change_sol(child, jnum, 2);
    }
}

template <typename Distances>
void Cross<Distances>::go_to_best(Individual& child)
{
    VertexId vertex_id_1, vertex_id_2, vertex_id_3, vertex_id_4;
    int jnum;

    for (int s = 0; s < number_of_best_applied_cycles_; ++s) {
        jnum = best_applied_cycle_[s];
        this->change_sol(child, jnum, 1);
    }

    for (int s = 0; s < number_of_best_modified_edges_; ++s) {
        vertex_id_1 = best_modified_edge_[s][0];
        vertex_id_2 = best_modified_edge_[s][1];
        vertex_id_3 = best_modified_edge_[s][2];
        vertex_id_4 = best_modified_edge_[s][3];

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

template <typename Distances>
void Cross<Distances>::increment_edge_freq(
        std::vector<std::vector<int>>& edge_frequency)
{
    int j, jnum, cycle_length;
    VertexId red_vertex_id_1, red_vertex_id_2, blue_vertex_id_1, blue_vertex_id_2;
    VertexId vertex_id_1, vertex_id_2, vertex_id_3, vertex_id_4;

    for (int s = 0; s < number_of_best_applied_cycles_; ++s) {
        jnum = best_applied_cycle_[s];

        cycle_length = ab_cycle_[jnum][0];
        c_[0] = ab_cycle_[jnum][0];

        for (j = 1; j <= cycle_length + 3; ++j)
            c_[j] = ab_cycle_[jnum][j];

        for (j = 0; j < cycle_length / 2; ++j) {
            red_vertex_id_1 = c_[2 + 2 * j];
            red_vertex_id_2 = c_[3 + 2 * j];
            blue_vertex_id_1 = c_[1 + 2 * j];
            blue_vertex_id_2 = c_[4 + 2 * j];

            ++edge_frequency[red_vertex_id_1][blue_vertex_id_1];
            --edge_frequency[red_vertex_id_1][red_vertex_id_2];
            --edge_frequency[red_vertex_id_2][red_vertex_id_1];
            ++edge_frequency[red_vertex_id_2][blue_vertex_id_2];
        }
    }
    for (int s = 0; s < number_of_best_modified_edges_; ++s) {
        vertex_id_1 = best_modified_edge_[s][0];
        vertex_id_2 = best_modified_edge_[s][1];
        vertex_id_3 = best_modified_edge_[s][2];
        vertex_id_4 = best_modified_edge_[s][3];

        --edge_frequency[vertex_id_1][vertex_id_2];
        --edge_frequency[vertex_id_3][vertex_id_4];
        ++edge_frequency[vertex_id_1][vertex_id_3];
        ++edge_frequency[vertex_id_2][vertex_id_4];
        --edge_frequency[vertex_id_2][vertex_id_1];
        --edge_frequency[vertex_id_4][vertex_id_3];
        ++edge_frequency[vertex_id_3][vertex_id_1];
        ++edge_frequency[vertex_id_4][vertex_id_2];
    }
}

template <typename Distances>
int Cross<Distances>::calc_adaptive_loss(
        std::vector<std::vector<int>>& edge_frequency)
{
    int j, jnum, cycle_length;
    VertexId red_vertex_id_1, red_vertex_id_2, blue_vertex_id_1, blue_vertex_id_2;
    VertexId vertex_id_1, vertex_id_2, vertex_id_3, vertex_id_4;
    double loss;

    loss = 0;
    for (int s = 0; s < number_of_applied_cycles_; ++s) {
        jnum = applied_cycle_[s];

        cycle_length = ab_cycle_[jnum][0];
        c_[0] = ab_cycle_[jnum][0];

        for (j = 1; j <= cycle_length + 3; ++j)
            c_[j] = ab_cycle_[jnum][j];

        for (j = 0; j < cycle_length / 2; ++j) {
            red_vertex_id_1 = c_[2 + 2 * j];
            red_vertex_id_2 = c_[3 + 2 * j];
            blue_vertex_id_1 = c_[1 + 2 * j];
            blue_vertex_id_2 = c_[4 + 2 * j];

            loss -= (edge_frequency[red_vertex_id_1][red_vertex_id_2] - 1);
            loss -= (edge_frequency[red_vertex_id_2][red_vertex_id_1] - 1);
            loss += edge_frequency[red_vertex_id_2][blue_vertex_id_2];
            loss += edge_frequency[blue_vertex_id_2][red_vertex_id_2];

            --edge_frequency[red_vertex_id_1][red_vertex_id_2];
            --edge_frequency[red_vertex_id_2][red_vertex_id_1];
            ++edge_frequency[red_vertex_id_2][blue_vertex_id_2];
            ++edge_frequency[blue_vertex_id_2][red_vertex_id_2];
        }
    }
    for (int s = 0; s < number_of_modified_edges_; ++s) {
        vertex_id_1 = modified_edge_[s][0];
        vertex_id_2 = modified_edge_[s][1];
        vertex_id_3 = modified_edge_[s][2];
        vertex_id_4 = modified_edge_[s][3];

        loss -= (edge_frequency[vertex_id_1][vertex_id_2] - 1);
        loss -= (edge_frequency[vertex_id_2][vertex_id_1] - 1);
        loss -= (edge_frequency[vertex_id_3][vertex_id_4] - 1);
        loss -= (edge_frequency[vertex_id_4][vertex_id_3] - 1);

        loss += edge_frequency[vertex_id_1][vertex_id_3];
        loss += edge_frequency[vertex_id_3][vertex_id_1];
        loss += edge_frequency[vertex_id_2][vertex_id_4];
        loss += edge_frequency[vertex_id_4][vertex_id_2];

        --edge_frequency[vertex_id_1][vertex_id_2];
        --edge_frequency[vertex_id_2][vertex_id_1];
        --edge_frequency[vertex_id_3][vertex_id_4];
        --edge_frequency[vertex_id_4][vertex_id_3];

        ++edge_frequency[vertex_id_1][vertex_id_3];
        ++edge_frequency[vertex_id_3][vertex_id_1];
        ++edge_frequency[vertex_id_2][vertex_id_4];
        ++edge_frequency[vertex_id_4][vertex_id_2];
    }
    for (int s = 0; s < number_of_applied_cycles_; ++s) {
        jnum = applied_cycle_[s];
        cycle_length = ab_cycle_[jnum][0];
        c_[0] = ab_cycle_[jnum][0];
        for (j = 1; j <= cycle_length + 3; ++j)
            c_[j] = ab_cycle_[jnum][j];

        for (j = 0; j < cycle_length / 2; ++j) {
            red_vertex_id_1 = c_[2 + 2 * j];
            red_vertex_id_2 = c_[3 + 2 * j];
            blue_vertex_id_1 = c_[1 + 2 * j];
            blue_vertex_id_2 = c_[4 + 2 * j];

            ++edge_frequency[red_vertex_id_1][red_vertex_id_2];
            ++edge_frequency[red_vertex_id_2][red_vertex_id_1];
            --edge_frequency[red_vertex_id_2][blue_vertex_id_2];
            --edge_frequency[blue_vertex_id_2][red_vertex_id_2];
        }
    }
    for (int s = 0; s < number_of_modified_edges_; ++s) {
        vertex_id_1 = modified_edge_[s][0];
        vertex_id_2 = modified_edge_[s][1];
        vertex_id_3 = modified_edge_[s][2];
        vertex_id_4 = modified_edge_[s][3];

        ++edge_frequency[vertex_id_1][vertex_id_2];
        ++edge_frequency[vertex_id_2][vertex_id_1];
        ++edge_frequency[vertex_id_3][vertex_id_4];
        ++edge_frequency[vertex_id_4][vertex_id_3];

        --edge_frequency[vertex_id_1][vertex_id_3];
        --edge_frequency[vertex_id_3][vertex_id_1];
        --edge_frequency[vertex_id_2][vertex_id_4];
        --edge_frequency[vertex_id_4][vertex_id_2];
    }
    return int(loss / 2);
}

template <typename Distances>
double Cross<Distances>::calc_entropy_loss(
        std::vector<std::vector<int>>& edge_frequency)
{
    int j, jnum, cycle_length;
    VertexId red_vertex_id_1, red_vertex_id_2, blue_vertex_id_1, blue_vertex_id_2;
    VertexId vertex_id_1, vertex_id_2, vertex_id_3, vertex_id_4;
    double loss;
    double h1, h2;

    loss = 0;  // AB-cycle
    for (int s = 0; s < number_of_applied_cycles_; ++s) {
        jnum = applied_cycle_[s];
        cycle_length = ab_cycle_[jnum][0];
        c_[0] = ab_cycle_[jnum][0];

        for (j = 1; j <= cycle_length + 3; ++j)
            c_[j] = ab_cycle_[jnum][j];

        for (j = 0; j < cycle_length / 2; ++j) {
            red_vertex_id_1 = c_[2 + 2 * j];
            red_vertex_id_2 = c_[3 + 2 * j];
            blue_vertex_id_1 = c_[1 + 2 * j];
            blue_vertex_id_2 = c_[4 + 2 * j];

            h1 = (double)(edge_frequency[red_vertex_id_1][red_vertex_id_2] - 1) / (double)population_size_;
            h2 = (double)(edge_frequency[red_vertex_id_1][red_vertex_id_2]) / (double)population_size_;
            if (edge_frequency[red_vertex_id_1][red_vertex_id_2] - 1 != 0)
                loss -= h1 * log(h1);
            loss += h2 * log(h2);
            --edge_frequency[red_vertex_id_1][red_vertex_id_2];
            --edge_frequency[red_vertex_id_2][red_vertex_id_1];

            h1 = (double)(edge_frequency[red_vertex_id_2][blue_vertex_id_2] + 1) / (double)population_size_;
            h2 = (double)(edge_frequency[red_vertex_id_2][blue_vertex_id_2]) / (double)population_size_;
            loss -= h1 * log(h1);
            if (edge_frequency[red_vertex_id_2][blue_vertex_id_2] != 0)
                loss += h2 * log(h2);
            ++edge_frequency[red_vertex_id_2][blue_vertex_id_2];
            ++edge_frequency[blue_vertex_id_2][red_vertex_id_2];
        }
    }

    for (int s = 0; s < number_of_modified_edges_; ++s) {
        vertex_id_1 = modified_edge_[s][0];
        vertex_id_2 = modified_edge_[s][1];
        vertex_id_3 = modified_edge_[s][2];
        vertex_id_4 = modified_edge_[s][3];

        h1 = (double)(edge_frequency[vertex_id_1][vertex_id_2] - 1) / (double)population_size_;
        h2 = (double)(edge_frequency[vertex_id_1][vertex_id_2]) / (double)population_size_;
        if (edge_frequency[vertex_id_1][vertex_id_2] - 1 != 0)
            loss -= h1 * log(h1);
        loss += h2 * log(h2);
        --edge_frequency[vertex_id_1][vertex_id_2];
        --edge_frequency[vertex_id_2][vertex_id_1];

        h1 = (double)(edge_frequency[vertex_id_3][vertex_id_4] - 1) / (double)population_size_;
        h2 = (double)(edge_frequency[vertex_id_3][vertex_id_4]) / (double)population_size_;
        if (edge_frequency[vertex_id_3][vertex_id_4] - 1 != 0)
            loss -= h1 * log(h1);
        loss += h2 * log(h2);
        --edge_frequency[vertex_id_3][vertex_id_4];
        --edge_frequency[vertex_id_4][vertex_id_3];

        h1 = (double)(edge_frequency[vertex_id_1][vertex_id_3] + 1) / (double)population_size_;
        h2 = (double)(edge_frequency[vertex_id_1][vertex_id_3]) / (double)population_size_;
        loss -= h1 * log(h1);
        if (edge_frequency[vertex_id_1][vertex_id_3] != 0)
            loss += h2 * log(h2);
        ++edge_frequency[vertex_id_1][vertex_id_3];
        ++edge_frequency[vertex_id_3][vertex_id_1];

        h1 = (double)(edge_frequency[vertex_id_2][vertex_id_4] + 1) / (double)population_size_;
        h2 = (double)(edge_frequency[vertex_id_2][vertex_id_4]) / (double)population_size_;
        loss -= h1 * log(h1);
        if (edge_frequency[vertex_id_2][vertex_id_4] != 0)
            loss += h2 * log(h2);
        ++edge_frequency[vertex_id_2][vertex_id_4];
        ++edge_frequency[vertex_id_4][vertex_id_2];
    }
    loss = -loss;

    // restores edge_frequency
    for (int s = 0; s < number_of_applied_cycles_; ++s) {
        jnum = applied_cycle_[s];

        cycle_length = ab_cycle_[jnum][0];
        c_[0] = ab_cycle_[jnum][0];

        for (j = 1; j <= cycle_length + 3; ++j)
            c_[j] = ab_cycle_[jnum][j];

        for (j = 0; j < cycle_length / 2; ++j) {
            red_vertex_id_1 = c_[2 + 2 * j];
            red_vertex_id_2 = c_[3 + 2 * j];
            blue_vertex_id_1 = c_[1 + 2 * j];
            blue_vertex_id_2 = c_[4 + 2 * j];

            ++edge_frequency[red_vertex_id_1][red_vertex_id_2];
            ++edge_frequency[red_vertex_id_2][red_vertex_id_1];
            --edge_frequency[red_vertex_id_2][blue_vertex_id_2];
            --edge_frequency[blue_vertex_id_2][red_vertex_id_2];
        }
    }
    for (int s = 0; s < number_of_modified_edges_; ++s) {
        vertex_id_1 = modified_edge_[s][0];
        vertex_id_2 = modified_edge_[s][1];
        vertex_id_3 = modified_edge_[s][2];
        vertex_id_4 = modified_edge_[s][3];

        ++edge_frequency[vertex_id_1][vertex_id_2];
        ++edge_frequency[vertex_id_2][vertex_id_1];
        ++edge_frequency[vertex_id_3][vertex_id_4];
        ++edge_frequency[vertex_id_4][vertex_id_3];

        --edge_frequency[vertex_id_1][vertex_id_3];
        --edge_frequency[vertex_id_3][vertex_id_1];
        --edge_frequency[vertex_id_2][vertex_id_4];
        --edge_frequency[vertex_id_4][vertex_id_2];
    }
    return loss;
}

template <typename Distances>
void Cross<Distances>::set_weight(
        const Individual& parent1,
        const Individual& parent2)
{
    int cycle_length;
    VertexId red_vertex_id_1, red_vertex_id_2, current_vertex_id, next_vertex_id, previous_vertex_id;
    int ab_number;

    for (VertexId vertex_id = 0; vertex_id < number_of_vertices_; ++vertex_id) {
        in_effect_node_[vertex_id][0] = -1;
        in_effect_node_[vertex_id][1] = -1;
    }

    // Step 1:
    for (int s = 0; s < number_of_ab_cycles_; ++s) {
        cycle_length = ab_cycle_[s][0];
        for (int j = 0; j < cycle_length / 2; ++j) {
            red_vertex_id_1 = ab_cycle_[s][2 * j + 2]; // red edge
            red_vertex_id_2 = ab_cycle_[s][2 * j + 3];

            if (in_effect_node_[red_vertex_id_1][0] == -1) {
                in_effect_node_[red_vertex_id_1][0] = s;
            } else if (in_effect_node_[red_vertex_id_1][1] == -1) {
                in_effect_node_[red_vertex_id_1][1] = s;
            }

            if (in_effect_node_[red_vertex_id_2][0] == -1) {
                in_effect_node_[red_vertex_id_2][0] = s;
            } else if (in_effect_node_[red_vertex_id_2][1] == -1) {
                in_effect_node_[red_vertex_id_2][1] = s;
            }
        }
    }

    // Step 2:
    for (VertexId vertex_id = 0; vertex_id < number_of_vertices_; ++vertex_id) {
        if (in_effect_node_[vertex_id][0] != -1 && in_effect_node_[vertex_id][1] == -1) {
            ab_number = in_effect_node_[vertex_id][0];
            current_vertex_id = vertex_id;

            if (parent1.neighbors[current_vertex_id][0] != parent2.neighbors[current_vertex_id][0] && parent1.neighbors[current_vertex_id][0] != parent2.neighbors[current_vertex_id][1]) {
                previous_vertex_id = parent1.neighbors[current_vertex_id][0];
            } else if (parent1.neighbors[current_vertex_id][1] != parent2.neighbors[current_vertex_id][0] && parent1.neighbors[current_vertex_id][1] != parent2.neighbors[current_vertex_id][1]) {
                previous_vertex_id = parent1.neighbors[current_vertex_id][1];
            }

            while (true) {
                in_effect_node_[current_vertex_id][1] = ab_number;

                if (parent1.neighbors[current_vertex_id][0] != previous_vertex_id) {
                    next_vertex_id = parent1.neighbors[current_vertex_id][0];
                } else if (parent1.neighbors[current_vertex_id][1] != previous_vertex_id) {
                    next_vertex_id = parent1.neighbors[current_vertex_id][1];
                }

                if (in_effect_node_[next_vertex_id][0] == -1) {
                    in_effect_node_[next_vertex_id][0] = ab_number;
                } else if (in_effect_node_[next_vertex_id][1] == -1) {
                    in_effect_node_[next_vertex_id][1] = ab_number;
                }

                if (in_effect_node_[next_vertex_id][1] != -1)
                    break;
                previous_vertex_id = current_vertex_id;
                current_vertex_id = next_vertex_id;
            }
        }
    }

    // Step 3:

    for (int s1 = 0; s1 < number_of_ab_cycles_; ++s1) {
        weight_c_[s1] = 0;
        for (int s2 = 0; s2 < number_of_ab_cycles_; ++s2)
            weight_rr_[s1][s2] = 0;
    }

    for (VertexId vertex_id = 0; vertex_id < number_of_vertices_; ++vertex_id) {
        if (in_effect_node_[vertex_id][0] != -1 && in_effect_node_[vertex_id][1] != -1) {
            ++weight_rr_[in_effect_node_[vertex_id][0]][in_effect_node_[vertex_id][1]];
            ++weight_rr_[in_effect_node_[vertex_id][1]][in_effect_node_[vertex_id][0]];
        }
        if (in_effect_node_[vertex_id][0] != in_effect_node_[vertex_id][1]) {
            ++weight_c_[in_effect_node_[vertex_id][0]];
            ++weight_c_[in_effect_node_[vertex_id][1]];
        }
    }
    for (int s1 = 0; s1 < number_of_ab_cycles_; ++s1)
        weight_rr_[s1][s1] = 0;
}

template <typename Distances>
int Cross<Distances>::calc_c_naive()
{
    int count_c_nodes;
    int tally;

    count_c_nodes = 0;

    for (VertexId vertex_id = 0; vertex_id < number_of_vertices_; ++vertex_id) {
        if (in_effect_node_[vertex_id][0] != -1 && in_effect_node_[vertex_id][1] != -1) {
            tally = 0;
            if (used_ab_cycle_[in_effect_node_[vertex_id][0]] == 1)
                ++tally;
            if (used_ab_cycle_[in_effect_node_[vertex_id][1]] == 1)
                ++tally;
            if (tally == 1)
                ++count_c_nodes;
        }
    }
    return count_c_nodes;
}

template <typename Distances>
void Cross<Distances>::search_eset(int center_ab)
{
    int iteration, stagnation;
    int delta_weight, min_delta_weight_non_tabu;
    int improving_change, non_tabu_change;
    int selected_ab_cycle, selected_ab_cycle_non_tabu;
    int jnum;

    number_of_c_nodes_ = 0; // Number of C nodes in E-set
    number_of_e_edges_ = 0; // Number of Edges in E-set

    number_of_used_ab_cycles_ = 0;
    for (int s1 = 0; s1 < number_of_ab_cycles_; ++s1) {
        used_ab_cycle_[s1] = 0;
        weight_sr_[s1] = 0;
        moved_ab_cycle_[s1] = 0;
    }

    for (int s = 0; s < number_of_ab_cycles_in_eset_; ++s) {
        jnum = ab_cycle_in_eset_[s];
        this->add_ab(jnum);
    }
    best_number_of_c_nodes_ = number_of_c_nodes_;
    best_number_of_e_edges_ = number_of_e_edges_;

    stagnation = 0;
    iteration = 0;
    while (1) {
        ++iteration;
        min_delta_weight_non_tabu = 99999999;
        improving_change = 0;
        non_tabu_change = 0;
        for (int s1 = 0; s1 < number_of_ab_cycles_; ++s1) {
            if (used_ab_cycle_[s1] == 0 && weight_sr_[s1] > 0) {
                delta_weight = weight_c_[s1] - 2 * weight_sr_[s1];
                if (number_of_c_nodes_ + delta_weight < best_number_of_c_nodes_) {
                    selected_ab_cycle = s1;
                    improving_change = 1;
                    best_number_of_c_nodes_ = number_of_c_nodes_ + delta_weight;
                }
                if (delta_weight < min_delta_weight_non_tabu && iteration > moved_ab_cycle_[s1]) {
                    selected_ab_cycle_non_tabu = s1;
                    non_tabu_change = 1;
                    min_delta_weight_non_tabu = delta_weight;
                }
            } else if (used_ab_cycle_[s1] == 1 && s1 != center_ab) {
                delta_weight = - weight_c_[s1] + 2 * weight_sr_[s1];
                if (number_of_c_nodes_ + delta_weight < best_number_of_c_nodes_) {
                    selected_ab_cycle = s1;
                    improving_change = -1;
                    best_number_of_c_nodes_ = number_of_c_nodes_ + delta_weight;
                }
                if (delta_weight < min_delta_weight_non_tabu && iteration > moved_ab_cycle_[s1]) {
                    selected_ab_cycle_non_tabu = s1;
                    non_tabu_change = -1;
                    min_delta_weight_non_tabu = delta_weight;
                }
            }
        }

        if (improving_change != 0) {
            if (improving_change == 1) {
                this->add_ab(selected_ab_cycle);
            } else if (improving_change == -1) {
                this->delete_ab(selected_ab_cycle);
            }

            moved_ab_cycle_[selected_ab_cycle] = iteration + random_integer(1, t_max_);

            best_number_of_e_edges_ = number_of_e_edges_;

            number_of_ab_cycles_in_eset_ = 0;
            for (int s1 = 0; s1 < number_of_ab_cycles_; ++s1)
                if (used_ab_cycle_[s1] == 1)
                    ab_cycle_in_eset_[number_of_ab_cycles_in_eset_++] = s1;

            stagnation = 0;
        } else if (non_tabu_change != 0) {
            if (non_tabu_change == 1) {
                this->add_ab(selected_ab_cycle_non_tabu);
            } else if (non_tabu_change == -1) {
                this->delete_ab(selected_ab_cycle_non_tabu);
            }
            moved_ab_cycle_[selected_ab_cycle_non_tabu] = iteration + random_integer(1, t_max_);
        }
        if (improving_change == 0)
            ++stagnation;
        if (stagnation == max_stagnation_)
            break;
    }
}

template <typename Distances>
void Cross<Distances>::add_ab(int num)
{
    number_of_c_nodes_ += weight_c_[num] - 2 * weight_sr_[num];
    number_of_e_edges_ += ab_cycle_[num][0] / 2;

    used_ab_cycle_[num] = 1;
    ++number_of_used_ab_cycles_;
    for (int s1 = 0; s1 < number_of_ab_cycles_; ++s1)
        weight_sr_[s1] += weight_rr_[s1][num];
}

template <typename Distances>
void Cross<Distances>::delete_ab(int num)
{
    number_of_c_nodes_ -= weight_c_[num] - 2 * weight_sr_[num];
    number_of_e_edges_ -= ab_cycle_[num][0] / 2;

    used_ab_cycle_[num] = 0;
    --number_of_used_ab_cycles_;
    for (int s1 = 0; s1 < number_of_ab_cycles_; ++s1)
        weight_sr_[s1] -= weight_rr_[s1][num];
}

/**
 * The genetic algorithm's main loop: maintains a population of individuals,
 * combining pairs of them with 'Cross' and locally optimizing the result
 * with 'KOpt' each generation.
 */
template <typename Distances>
class Environment
{

public:

    /** Constructor. */
    Environment(
            const Distances& distances,
            VertexId number_of_vertices,
            int population_size,
            int number_of_children,
            double time_limit);

    /** Run the genetic algorithm until termination. */
    void run();

    /** Best tour found, as a 0-indexed list of vertices. */
    std::vector<VertexId> get_best_tour() const { return evaluator_.get_tour(best_individual_); }

private:

    /** Reset the generation counters and eset-selection strategy at the start of a run. */
    void reset_state();

    /** Whether the algorithm should stop (time limit, population convergence, or stagnation). */
    bool termination_condition();

    /** Update 'average_value_'/'best_value_'/'best_individual_' from the current population. */
    void set_average_best();

    /** Set the population to random tours, then locally optimize each of them. */
    void init_population();

    /** Draw a random pairing of the population for crossover. */
    void select_for_mating();

    /** Cross the pair of individuals at 'index_for_mating_[s]'/'[s + 1]'. */
    void generate_kids(
            int s);

    /** Recompute 'edge_frequency_' from the current population. */
    void compute_edge_frequencies();

    /** Evaluator, used to look up distances and extract the final tour. */
    Evaluator<Distances> evaluator_;

    /** Crossover operator. */
    Cross<Distances> cross_;

    /** Local search operator. */
    KOpt<Distances> kopt_;

    /** Population size. */
    int population_size_;

    /** Number of children generated per generation. */
    int number_of_children_;

    /** Current population. */
    std::vector<Individual> population_;

    /** Best individual found so far. */
    Individual best_individual_;

    /** Number of generations elapsed for the current population. */
    int current_number_of_generations_ = 0;

    /** Accumulated number of children generated so far. */
    long int accumulated_number_of_children_ = 0;

    /** Edge frequency across the population. */
    std::vector<std::vector<int>> edge_frequency_;

    /** Average tour length in the population. */
    double average_value_ = 0.0;

    /** Tour length of the best individual in the population. */
    Distance best_value_ = 0;

    /** Index of the best individual in the population. */
    int best_index_ = 0;

    /** Random pairing of the population for the current generation's crossovers. */
    std::vector<int> index_for_mating_;

    /** Number of consecutive generations without an improvement to 'best_individual_'. */
    int stagnation_count_ = 0;

    /** EAX method and eset-selection strategy (see 'Cross::run()'): '[0]' 1:Greedy 3:Distance 4:Entropy, '[1]' 1:Single-AB 2:Block2. */
    int flags_[10];

    /** Current stage of the genetic algorithm (1: Single-AB, 2: Block2). */
    int stage_ = 1;

    /** 'stagnation_count_' threshold at which the algorithm moves to the next stage. */
    int max_stagnation_ = 0;

    /** Number of generations elapsed when stage 1 completed. */
    int number_of_generations_stage_1_ = 0;

    /** Time at which 'run()' started. */
    clock_t time_start_ = 0;

    /** Wall-clock time limit in seconds. */
    double time_limit_;

};

template <typename Distances>
Environment<Distances>::Environment(
        const Distances& distances,
        VertexId number_of_vertices,
        int population_size,
        int number_of_children,
        double time_limit):
    evaluator_(distances, number_of_vertices),
    cross_(evaluator_, number_of_vertices, population_size),
    kopt_(evaluator_, number_of_vertices),
    population_size_(population_size),
    number_of_children_(number_of_children),
    population_(population_size, Individual(number_of_vertices)),
    best_individual_(number_of_vertices),
    edge_frequency_(number_of_vertices, std::vector<int>(number_of_vertices)),
    index_for_mating_(population_size + 1),
    time_limit_(time_limit)
{
}

template <typename Distances>
void Environment<Distances>::run()
{
    time_start_ = clock();
    init_population();
    reset_state();

    compute_edge_frequencies();
    while (true) {
        set_average_best();
        if (termination_condition())
            break;

        select_for_mating();
        for (int s = 0; s < population_size_; ++s)
            generate_kids(s);

        ++current_number_of_generations_;
    }
}

template <typename Distances>
void Environment<Distances>::reset_state()
{
    accumulated_number_of_children_ = 0;
    current_number_of_generations_ = 0;
    stagnation_count_ = 0;
    max_stagnation_ = 0;
    stage_ = 1; // sets stage to 1
    flags_[0] = 4; // maintains population diversity  1:Greedy, 2:---, 3:Distance, 4:Entropy
    flags_[1] = 1; // the type of Eset: 1:Single-AB, 2:Block2
}

template <typename Distances>
bool Environment<Distances>::termination_condition()
{
    if ((double)(clock() - time_start_) / CLOCKS_PER_SEC >= time_limit_)
        return true;
    if (average_value_ - best_value_ < 0.001)
        return true;
    if (stage_ == 1) {
        if (stagnation_count_ == int(1500 / number_of_children_) && max_stagnation_ == 0) { // 1500/Nch
            max_stagnation_ = int(current_number_of_generations_ / 10); // max_stagnation_ = G/10
        } else if (max_stagnation_ != 0 && max_stagnation_ <= stagnation_count_) {
            stagnation_count_ = 0;
            max_stagnation_ = 0;
            number_of_generations_stage_1_ = current_number_of_generations_;
            flags_[1] = 2;
            stage_ = 2;
        }
        return false;
    }
    if (stage_ == 2) {
        if (stagnation_count_ == int(1500 / number_of_children_) && max_stagnation_ == 0) { // 1500/Nch
            max_stagnation_ = int((current_number_of_generations_ - number_of_generations_stage_1_) / 10); // max_stagnation_ = G/10
        } else if (max_stagnation_ != 0 && max_stagnation_ <= stagnation_count_) {
            return true;
        }
        return false;
    }

    return true;
}

template <typename Distances>
void Environment<Distances>::set_average_best()
{
    Distance stock_best = best_individual_.length;
    average_value_ = 0.0;
    best_index_ = 0;
    best_value_ = population_[0].length;
    for (int i = 0; i < population_size_; ++i) {
        average_value_ += population_[i].length;
        if (population_[i].length < best_value_) {
            best_index_ = i;
            best_value_ = population_[i].length;
        }
    }
    best_individual_ = population_[best_index_];
    average_value_ /= (double)population_size_;
    if (best_individual_.length < stock_best) {
        stagnation_count_ = 0;
    } else {
        ++stagnation_count_;
    }
}

template <typename Distances>
void Environment<Distances>::init_population()
{
    for (int i = 0; i < population_size_; ++i) {
        kopt_.make_random_solution(population_[i]); // randomly sets a route
        kopt_.run(population_[i]); // local search (2-opt neighborhood)
    }
}

template <typename Distances>
void Environment<Distances>::select_for_mating()
{
    random_permutation(index_for_mating_, population_size_, population_size_);
    index_for_mating_[population_size_] = index_for_mating_[0];
}

template <typename Distances>
void Environment<Distances>::generate_kids(
        int s)
{
    // population_[index_for_mating_[s]] gets replaced by the best solution found by 'cross_.run()'
    // 'edge_frequency_' gets updated at the same time
    cross_.set_parents(population_[index_for_mating_[s]], population_[index_for_mating_[s + 1]], flags_, number_of_children_);
    cross_.run(population_[index_for_mating_[s]], population_[index_for_mating_[s + 1]], number_of_children_, 1, flags_, edge_frequency_);
    accumulated_number_of_children_ += cross_.number_of_generated_children;
}

template <typename Distances>
void Environment<Distances>::compute_edge_frequencies()
{
    VertexId number_of_vertices = evaluator_.number_of_vertices;
    for (VertexId vertex_id_1 = 0; vertex_id_1 < number_of_vertices; ++vertex_id_1)
        for (VertexId vertex_id_2 = 0; vertex_id_2 < number_of_vertices; ++vertex_id_2)
            edge_frequency_[vertex_id_1][vertex_id_2] = 0;

    for (int i = 0; i < population_size_; ++i)
        for (VertexId vertex_id = 0; vertex_id < number_of_vertices; ++vertex_id) {
            VertexId neighbor_vertex_id_1 = population_[i].neighbors[vertex_id][0];
            VertexId neighbor_vertex_id_2 = population_[i].neighbors[vertex_id][1];
            ++edge_frequency_[vertex_id][neighbor_vertex_id_1];
            ++edge_frequency_[vertex_id][neighbor_vertex_id_2];
        }
}

}

template <typename Distances>
const Output eax(
        const Distances& distances,
        const Instance& instance,
        const EaxParameters& parameters)
{
    Output output(instance);
    AlgorithmFormatter algorithm_formatter(parameters, output);
    algorithm_formatter.start("EAX");
    algorithm_formatter.print_header();

    VertexId number_of_vertices = instance.number_of_vertices();

    // The population/local-search RNG is seeded here (matching what
    // 'InitURandom(seed)' did in the original vendored code); the sort/RNG
    // helper objects the original code lazily allocated as process-wide
    // globals ('tRand'/'tSort') held no state of their own, so this
    // integration replaced them with plain functions ('random_*'/'sort_*'
    // above) instead of carrying that non-reentrancy hazard forward.
    seed_random(parameters.seed);

    Environment<Distances> environment(
            distances,
            number_of_vertices,
            parameters.population_size,
            parameters.number_of_children,
            parameters.timer.remaining_time());
    environment.run();

    std::vector<VertexId> tour = environment.get_best_tour();

    // 'tour[0]' is always city 0 (the traversal in 'get_tour' starts there),
    // and the 'Solution' constructor already starts with vertex 0, so it
    // must be skipped here to avoid visiting it twice.
    Solution solution(instance);
    for (std::size_t i = 1; i < tour.size(); ++i)
        solution.add_vertex(distances, tour[i]);

    algorithm_formatter.update_solution(solution, "final solution");

    algorithm_formatter.end();
    return output;
}

}
