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

template <typename Distances>
class Eax
{

public:

    Eax(
            const Distances& distances,
            const Instance& instance,
            std::mt19937_64& generator,
            const EaxParameters& parameters,
            AlgorithmFormatter& algorithm_formatter,
            const Output& output):
        distances_(distances),
        instance_(instance),
        generator_(generator),
        parameters_(parameters),
        algorithm_formatter_(algorithm_formatter),
        output_(output),
        number_of_vertices_(instance.number_of_vertices()),
        edge_frequencies_(number_of_vertices_),
        remaining_a_(number_of_vertices_),
        remaining_b_(number_of_vertices_),
        trace_stamps_(number_of_vertices_, 0),
        occurrences_(number_of_vertices_),
        is_cut_(number_of_vertices_, 0),
        links_(number_of_vertices_, {-1, -1}),
        in_subtour_(number_of_vertices_, 0),
        cached_neighbors_(number_of_vertices_),
        cached_neighbor_distances_(number_of_vertices_),
        neighbor_stamps_(number_of_vertices_, 0),
        ineffective_(number_of_vertices_, 0),
        visits_(number_of_vertices_),
        is_linked_(number_of_vertices_, 0)
    {
    }

    void run()
    {
        initialize_population();
        if (number_of_vertices_ < 8)
            return;
        run_stage(Strategy::single);
        if (parameters_.stage_2)
            run_stage(Strategy::block2);
    }

private:

    /*
     * Tours
     */

    inline Distance distance(VertexId vertex_id_1, VertexId vertex_id_2) const
    {
        return distances_.distance(vertex_id_1, vertex_id_2);
    }

    Tour tour_from_solution(const Solution& solution) const
    {
        Tour tour;
        tour.neighbors.resize(number_of_vertices_);
        VertexId n = number_of_vertices_;
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
     * Make 'order_' and 'positions_' the order of the vertices along a tour
     * (computed if the tour has changed since the last time).
     */
    void compute_order(Tour& tour)
    {
        if (!tour.order_valid) {
            tour.order.resize(number_of_vertices_);
            tour.positions.resize(number_of_vertices_);
            VertexId previous_vertex_id = tour.neighbors[0][0];
            VertexId vertex_id = 0;
            for (VertexId pos = 0; pos < number_of_vertices_; ++pos) {
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
        order_ = tour.order.data();
        positions_ = tour.positions.data();
    }

    Solution solution_from_tour(const Tour& tour)
    {
        SolutionBuilder solution_builder(instance_);
        VertexId previous_vertex_id = tour.neighbors[0][0];
        VertexId vertex_id = 0;
        for (VertexId pos = 0; pos < number_of_vertices_; ++pos) {
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

    void initialize_population()
    {
        // Nearest neighbors, for the 2-opt and the merging of the subtours.
        VertexId number_of_two_opt_candidates = 50;
        CandidateLists two_opt_candidates = nearest_neighbor_candidates(
                instance_.distances(),
                number_of_two_opt_candidates);
        near_vertices_.resize(number_of_vertices_);
        near_distances_.resize(number_of_vertices_);
        for (VertexId vertex_id = 0; vertex_id < number_of_vertices_; ++vertex_id) {
            const auto& list = two_opt_candidates[vertex_id];
            near_vertices_[vertex_id] = list;
            near_distances_[vertex_id].resize(list.size());
            for (size_t pos = 0; pos < list.size(); ++pos)
                near_distances_[vertex_id][pos] = distance(vertex_id, list[pos]);
        }

        const std::string& initial_tours = parameters_.initial_tours;
        const std::string& initial_local_search = parameters_.initial_local_search;
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
            alpha_nearness_parameters.number_of_candidates = parameters_.number_of_candidates;
            if (number_of_vertices_ >= 10000)
                alpha_nearness_parameters.initial_period = 100;
            alpha_nearness = alpha_nearness_candidates(
                    instance_.distances(),
                    alpha_nearness_parameters).candidates;
        }

        for (int i = 0; i < parameters_.population_size; ++i) {
            if (i > 0 && parameters_.timer.needs_to_end())
                break;
            Solution solution(instance_);
            if (initial_tours == "random-permutation") {
                RandomPermutationParameters random_permutation_parameters;
                random_permutation_parameters.verbosity_level = 0;
                solution = random_permutation(instance_, generator_, random_permutation_parameters).solution;
            } else {
                RandomWalkParameters random_walk_parameters;
                random_walk_parameters.verbosity_level = 0;
                solution = random_walk(instance_, generator_, random_walk_parameters, &alpha_nearness).solution;
            }
            if (initial_local_search == "2-opt") {
                TwoOptParameters two_opt_parameters;
                two_opt_parameters.verbosity_level = 0;
                solution = two_opt(distances_, instance_, generator_, two_opt_parameters, &solution, &two_opt_candidates).solution;
            } else if (initial_local_search == "lin-kernighan") {
                LinKernighanParameters lin_kernighan_parameters;
                lin_kernighan_parameters.verbosity_level = 0;
                lin_kernighan_parameters.maximum_number_of_trials = 0;
                solution = lin_kernighan(distances_, instance_, generator_, lin_kernighan_parameters, &solution, &alpha_nearness).solution;
            }
            population_.push_back(tour_from_solution(solution));
            if (!output_.solution.feasible()
                    || solution.distance() < output_.solution.distance()) {
                algorithm_formatter_.update_solution(solution, "initial solution " + std::to_string(i));
            }
        }
        // Each edge is counted from its smaller end.
        for (const Tour& tour: population_)
            for (VertexId vertex_id = 0; vertex_id < number_of_vertices_; ++vertex_id)
                for (VertexId neighbor_id: tour.neighbors[vertex_id])
                    if (vertex_id < neighbor_id)
                        edge_frequencies_.add(vertex_id, neighbor_id, 1);
    }

    /*
     * AB-cycles
     */

    /**
     * Compute the AB-cycles of parents A (whose order is in 'order_') and B:
     * the edges of the multigraph G_AB (the edges of A and the edges of B,
     * common edges included twice) are partitioned into cycles alternating
     * edges of A and B, by random tracing.
     *
     * Each effective AB-cycle (with at least one edge of A not in B) is
     * stored in 'ab_cycles_' as its vertices c[0], ..., c[2m - 1]: the edges
     * (c[2i], c[2i + 1]) are edges of A, the edges (c[2i + 1], c[2i + 2])
     * (indices modulo 2m) edges of B. The edges of A of the ineffective
     * AB-cycles (made of common edges) are marked in 'ineffective_'.
     *
     * The tracing stops once 'maximum_number_of_ab_cycles' effective
     * AB-cycles have been formed ('-1': no limit). Each new tracing starting
     * from a random vertex, the AB-cycles formed are then biased towards the
     * largest ones.
     */
    void compute_ab_cycles(
            const Tour& tour_a,
            const Tour& tour_b,
            VertexId maximum_number_of_ab_cycles)
    {
        number_of_ab_cycles_ = 0;
        ineffective_stamp_++;
        // The remaining edges of the vertices are initialized when they're
        // first reached (the tracing may stop early).
        tracing_tour_a_ = &tour_a;
        tracing_tour_b_ = &tour_b;
        trace_stamp_++;
        number_of_remaining_edges_ = 2 * number_of_vertices_;

        std::vector<VertexId>& path = path_;
        path.clear();
        while (number_of_remaining_edges_ > 0) {
            if (path.empty()) {
                // Random vertex with remaining edges.
                VertexId vertex_id = -1;
                std::uniform_int_distribution<VertexId> distribution(0, number_of_vertices_ - 1);
                for (int attempt = 0; attempt < 32 && vertex_id == -1; ++attempt) {
                    VertexId candidate_id = distribution(generator_);
                    reach(candidate_id);
                    if (remaining_a_[candidate_id].size > 0)
                        vertex_id = candidate_id;
                }
                if (vertex_id == -1) {
                    // Few vertices left: the first one from a random position.
                    VertexId start = distribution(generator_);
                    for (VertexId offset = 0; offset < number_of_vertices_; ++offset) {
                        VertexId candidate_id = (start + offset) % number_of_vertices_;
                        reach(candidate_id);
                        if (remaining_a_[candidate_id].size > 0) {
                            vertex_id = candidate_id;
                            break;
                        }
                    }
                }
                path.push_back(vertex_id);
                occurrences_[vertex_id].size = 0;
                occurrences_[vertex_id].vertices[occurrences_[vertex_id].size++] = 0;
            }
            VertexId vertex_id = path.back();
            // Edge index of the next edge: even for an edge of A.
            bool edge_a = ((path.size() - 1) % 2 == 0);
            NeighborList& remaining = (edge_a)? remaining_a_[vertex_id]: remaining_b_[vertex_id];
            if (remaining.size == 0) {
                // Only possible for a path reduced to its start vertex.
                clear_path();
                continue;
            }
            int pick = (remaining.size == 1)? 0: random_bit();
            VertexId next_vertex_id = remaining.vertices[pick];
            reach(next_vertex_id);
            remove_remaining_edge(vertex_id, next_vertex_id, edge_a);
            number_of_remaining_edges_--;
            path.push_back(next_vertex_id);
            VertexId k = path.size() - 1;
            // Does the path close an AB-cycle at 'next_vertex_id'?
            NeighborList& occurrences = occurrences_[next_vertex_id];
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
                    if (number_of_ab_cycles_ == (VertexId)ab_cycles_.size())
                        ab_cycles_.emplace_back();
                    std::vector<VertexId>& cycle = ab_cycles_[number_of_ab_cycles_++];
                    cycle.clear();
                    for (VertexId i = first_a; i < k; ++i)
                        cycle.push_back(path[i]);
                    if (first_a != j)
                        cycle.push_back(path[j]);
                } else {
                    for (VertexId i = first_a; i < k; i += 2)
                        ineffective_[cut_position(path[i], path[i + 1])] = ineffective_stamp_;
                }
                // Truncate the path to path[0], ..., path[j].
                for (VertexId i = j + 1; i < k; ++i)
                    remove_occurrence(path[i], i);
                path.resize(j + 1);
                closed = true;
                break;
            }
            if (number_of_ab_cycles_ == maximum_number_of_ab_cycles && maximum_number_of_ab_cycles != -1)
                break;
            if (!closed)
                occurrences.vertices[occurrences.size++] = k;
        }
        clear_path();
    }

    inline int random_bit()
    {
        if (number_of_random_bits_ == 0) {
            random_bits_ = generator_();
            number_of_random_bits_ = 64;
        }
        int bit = random_bits_ & 1;
        random_bits_ >>= 1;
        number_of_random_bits_--;
        return bit;
    }

    struct NeighborList
    {
        std::array<VertexId, 2> vertices;
        int size = 0;
    };

    void remove_occurrence(VertexId vertex_id, VertexId index)
    {
        NeighborList& occurrences = occurrences_[vertex_id];
        for (int o = 0; o < occurrences.size; ++o) {
            if (occurrences.vertices[o] == index) {
                occurrences.vertices[o] = occurrences.vertices[occurrences.size - 1];
                occurrences.size--;
                return;
            }
        }
    }

    void clear_path()
    {
        for (VertexId vertex_id: path_)
            occurrences_[vertex_id].size = 0;
        path_.clear();
    }

    void remove_from_list(NeighborList& list, VertexId vertex_id)
    {
        for (int i = 0; i < list.size; ++i) {
            if (list.vertices[i] == vertex_id) {
                list.vertices[i] = list.vertices[list.size - 1];
                list.size--;
                return;
            }
        }
    }

    void remove_remaining_edge(VertexId vertex_id_1, VertexId vertex_id_2, bool edge_a)
    {
        if (edge_a) {
            remove_from_list(remaining_a_[vertex_id_1], vertex_id_2);
            remove_from_list(remaining_a_[vertex_id_2], vertex_id_1);
        } else {
            remove_from_list(remaining_b_[vertex_id_1], vertex_id_2);
            remove_from_list(remaining_b_[vertex_id_2], vertex_id_1);
        }
    }

    /** Initialize the remaining edges of a vertex, the first time it's reached. */
    inline void reach(VertexId vertex_id)
    {
        if (trace_stamps_[vertex_id] == trace_stamp_)
            return;
        trace_stamps_[vertex_id] = trace_stamp_;
        remaining_a_[vertex_id].vertices = tracing_tour_a_->neighbors[vertex_id];
        remaining_a_[vertex_id].size = 2;
        remaining_b_[vertex_id].vertices = tracing_tour_b_->neighbors[vertex_id];
        remaining_b_[vertex_id].size = 2;
        occurrences_[vertex_id].size = 0;
    }

    /*
     * Intermediate solutions
     *
     * An intermediate solution is pA (whose order along the tour is
     * 'order_') whose edges at the cut positions are removed ('cuts_': the
     * edge between the positions c and c + 1), and with additional edges
     * ('links_') between the ends of the resulting paths (the segments).
     */

    void reset_intermediate()
    {
        for (VertexId cut: cuts_)
            is_cut_[cut] = 0;
        cuts_.clear();
        for (VertexId vertex_id: linked_vertices_) {
            links_[vertex_id] = {-1, -1};
            is_linked_[vertex_id] = 0;
        }
        linked_vertices_.clear();
    }

    /** Position of the edge between consecutive vertices of pA (-1 if not consecutive). */
    inline VertexId cut_position(VertexId vertex_id_1, VertexId vertex_id_2) const
    {
        VertexId pos_1 = positions_[vertex_id_1];
        VertexId pos_2 = positions_[vertex_id_2];
        if (pos_1 + 1 == pos_2 || (pos_1 == number_of_vertices_ - 1 && pos_2 == 0))
            return pos_1;
        if (pos_2 + 1 == pos_1 || (pos_2 == number_of_vertices_ - 1 && pos_1 == 0))
            return pos_2;
        return -1;
    }

    void add_link(VertexId vertex_id_1, VertexId vertex_id_2)
    {
        for (auto [a, b]: {std::pair{vertex_id_1, vertex_id_2}, std::pair{vertex_id_2, vertex_id_1}}) {
            auto& links = links_[a];
            if (!is_linked_[a]) {
                is_linked_[a] = 1;
                linked_vertices_.push_back(a);
            }
            if (links[0] == -1) {
                links[0] = b;
            } else {
                links[1] = b;
            }
        }
    }

    void remove_link(VertexId vertex_id_1, VertexId vertex_id_2)
    {
        for (auto [a, b]: {std::pair{vertex_id_1, vertex_id_2}, std::pair{vertex_id_2, vertex_id_1}}) {
            auto& links = links_[a];
            if (links[0] == b) {
                links[0] = -1;
            } else {
                links[1] = -1;
            }
        }
    }

    /** Remove an edge of the intermediate solution. */
    void remove_edge(VertexId vertex_id_1, VertexId vertex_id_2)
    {
        VertexId cut = cut_position(vertex_id_1, vertex_id_2);
        if (cut != -1 && !is_cut_[cut]) {
            is_cut_[cut] = 1;
            cuts_.insert(std::upper_bound(cuts_.begin(), cuts_.end(), cut), cut);
        } else {
            remove_link(vertex_id_1, vertex_id_2);
        }
    }

    /** Get the two neighbors of a vertex in the intermediate solution. */
    inline std::array<VertexId, 2> neighbors(VertexId vertex_id) const
    {
        std::array<VertexId, 2> result;
        int size = 0;
        VertexId pos = positions_[vertex_id];
        VertexId previous_pos = (pos == 0)? number_of_vertices_ - 1: pos - 1;
        if (!is_cut_[previous_pos])
            result[size++] = order_[previous_pos];
        if (!is_cut_[pos])
            result[size++] = order_[(pos + 1 == number_of_vertices_)? 0: pos + 1];
        for (VertexId link: links_[vertex_id])
            if (link != -1 && size < 2)
                result[size++] = link;
        return result;
    }

    /** Segment containing a position: the index of the cut ending it. */
    inline VertexId segment(VertexId pos) const
    {
        VertexId index = std::lower_bound(cuts_.begin(), cuts_.end(), pos) - cuts_.begin();
        return (index == (VertexId)cuts_.size())? 0: index;
    }

    inline VertexId segment_start(VertexId segment_id) const
    {
        VertexId pos = (segment_id == 0)? cuts_.back() + 1: cuts_[segment_id - 1] + 1;
        return (pos == number_of_vertices_)? 0: pos;
    }

    inline VertexId segment_end(VertexId segment_id) const
    {
        return cuts_[segment_id];
    }

    inline VertexId segment_size(VertexId segment_id) const
    {
        VertexId start = segment_start(segment_id);
        VertexId end = segment_end(segment_id);
        return (end >= start)? end - start + 1: end - start + 1 + number_of_vertices_;
    }

    /** Compute the subtours of the intermediate solution. */
    void compute_subtours()
    {
        VertexId number_of_segments = cuts_.size();
        segment_subtours_.assign(number_of_segments, -1);
        subtour_sizes_.clear();
        for (VertexId segment_id_0 = 0; segment_id_0 < number_of_segments; ++segment_id_0) {
            if (segment_subtours_[segment_id_0] != -1)
                continue;
            VertexId subtour_id = subtour_sizes_.size();
            VertexId size = 0;
            VertexId segment_id = segment_id_0;
            VertexId entry_vertex_id = order_[segment_start(segment_id)];
            VertexId previous_vertex_id = -1;
            for (;;) {
                segment_subtours_[segment_id] = subtour_id;
                size += segment_size(segment_id);
                VertexId start_vertex_id = order_[segment_start(segment_id)];
                VertexId end_vertex_id = order_[segment_end(segment_id)];
                VertexId exit_vertex_id = (entry_vertex_id == start_vertex_id)?
                    end_vertex_id:
                    start_vertex_id;
                // Link leaving the segment.
                const auto& links = links_[exit_vertex_id];
                VertexId next_vertex_id;
                if (start_vertex_id == end_vertex_id) {
                    // Single vertex: leave by the link it wasn't entered by.
                    next_vertex_id = (links[0] != previous_vertex_id)? links[0]: links[1];
                } else {
                    next_vertex_id = (links[0] != -1)? links[0]: links[1];
                }
                previous_vertex_id = exit_vertex_id;
                entry_vertex_id = next_vertex_id;
                segment_id = segment(positions_[next_vertex_id]);
                if (segment_id == segment_id_0)
                    break;
            }
            subtour_sizes_.push_back(size);
        }
    }

    /**
     * Connect the subtours of the intermediate solution into a tour (step 5
     * of EAX): repeatedly, the smallest subtour is connected to another one
     * by the best exchange of an edge of it and an edge of the other
     * subtour, at least one end of which is among the nearest neighbors of
     * an end of the first edge.
     */
    void merge_subtours()
    {
        compute_subtours();
        while (subtour_sizes_.size() > 1) {
            VertexId subtour_id = std::min_element(subtour_sizes_.begin(), subtour_sizes_.end()) - subtour_sizes_.begin();
            Distance best_delta = std::numeric_limits<Distance>::max();
            std::array<VertexId, 4> best = {-1, -1, -1, -1};
            // Mark the vertices of the subtour.
            subtour_stamp_++;
            for (VertexId segment_id = 0; segment_id < (VertexId)cuts_.size(); ++segment_id) {
                if (segment_subtours_[segment_id] != subtour_id)
                    continue;
                VertexId pos = segment_start(segment_id);
                VertexId end = segment_end(segment_id);
                for (;;) {
                    in_subtour_[order_[pos]] = subtour_stamp_;
                    if (pos == end)
                        break;
                    pos = (pos + 1 == number_of_vertices_)? 0: pos + 1;
                }
            }
            // The nearest neighbors first, then more of them if none is
            // outside of the subtour.
            for (VertexId near_limit: {parameters_.number_of_near_vertices, (VertexId)50}) {
                for (VertexId segment_id = 0; segment_id < (VertexId)cuts_.size(); ++segment_id) {
                    if (segment_subtours_[segment_id] != subtour_id)
                        continue;
                    VertexId pos = segment_start(segment_id);
                    VertexId end = segment_end(segment_id);
                    for (;;) {
                        // Edges (u, u_1) and (u, u_2) of the subtour; edges (w,
                        // w_1) and (w, w_2) of another subtour, w being near u.
                        VertexId u = order_[pos];
                        std::array<VertexId, 2> u_neighbors = neighbors(u);
                        Distance d_u_0 = distance(u, u_neighbors[0]);
                        Distance d_u_1 = distance(u, u_neighbors[1]);
                        const std::vector<VertexId>& near_vertices = near_vertices_[u];
                        const std::vector<Distance>& near_distances = near_distances_[u];
                        VertexId number_of_near_vertices = std::min(near_limit, (VertexId)near_vertices.size());
                        for (VertexId near_pos = 0; near_pos < number_of_near_vertices; ++near_pos) {
                            VertexId w = near_vertices[near_pos];
                            if (in_subtour_[w] == subtour_stamp_)
                                continue;
                            if (neighbor_stamps_[w] != subtour_stamp_) {
                                neighbor_stamps_[w] = subtour_stamp_;
                                cached_neighbors_[w] = neighbors(w);
                                cached_neighbor_distances_[w] = {
                                    distance(w, cached_neighbors_[w][0]),
                                    distance(w, cached_neighbors_[w][1])};
                            }
                            const std::array<VertexId, 2>& w_neighbors = cached_neighbors_[w];
                            const std::array<Distance, 2>& w_distances = cached_neighbor_distances_[w];
                            Distance d_uw = near_distances[near_pos];
                            for (int i = 0; i < 2; ++i) {
                                VertexId u_2 = u_neighbors[i];
                                Distance d_u = (i == 0)? d_u_0: d_u_1;
                                Distance d_u2_w = distance(u_2, w);
                                for (int j = 0; j < 2; ++j) {
                                    VertexId w_2 = w_neighbors[j];
                                    Distance d_removed = d_u + w_distances[j];
                                    // Add (u, w) and (u_2, w_2). (The length of
                                    // the second added edge is only computed
                                    // if the move can be better than the best
                                    // one, and the same below.)
                                    if (best_delta > d_uw - d_removed) {
                                        Distance delta_1 = d_uw + distance(u_2, w_2) - d_removed;
                                        if (best_delta > delta_1) {
                                            best_delta = delta_1;
                                            best = {u, u_2, w, w_2};
                                        }
                                    }
                                    // Add (u, w_2) and (u_2, w).
                                    if (best_delta > d_u2_w - d_removed) {
                                        Distance delta_2 = distance(u, w_2) + d_u2_w - d_removed;
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
                        pos = (pos + 1 == number_of_vertices_)? 0: pos + 1;
                    }
                }
                if (best[0] != -1)
                    break;
            }
            if (best[0] == -1) {
                // No near vertex in another subtour: connect a vertex of the
                // subtour to the nearest vertex outside of it.
                VertexId vertex_id = -1;
                for (VertexId segment_id = 0; segment_id < (VertexId)cuts_.size(); ++segment_id)
                    if (segment_subtours_[segment_id] == subtour_id)
                        vertex_id = order_[segment_start(segment_id)];
                VertexId vertex_id_2 = neighbors(vertex_id)[0];
                for (VertexId vertex_id_3 = 0; vertex_id_3 < number_of_vertices_; ++vertex_id_3) {
                    if (in_subtour_[vertex_id_3] == subtour_stamp_)
                        continue;
                    for (VertexId vertex_id_4: neighbors(vertex_id_3)) {
                        Distance delta = distance(vertex_id, vertex_id_3)
                            + distance(vertex_id_2, vertex_id_4)
                            - distance(vertex_id, vertex_id_2)
                            - distance(vertex_id_3, vertex_id_4);
                        if (best_delta > delta) {
                            best_delta = delta;
                            best = {vertex_id, vertex_id_2, vertex_id_3, vertex_id_4};
                        }
                    }
                }
            }
            // Remove (0, 1) and (2, 3), add (0, 2) and (1, 3).
            remove_edge(best[0], best[1]);
            remove_edge(best[2], best[3]);
            add_link(best[0], best[2]);
            add_link(best[1], best[3]);
            compute_subtours();
        }
    }

    /*
     * Offspring
     */

    inline double entropy_term(int64_t frequency) const
    {
        if (frequency <= 0)
            return 0;
        double x = (double)frequency / parameters_.population_size;
        return -x * std::log(x);
    }

    /**
     * Build the offspring solution of the current intermediate solution
     * (after 'merge_subtours'): its edges removed from and added to pA, and
     * its evaluation.
     */
    void build_offspring(Offspring& offspring)
    {
        offspring.removed_edges.clear();
        offspring.added_edges.clear();
        for (VertexId cut: cuts_) {
            offspring.removed_edges.push_back(Edge(
                        order_[cut],
                        order_[(cut + 1 == number_of_vertices_)? 0: cut + 1]));
        }
        for (VertexId vertex_id: linked_vertices_)
            for (VertexId link: links_[vertex_id])
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
            length_difference -= distance(edge.vertex_id_1, edge.vertex_id_2);
            int64_t frequency = edge_frequencies_.get(edge.vertex_id_1, edge.vertex_id_2);
            entropy_difference += entropy_term(frequency - 1) - entropy_term(frequency);
        }
        for (const Edge& edge: added) {
            length_difference += distance(edge.vertex_id_1, edge.vertex_id_2);
            int64_t frequency = edge_frequencies_.get(edge.vertex_id_1, edge.vertex_id_2);
            entropy_difference += entropy_term(frequency + 1) - entropy_term(frequency);
        }
        offspring.length_difference = length_difference;
        // Entropy-preserving selection (the average length difference is
        // proportional to the length difference).
        double average_length_difference = (double)length_difference / parameters_.population_size;
        if (length_difference >= 0) {
            offspring.evaluation = -average_length_difference;
        } else if (entropy_difference < 0) {
            offspring.evaluation = average_length_difference / entropy_difference;
        } else {
            offspring.evaluation = -average_length_difference / 1e-9;
        }
    }

    /** Replace a tour of the population by an offspring solution of it. */
    void apply_offspring(Tour& tour, const Offspring& offspring)
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
            edge_frequencies_.add(edge.vertex_id_1, edge.vertex_id_2, -1);
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
            edge_frequencies_.add(edge.vertex_id_1, edge.vertex_id_2, 1);
        }
        tour.length += offspring.length_difference;
        tour.order_valid = false;
#ifndef NDEBUG
        // The tour must be a Hamiltonian cycle of the recorded length.
        Distance length = 0;
        VertexId previous_vertex_id = tour.neighbors[0][0];
        VertexId vertex_id = 0;
        for (VertexId pos = 0; pos < number_of_vertices_; ++pos) {
            VertexId next_vertex_id = (tour.neighbors[vertex_id][0] != previous_vertex_id)?
                tour.neighbors[vertex_id][0]:
                tour.neighbors[vertex_id][1];
            length += distance(vertex_id, next_vertex_id);
            previous_vertex_id = vertex_id;
            vertex_id = next_vertex_id;
            assert(pos == number_of_vertices_ - 1 || vertex_id != 0);
        }
        assert(vertex_id == 0);
        assert(length == tour.length);
#endif
    }

    /*
     * Crossover
     */

    /** Apply the AB-cycles of an E-set to pA, connect the subtours, and build the offspring. */
    void generate_offspring(
            const std::vector<VertexId>& e_set,
            Offspring& offspring)
    {
        reset_intermediate();
        for (VertexId ab_cycle_id: e_set) {
            const std::vector<VertexId>& cycle = ab_cycles_[ab_cycle_id];
            VertexId size = cycle.size();
            for (VertexId i = 0; i < size; i += 2) {
                VertexId cut = cut_position(cycle[i], cycle[i + 1]);
                is_cut_[cut] = 1;
                cuts_.push_back(cut);
                add_link(cycle[i + 1], cycle[(i + 2 == size)? 0: i + 2]);
            }
        }
        std::sort(cuts_.begin(), cuts_.end());
        merge_subtours();
        build_offspring(offspring);
    }

    /**
     * Generate offspring solutions of pA and pB, and replace pA by the best
     * one if it improves it (entropy-preserving selection), unless its
     * length is the one of pB (it's then very likely pB itself, and
     * replacing pA by it would decrease the diversity). Return 'true' iff pA
     * has been replaced.
     */
    bool crossover(
            Tour& tour_a,
            const Tour& tour_b,
            Strategy strategy)
    {
        compute_order(tour_a);
        compute_ab_cycles(
                tour_a,
                tour_b,
                (strategy == Strategy::single)? parameters_.number_of_children: -1);
        if (number_of_ab_cycles_ == 0)
            return false;

        best_offspring_.evaluation = 0;
        bool found = false;
        std::vector<VertexId> e_set;
        if (strategy == Strategy::single) {
            // Each offspring from a single AB-cycle, without repetition.
            std::vector<VertexId> ab_cycle_ids(number_of_ab_cycles_);
            std::iota(ab_cycle_ids.begin(), ab_cycle_ids.end(), 0);
            std::shuffle(ab_cycle_ids.begin(), ab_cycle_ids.end(), generator_);
            VertexId number_of_children = std::min(
                    (VertexId)parameters_.number_of_children,
                    number_of_ab_cycles_);
            for (VertexId child_id = 0; child_id < number_of_children; ++child_id) {
                e_set = {ab_cycle_ids[child_id]};
                generate_offspring(e_set, offspring_);
                if (offspring_.evaluation > best_offspring_.evaluation
                        && tour_a.length + offspring_.length_difference != tour_b.length) {
                    std::swap(best_offspring_, offspring_);
                    found = true;
                }
            }
        } else {
            for (VertexId child_id = 0; child_id < parameters_.number_of_children; ++child_id) {
                if (!block2_e_set(child_id, e_set))
                    break;
                generate_offspring(e_set, offspring_);
                if (offspring_.evaluation > best_offspring_.evaluation
                        && tour_a.length + offspring_.length_difference != tour_b.length) {
                    std::swap(best_offspring_, offspring_);
                    found = true;
                }
            }
        }
        reset_intermediate();
        if (found)
            apply_offspring(tour_a, best_offspring_);
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
    void prepare_block2()
    {
        VertexId number_of_ab_cycles = number_of_ab_cycles_;
        // Visits of the vertices by the AB-cycles (at most two).
        for (VertexId vertex_id: visited_vertices_)
            visits_[vertex_id].size = 0;
        visited_vertices_.clear();
        auto add_visit = [this](VertexId vertex_id, VertexId ab_cycle_id)
        {
            NeighborList& visits = visits_[vertex_id];
            if (visits.size == 0)
                visited_vertices_.push_back(vertex_id);
            visits.vertices[visits.size++] = ab_cycle_id;
        };
        for (VertexId ab_cycle_id = 0; ab_cycle_id < number_of_ab_cycles; ++ab_cycle_id)
            for (VertexId vertex_id: ab_cycles_[ab_cycle_id])
                add_visit(vertex_id, ab_cycle_id);

        // Sequences of successive ineffective AB-cycles (of common edges),
        // along pA: each one is integrated into an effective AB-cycle
        // visiting one of its ends (each end is visited once by an effective
        // AB-cycle).
        auto is_common = [this](VertexId pos)
        {
            return ineffective_[pos] == ineffective_stamp_;
        };
        // Start after a non-common edge.
        VertexId start_pos = 0;
        while (is_common(start_pos))
            start_pos++;
        start_pos = (start_pos + 1 == number_of_vertices_)? 0: start_pos + 1;
        VertexId pos = start_pos;
        for (VertexId count = 0; count < number_of_vertices_;) {
            if (!is_common(pos)) {
                pos = (pos + 1 == number_of_vertices_)? 0: pos + 1;
                count++;
                continue;
            }
            // Sequence of common edges from the vertex at 'pos'.
            VertexId first_vertex_id = order_[pos];
            while (count < number_of_vertices_ && is_common(pos)) {
                pos = (pos + 1 == number_of_vertices_)? 0: pos + 1;
                count++;
            }
            VertexId last_vertex_id = order_[pos];
            if (visits_[first_vertex_id].size != 1 || visits_[last_vertex_id].size != 1)
                continue;
            VertexId first_ab_cycle_id = visits_[first_vertex_id].vertices[0];
            VertexId last_ab_cycle_id = visits_[last_vertex_id].vertices[0];
            VertexId ab_cycle_id = (std::uniform_int_distribution<int>(0, 1)(generator_) == 0)?
                first_ab_cycle_id:
                last_ab_cycle_id;
            add_visit(first_vertex_id, ab_cycle_id);
            add_visit(last_vertex_id, ab_cycle_id);
        }

        // n_i and w_ij.
        c_vertices_.assign(number_of_ab_cycles, 0);
        contacts_.resize(number_of_ab_cycles);
        for (VertexId ab_cycle_id = 0; ab_cycle_id < number_of_ab_cycles; ++ab_cycle_id)
            contacts_[ab_cycle_id].clear();
        for (VertexId vertex_id: visited_vertices_) {
            const NeighborList& visits = visits_[vertex_id];
            if (visits.size == 1) {
                c_vertices_[visits.vertices[0]]++;
            } else if (visits.vertices[0] != visits.vertices[1]) {
                VertexId ab_cycle_id_1 = visits.vertices[0];
                VertexId ab_cycle_id_2 = visits.vertices[1];
                c_vertices_[ab_cycle_id_1]++;
                c_vertices_[ab_cycle_id_2]++;
                add_contact(ab_cycle_id_1, ab_cycle_id_2);
                add_contact(ab_cycle_id_2, ab_cycle_id_1);
            }
        }

        // AB-cycles by decreasing size (the central AB-cycles).
        ab_cycles_by_size_.resize(number_of_ab_cycles);
        std::iota(ab_cycles_by_size_.begin(), ab_cycles_by_size_.end(), 0);
        std::sort(
                ab_cycles_by_size_.begin(),
                ab_cycles_by_size_.end(),
                [this](VertexId ab_cycle_id_1, VertexId ab_cycle_id_2)
                {
                    return ab_cycles_[ab_cycle_id_1].size() > ab_cycles_[ab_cycle_id_2].size();
                });
        in_e_set_.assign(number_of_ab_cycles, 0);
        contact_sums_.assign(number_of_ab_cycles, 0);
        tabu_.assign(number_of_ab_cycles, 0);
        candidate_positions_.assign(number_of_ab_cycles, -1);
        candidates_.clear();
        e_set_members_.clear();
    }

    void add_contact(VertexId ab_cycle_id_1, VertexId ab_cycle_id_2)
    {
        for (auto& contact: contacts_[ab_cycle_id_1]) {
            if (contact.first == ab_cycle_id_2) {
                contact.second++;
                return;
            }
        }
        contacts_[ab_cycle_id_1].push_back({ab_cycle_id_2, 1});
    }

    /** Add an AB-cycle to the E-set of the tabu search (or remove it). */
    void block2_move(VertexId ab_cycle_id, bool add)
    {
        int64_t sign = (add)? 1: -1;
        number_of_c_vertices_ += sign * (c_vertices_[ab_cycle_id] - 2 * contact_sums_[ab_cycle_id]);
        in_e_set_[ab_cycle_id] = add;
        if (add) {
            e_set_members_.push_back(ab_cycle_id);
            remove_candidate(ab_cycle_id);
        } else {
            e_set_members_.erase(std::find(e_set_members_.begin(), e_set_members_.end(), ab_cycle_id));
            if (contact_sums_[ab_cycle_id] > 0)
                add_candidate(ab_cycle_id);
        }
        for (const auto& contact: contacts_[ab_cycle_id]) {
            VertexId other_ab_cycle_id = contact.first;
            contact_sums_[other_ab_cycle_id] += sign * contact.second;
            if (!in_e_set_[other_ab_cycle_id]) {
                if (contact_sums_[other_ab_cycle_id] > 0) {
                    add_candidate(other_ab_cycle_id);
                } else {
                    remove_candidate(other_ab_cycle_id);
                }
            }
        }
    }

    void add_candidate(VertexId ab_cycle_id)
    {
        if (candidate_positions_[ab_cycle_id] != -1)
            return;
        candidate_positions_[ab_cycle_id] = candidates_.size();
        candidates_.push_back(ab_cycle_id);
    }

    void remove_candidate(VertexId ab_cycle_id)
    {
        VertexId pos = candidate_positions_[ab_cycle_id];
        if (pos == -1)
            return;
        VertexId last_ab_cycle_id = candidates_.back();
        candidates_[pos] = last_ab_cycle_id;
        candidate_positions_[last_ab_cycle_id] = pos;
        candidates_.pop_back();
        candidate_positions_[ab_cycle_id] = -1;
    }

    /**
     * E-set of the child 'child_id' with the block2 strategy: a tabu search
     * minimizing the number of C-vertices of the intermediate solution (an
     * estimate of its number of subtours), from an E-set made of the central
     * AB-cycle (the 'child_id'-th largest one) and of some of the smaller
     * AB-cycles in contact with it. Return 'false' if there's no AB-cycle
     * left to be central.
     */
    bool block2_e_set(VertexId child_id, std::vector<VertexId>& e_set)
    {
        if (child_id == 0)
            prepare_block2();
        if (child_id >= (VertexId)ab_cycles_by_size_.size())
            return false;
        // Reset the E-set.
        for (VertexId ab_cycle_id: std::vector<VertexId>(e_set_members_))
            block2_move(ab_cycle_id, false);
        for (VertexId ab_cycle_id: std::vector<VertexId>(candidates_))
            remove_candidate(ab_cycle_id);
        for (VertexId ab_cycle_id = 0; ab_cycle_id < number_of_ab_cycles_; ++ab_cycle_id)
            tabu_[ab_cycle_id] = 0;
        number_of_c_vertices_ = 0;

        VertexId central_ab_cycle_id = ab_cycles_by_size_[child_id];
        VertexId central_size = ab_cycles_[central_ab_cycle_id].size();
        block2_move(central_ab_cycle_id, true);
        for (const auto& contact: std::vector<std::pair<VertexId, int64_t>>(contacts_[central_ab_cycle_id])) {
            VertexId ab_cycle_id = contact.first;
            if ((VertexId)ab_cycles_[ab_cycle_id].size() < central_size
                    && std::uniform_int_distribution<int>(0, 1)(generator_) == 0) {
                block2_move(ab_cycle_id, true);
            }
        }

        int64_t best_number_of_c_vertices = number_of_c_vertices_;
        std::vector<VertexId> best_e_set = e_set_members_;
        int64_t iteration = 0;
        int64_t iterations_without_improvement = 0;
        while (iterations_without_improvement < parameters_.tabu_maximum_number_of_iterations_without_improvement) {
            // Best move (random tie-breaking).
            VertexId best_ab_cycle_id = -1;
            bool best_add = false;
            int64_t best_delta = std::numeric_limits<int64_t>::max();
            int64_t number_of_ties = 0;
            auto consider = [&](VertexId ab_cycle_id, bool add)
            {
                int64_t delta = (add)?
                    c_vertices_[ab_cycle_id] - 2 * contact_sums_[ab_cycle_id]:
                    -c_vertices_[ab_cycle_id] + 2 * contact_sums_[ab_cycle_id];
                bool tabu = (iteration <= tabu_[ab_cycle_id]);
                if (tabu && number_of_c_vertices_ + delta >= best_number_of_c_vertices)
                    return;
                if (delta < best_delta) {
                    best_delta = delta;
                    best_ab_cycle_id = ab_cycle_id;
                    best_add = add;
                    number_of_ties = 1;
                } else if (delta == best_delta) {
                    number_of_ties++;
                    if (std::uniform_int_distribution<int64_t>(0, number_of_ties - 1)(generator_) == 0) {
                        best_ab_cycle_id = ab_cycle_id;
                        best_add = add;
                    }
                }
            };
            for (VertexId ab_cycle_id: candidates_)
                consider(ab_cycle_id, true);
            for (VertexId ab_cycle_id: e_set_members_)
                if (ab_cycle_id != central_ab_cycle_id)
                    consider(ab_cycle_id, false);
            if (best_ab_cycle_id == -1)
                break;
            block2_move(best_ab_cycle_id, best_add);
            tabu_[best_ab_cycle_id] = iteration
                + std::uniform_int_distribution<int64_t>(0, parameters_.tabu_maximum_tenure)(generator_);
            iteration++;
            if (number_of_c_vertices_ < best_number_of_c_vertices) {
                best_number_of_c_vertices = number_of_c_vertices_;
                best_e_set = e_set_members_;
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

    Distance best_length() const
    {
        Distance best = population_[0].length;
        for (const Tour& tour: population_)
            best = std::min(best, tour.length);
        return best;
    }

    void report_best(const std::string& comment)
    {
        VertexId best_id = 0;
        for (VertexId i = 0; i < (VertexId)population_.size(); ++i)
            if (population_[i].length < population_[best_id].length)
                best_id = i;
        if (!output_.solution.feasible()
                || population_[best_id].length < output_.solution.distance()) {
            algorithm_formatter_.update_solution(solution_from_tour(population_[best_id]), comment);
        }
    }

    /**
     * Run a stage: generations until the best solution doesn't improve over
     * G / 10 generations, G being the number of generations at which it
     * first didn't improve over 1500 / N_ch generations.
     */
    void run_stage(Strategy strategy)
    {
        VertexId population_size = population_.size();
        std::vector<VertexId> permutation(population_size);
        std::iota(permutation.begin(), permutation.end(), 0);
        Distance best = best_length();
        int64_t generation = 0;
        int64_t stagnation = 0;
        int64_t g = -1;
        int64_t stagnation_limit = 1500 / parameters_.number_of_children;
        std::string stage_name = (strategy == Strategy::single)? "stage 1": "stage 2";
        while (!parameters_.timer.needs_to_end()) {
            std::shuffle(permutation.begin(), permutation.end(), generator_);
            for (VertexId i = 0; i < population_size; ++i) {
                if (parameters_.timer.needs_to_end())
                    break;
                Tour& tour_a = population_[permutation[i]];
                const Tour& tour_b = population_[permutation[(i + 1) % population_size]];
                crossover(tour_a, tour_b, strategy);
            }
            generation++;
            Distance new_best = best_length();
            if (new_best < best) {
                best = new_best;
                stagnation = 0;
                report_best(stage_name + ", generation " + std::to_string(generation));
            } else {
                stagnation++;
            }
            if (g == -1 && stagnation >= stagnation_limit)
                g = generation;
            // Converged population.
            double average_length = 0;
            for (const Tour& tour: population_)
                average_length += tour.length;
            average_length /= population_size;
            if (average_length - best < 0.001)
                break;
            if (g != -1 && stagnation >= g / 10)
                break;
        }
        report_best(stage_name + ", end");
    }

    /*
     * Attributes
     */

    const Distances& distances_;

    const Instance& instance_;

    std::mt19937_64& generator_;

    const EaxParameters& parameters_;

    AlgorithmFormatter& algorithm_formatter_;

    const Output& output_;

    VertexId number_of_vertices_;

    /** Nearest neighbors, for the connection of the subtours, and their distances. */
    std::vector<std::vector<VertexId>> near_vertices_;
    std::vector<std::vector<Distance>> near_distances_;

    /**
     * Neighbors in the intermediate solution of the vertices near the
     * subtour being connected, and the lengths of their edges (valid if
     * 'neighbor_stamps_' is 'subtour_stamp_').
     */
    std::vector<std::array<VertexId, 2>> cached_neighbors_;
    std::vector<std::array<Distance, 2>> cached_neighbor_distances_;
    std::vector<int64_t> neighbor_stamps_;

    std::vector<Tour> population_;

    EdgeFrequencies edge_frequencies_;

    /** Order of the vertices along pA ('Tour::order' of pA). */
    const int32_t* order_ = nullptr;

    /** Position of each vertex along pA ('Tour::positions' of pA). */
    const int32_t* positions_ = nullptr;

    /** Edges of pA (by position) in ineffective AB-cycles. */
    std::vector<int64_t> ineffective_;
    int64_t ineffective_stamp_ = 0;

    /** Edges not yet in an AB-cycle. */
    std::vector<NeighborList> remaining_a_;
    std::vector<NeighborList> remaining_b_;

    /** Tracing: parents, and vertices reached (stamps). */
    const Tour* tracing_tour_a_ = nullptr;
    const Tour* tracing_tour_b_ = nullptr;
    std::vector<int64_t> trace_stamps_;
    int64_t trace_stamp_ = 0;
    VertexId number_of_remaining_edges_ = 0;

    /** Traced path. */
    std::vector<VertexId> path_;

    /** Indices of each vertex in the traced path. */
    std::vector<NeighborList> occurrences_;

    /** Effective AB-cycles ('number_of_ab_cycles_' first ones). */
    std::vector<std::vector<VertexId>> ab_cycles_;
    VertexId number_of_ab_cycles_ = 0;

    /** Random bits, drawn 64 at a time. */
    uint64_t random_bits_ = 0;
    int number_of_random_bits_ = 0;

    /** Intermediate solution. */
    std::vector<char> is_cut_;
    std::vector<VertexId> cuts_;
    std::vector<std::array<VertexId, 2>> links_;
    /** Vertices which have had links (each once). */
    std::vector<VertexId> linked_vertices_;
    std::vector<char> is_linked_;
    std::vector<VertexId> segment_subtours_;
    std::vector<VertexId> subtour_sizes_;

    /** Marks of the vertices of the subtour being connected. */
    std::vector<int64_t> in_subtour_;
    int64_t subtour_stamp_ = 0;

    Offspring offspring_;
    Offspring best_offspring_;

    /** Block2 strategy. */
    std::vector<NeighborList> visits_;
    std::vector<VertexId> visited_vertices_;
    std::vector<int64_t> c_vertices_;
    std::vector<std::vector<std::pair<VertexId, int64_t>>> contacts_;
    std::vector<VertexId> ab_cycles_by_size_;
    std::vector<char> in_e_set_;
    std::vector<int64_t> contact_sums_;
    std::vector<int64_t> tabu_;
    std::vector<VertexId> candidates_;
    std::vector<VertexId> candidate_positions_;
    std::vector<VertexId> e_set_members_;
    int64_t number_of_c_vertices_ = 0;
};

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

    Eax<Distances> eax(distances, instance, generator, parameters, algorithm_formatter, output);
    eax.run();

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
