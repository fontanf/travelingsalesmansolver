#pragma once

#include "travelingsalesmansolver/candidates/candidate_lists.hpp"
#include "travelingsalesmansolver/distances/distances.hpp"

namespace travelingsalesmansolver
{

struct AlphaNearnessParameters
{
    /** Number of candidates per vertex. */
    VertexId number_of_candidates = 5;

    /**
     * Precision of the transformed costs: an edge (i, j) costs
     * 'precision * d(i, j) + pi[i] + pi[j]', the penalties 'pi' being
     * integers.
     */
    Distance precision = 100;

    /**
     * Graph on which the 1-trees and the alpha values are computed: each
     * vertex is linked to its 'number_of_nearest_neighbors' nearest neighbors
     * and, for coordinate-based distances, to its
     * 'number_of_octant_neighbors' nearest neighbors in each of the 8 octants
     * around it (a graph with the nearest neighbor of each octant contains a
     * Euclidean minimum spanning tree).
     *
     * If 'number_of_nearest_neighbors' is negative or at least the number of
     * vertices minus one, the graph is complete (O(n^2) time per
     * iteration).
     */
    VertexId number_of_nearest_neighbors = 10;

    /** See 'number_of_nearest_neighbors'. */
    VertexId number_of_octant_neighbors = 1;

    /**
     * If the graph is complete, the ascent is done on the graph of the
     * 'number_of_ascent_candidates' alpha-nearest edges of each vertex in the
     * first 1-tree (the final 1-tree and the alpha values are computed on the
     * complete graph); '-1': on the complete graph.
     */
    VertexId number_of_ascent_candidates = 50;

    /**
     * Number of iterations of the first period of the subgradient ascent
     * ('-1': half the number of vertices, at least 100).
     */
    int64_t initial_period = -1;

    /** Initial step size of the subgradient ascent. */
    Distance initial_step_size = 1;

    /** Maximum number of iterations of the subgradient ascent ('-1': no limit). */
    int64_t maximum_number_of_iterations = -1;
};

struct AlphaNearnessOutput
{
    /** Candidates of each vertex, by increasing alpha value. */
    CandidateLists candidates;

    /** Alpha value of each candidate (in transformed cost units). */
    std::vector<std::vector<Distance>> alphas;

    /** Penalties of the vertices, from the best iteration of the ascent. */
    std::vector<Distance> pi;

    /**
     * 1-tree bound (in distance units): if the graph is complete, the bound
     * of the minimum 1-tree for the final penalties, a lower bound of the
     * length of a tour; otherwise, the largest bound of the ascent, whose
     * 1-trees are only minimum on the graph.
     */
    double lower_bound = 0;

    /** Number of iterations of the ascent. */
    int64_t number_of_iterations = 0;

    /** Whether a 1-tree of the ascent was a tour (then an optimal one). */
    bool tour_found = false;
};

/**
 * Alpha-nearness candidates (Helsgaun, 2000).
 *
 * The alpha value of an edge is the increase of the length of a minimum
 * 1-tree forced to contain it. The 1-trees are computed on transformed costs
 * 'd(i, j) + pi[i] + pi[j]', the penalties 'pi' being those which maximize
 * the 1-tree lower bound, approximated by a subgradient ascent (Held and
 * Karp, 1970-1971).
 *
 * The 1-trees are minimum spanning trees with an extra edge between a leaf
 * (the special vertex) and its second nearest neighbor, the leaf being
 * chosen to maximize the length of that edge.
 */
AlphaNearnessOutput alpha_nearness_candidates(
        const Distances& distances,
        const AlphaNearnessParameters& parameters = {});

}
