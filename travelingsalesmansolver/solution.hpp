#pragma once

#include "travelingsalesmansolver/instance.hpp"

#include <string>

namespace travelingsalesmansolver
{

/**
 * Solution class for a 'travelingsalesman' problem.
 */
class Solution
{

public:

    /*
     * Constructors and destructor.
     */

    /** Create an empty solution. */
    Solution(const Instance& instance);

    /** Create a solution from a certificate file. */
    Solution(
            const Instance& instance,
            std::string certificate_path);

    /*
     * Getters
     */

    /** Get the instance. */
    inline const Instance& instance() const { return *instance_; }

    /** Get the number of vertices in the solution. */
    inline VertexId number_of_vertices() const { return vertex_ids_.size(); }

    /** Get the full distance of the solution. */
    inline Distance distance() const { return distance_; }

    /**
     * Return 'true' iff the solution is feasible.
     * - all vertices have been visited
     */
    bool feasible() const;

    /*
     * Setters
     */

    /** Append a vertex at the end of the solution. */
    void add_vertex(VertexId vertex_id);

    /*
     * Export.
     */

    /** Write the solution to a file. */
    void write(std::string certificate_path) const;

private:

    /*
     * Private attributes.
     */

    /** Instance. */
    const Instance* instance_;

    /** List of visited vertices. */
    std::vector<VertexId> vertex_ids_;

    /**
     * Array indexed by vertices indicating whether of not they have been
     * visited.
     */
    std::vector<uint8_t> vertices_is_visited_;

    /**
     * Travelled distance.
     *
     * This time doesn't take into account the return to the depot.
     */
    Distance distance_cur_ = 0;

    /**
     * Full travelled distance.
     *
     * This time takes into account the return to the depot.
     */
    Distance distance_ = 0;

};

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////// Output ////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/**
 * Output structure for a 'travelingsalesman' problem.
 */
struct Output
{
    /** Constructor. */
    Output(
            const Instance& instance,
            optimizationtools::Info& info);

    /** Solution. */
    Solution solution;

    /** Bound. */
    Distance bound = 0;

    /** Elapsed time. */
    double time = -1;

    /** Print current state. */
    void print(
            optimizationtools::Info& info,
            const std::stringstream& s) const;

    /** Update the solution. */
    void update_solution(
            const Solution& solution_new,
            const std::stringstream& s,
            optimizationtools::Info& info);

    /** Method to call at the end of the algorithm. */
    Output& algorithm_end(optimizationtools::Info& info);
};

}

