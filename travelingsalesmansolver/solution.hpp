#pragma once

#include "travelingsalesmansolver/instance.hpp"

#include "optimizationtools/utils/output.hpp"
#include "optimizationtools/utils/utils.hpp"

#include "nlohmann/json.hpp"

#include <iomanip>
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
     * Constructors and destructor
     */

    /** Create an empty solution. */
    Solution(const Instance& instance);

    /** Create a solution from a certificate file. */
    Solution(
            const Instance& instance,
            const std::string& certificate_path);

    /*
     * Getters
     */

    /** Get the instance. */
    inline const Instance& instance() const { return *instance_; }

    /** Get the number of vertices in the solution. */
    inline VertexId number_of_vertices() const { return vertex_ids_.size(); }

    /** Get the id of the vertex at a given position. */
    inline VertexId vertex_id(VertexPos pos) const { return vertex_ids_[pos]; }

    /** Get the full distance of the solution. */
    inline Distance distance() const { return distance_; }

    /** Get the total cost of the solution. */
    inline Distance objective_value() const { return distance(); }

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
     * Export
     */

    /** Write the solution to a file. */
    void write(
            const std::string& certificate_path) const;

    /** Export solution characteristics to a JSON structure. */
    nlohmann::json to_json() const;

    /** Write a formatted output of the instance to a stream. */
    void format(
            std::ostream& os,
            int verbosity_level = 1) const;

private:

    /*
     * Private attributes
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

inline optimizationtools::ObjectiveDirection objective_direction()
{
    return optimizationtools::ObjectiveDirection::Minimize;
}

/**
 * Output structure for a 'travelingsalesman' problem.
 */
struct Output: optimizationtools::Output
{
    /** Constructor. */
    Output(const Instance& instance): solution(instance) { }


    /** Solution. */
    Solution solution;

    /** Bound. */
    Distance bound = 0;

    /** Elapsed time. */
    double time = -1;


    std::string solution_value() const
    {
        return optimizationtools::solution_value(
            objective_direction(),
            solution.feasible(),
            solution.objective_value());
    }

    double absolute_optimality_gap() const
    {
        return optimizationtools::absolute_optimality_gap(
                objective_direction(),
                solution.feasible(),
                solution.objective_value(),
                bound);
    }

    double relative_optimality_gap() const
    {
       return optimizationtools::relative_optimality_gap(
            objective_direction(),
            solution.feasible(),
            solution.objective_value(),
            bound);
    }

    virtual nlohmann::json to_json() const
    {
        return nlohmann::json {
            {"Solution", solution.to_json()},
            {"Value", solution_value()},
            {"Bound", bound},
            {"AbsoluteOptimalityGap", absolute_optimality_gap()},
            {"RelativeOptimalityGap", relative_optimality_gap()},
            {"Time", time}
        };
    }

    virtual int format_width() const { return 30; }

    virtual void format(std::ostream& os) const
    {
        int width = format_width();
        os
            << std::setw(width) << std::left << "Value: " << solution_value() << std::endl
            << std::setw(width) << std::left << "Bound: " << bound << std::endl
            << std::setw(width) << std::left << "Absolute optimality gap: " << absolute_optimality_gap() << std::endl
            << std::setw(width) << std::left << "Relative optimality gap (%): " << relative_optimality_gap() * 100 << std::endl
            << std::setw(width) << std::left << "Time (s): " << time << std::endl
            ;
    }
};

using NewSolutionCallback = std::function<void(const Output&, const std::string&)>;

struct Parameters: optimizationtools::Parameters
{
    /** Callback function called when a new best solution is found. */
    NewSolutionCallback new_solution_callback = [](const Output&, const std::string&) { };


    virtual nlohmann::json to_json() const override
    {
        nlohmann::json json = optimizationtools::Parameters::to_json();
        json.merge_patch({});
        return json;
    }

    virtual int format_width() const override { return 23; }

    virtual void format(std::ostream& os) const override
    {
        optimizationtools::Parameters::format(os);
        //int width = format_width();
    }
};

}

