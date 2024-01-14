#pragma once

#include "travelingsalesmansolver/distances.hpp"

#include <memory>

namespace travelingsalesmansolver
{

using VertexId = int64_t;
using VertexPos = int64_t;
using Distance = int64_t;
using Counter = int64_t;
using Seed = int64_t;

/**
 * Instance class for a 'travelingsalesman' problem.
 */
class Instance
{

public:

    /*
     * Constructors and destructor
     */

    /** Build an instance from a distances object. */
    Instance(
            const std::shared_ptr<const Distances>& distances):
        distances_(distances) { }

    Instance(
            const std::string& instance_path,
            const std::string& format = "");

    /*
     * Getters
     */

    /** Get the number of vertices. */
    inline VertexId number_of_vertices() const { return distances().number_of_vertices(); }

    /** Get distances. */
    const Distances& distances() const { return *distances_; }

    /*
     * Export
     */

    /** Print the instance. */
    std::ostream& format(
            std::ostream& os,
            int verbosity_level = 1) const;

    /** Write the instance to a file. */
    void write(
            const std::string& instance_path) const;

    /** Check a certificate. */
    std::pair<bool, Distance> check(
            const std::string& certificate_path,
            std::ostream& os,
            int verbosity_level = 1) const;

private:

    /*
     * Private methods
     */

    /** Create an instance manually. */
    Instance() { }

    /*
     * Read input file
     */

    /** Read an instance from a file in 'polyakovskiy2014' format. */
    void read_tsplib(std::ifstream& file);

    /*
     * Private attributes
     */

    /** Distances. */
    std::shared_ptr<const Distances> distances_;

    friend class InstanceBuilder;

};

}
