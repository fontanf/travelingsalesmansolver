#pragma once

#include "travelingsalesmansolver/distances.hpp"

#include "optimizationtools/containers/indexed_set.hpp"
#include "optimizationtools/utils/utils.hpp"
#include "optimizationtools/utils/info.hpp"

#include <fstream>
#include <iostream>
#include <iomanip>

namespace travelingsalesmansolver
{

using VertexId = int64_t;
using VertexPos = int64_t;
using Distance = int64_t;
using Counter = int64_t;

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
    Instance(const std::shared_ptr<const Distances>& distances): distances_(distances) { }

    Instance(
            std::string instance_path,
            std::string format = "");

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
    std::ostream& print(
            std::ostream& os,
            int verbose = 1) const;

    /** Write the instance to a file. */
    void write(std::string instance_path) const;

    /** Check a certificate. */
    std::pair<bool, Distance> check(
            std::string certificate_path,
            std::ostream& os,
            int verbose = 1) const;

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

void init_display(
        const Instance& instance,
        optimizationtools::Info& info);

}

