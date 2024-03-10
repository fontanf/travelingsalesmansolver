#pragma once

#include "travelingsalesmansolver/distances/distances.hpp"

#include <iomanip>

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
    inline VertexId number_of_vertices() const { return distances_->number_of_vertices(); }

    /** Get distances. */
    const Distances& distances() const { return *distances_; }

    /*
     * Export
     */

    /** Print the instance. */
    std::ostream& format(
            std::ostream& os,
            int verbosity_level = 1) const;

    /** Print the instance. */
    template <typename Distances>
    std::ostream& format(
            const Distances& distances,
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

    /** Check a certificate. */
    template <typename Distances>
    std::pair<bool, Distance> check(
            const Distances& distances,
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

template <typename Distances>
std::ostream& Instance::format(
        const Distances& distances,
        std::ostream& os,
        int verbosity_level) const
{
    if (verbosity_level >= 1) {
        os << "Number of vertices:  " << number_of_vertices() << std::endl;
    }

    if (verbosity_level >= 2) {
        os << std::endl
            << std::setw(12) << "Loc. 1"
            << std::setw(12) << "Loc. 2"
            << std::setw(12) << "Distance"
            << std::endl
            << std::setw(12) << "------"
            << std::setw(12) << "------"
            << std::setw(12) << "--------"
            << std::endl;
        for (VertexId vertex_id_1 = 0;
                vertex_id_1 < number_of_vertices();
                ++vertex_id_1) {
            for (VertexId vertex_id_2 = vertex_id_1 + 1;
                    vertex_id_2 < number_of_vertices();
                    ++vertex_id_2) {
                os
                    << std::setw(12) << vertex_id_1
                    << std::setw(12) << vertex_id_2
                    << std::setw(12) << distances.distance(vertex_id_1, vertex_id_2)
                    << std::endl;
            }
        }
    }

    return os;
}

}
