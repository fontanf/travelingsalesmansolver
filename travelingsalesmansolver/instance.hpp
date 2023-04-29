#pragma once

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

/*
 * Structure for a location.
 */
struct Location
{
    /** x-coordinate of the location. */
    double x;

    /** y-coordinate of the location. */
    double y;

    /** z-coordinate of the location. */
    double z;
};

/**
 * Instance class for a 'travelingsalesman' problem.
 */
class Instance
{

public:

    /*
     * Constructors and destructor
     */

    /** Constructor to build an instance manually. */
    Instance(VertexId number_of_vertices):
        locations_(number_of_vertices),
        distances_(number_of_vertices, std::vector<Distance>(number_of_vertices, -1))
    {
        for (VertexId vertex_id = 0; vertex_id < number_of_vertices; ++vertex_id)
            distances_[vertex_id][vertex_id] = 0;
    }

    /** Set the coordinates of a vertex. */
    void set_xy(
            VertexId vertex_id,
            double x,
            double y,
            double z = -1)
    {
        locations_[vertex_id].x = x;
        locations_[vertex_id].y = y;
        locations_[vertex_id].z = z;
    }

    /** Set the distance between two vertices. */
    void set_distance(
            VertexId vertex_id_1,
            VertexId vertex_id_2,
            Distance distance)
    {
        distances_[vertex_id_1][vertex_id_2] = distance;
        distances_[vertex_id_2][vertex_id_1] = distance;
        distance_max_ = std::max(distance_max_, distance);
    }

    /** Build an instance from a file. */
    Instance(
            std::string instance_path,
            std::string format = "");

    /*
     * Getters
     */

    /** Get the number of vertices. */
    inline VertexId number_of_vertices() const { return locations_.size(); }

    /** Get the x-coordinate of a vertex. */
    inline double x(VertexId vertex_id) const { return locations_[vertex_id].x; }

    /** Get the y-coordinate of a vertex. */
    inline double y(VertexId vertex_id) const { return locations_[vertex_id].y; }

    /** Get the distance between two vertices. */
    inline Distance distance(
            VertexId vertex_id_1,
            VertexId vertex_id_2) const
    {
        return distances_[vertex_id_1][vertex_id_2];
    }

    /** Get the maximum distance between two vertices. */
    inline Distance maximum_distance() const { return distance_max_; }

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

    /** Read an instance from a file in 'tsplib' format. */
    void read_tsplib(std::ifstream& file);

    /*
     * Private attributes
     */

    /** Locations. */
    std::vector<Location> locations_;

    /** Distances between locations. */
    std::vector<std::vector<Distance>> distances_;

    /*
     * Computed attributes
     */

    /** Maximum distance. */
    Distance distance_max_ = 0;

};

void init_display(
        const Instance& instance,
        optimizationtools::Info& info);

}

