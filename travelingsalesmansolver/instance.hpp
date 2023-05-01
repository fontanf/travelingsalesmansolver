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
struct Vertex
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
        vertices_(number_of_vertices),
        distances_(number_of_vertices, std::vector<Distance>(number_of_vertices, -1))
    {
        for (VertexId vertex_id = 0; vertex_id < number_of_vertices; ++vertex_id)
            distances_[vertex_id][vertex_id] = 0;
    }

    /** Set the distance between two vertices. */
    inline void set_distance(
            VertexId vertex_id_1,
            VertexId vertex_id_2,
            Distance distance)
    {
        distances_[vertex_id_1][vertex_id_2] = distance;
        distances_[vertex_id_2][vertex_id_1] = distance;
        distance_max_ = std::max(distance_max_, distance);
    }

    /** Set the coordinates of a vertex. */
    void set_coordinates(
            VertexId vertex_id,
            double x,
            double y,
            double z = -1)
    {
        vertices_[vertex_id].x = x;
        vertices_[vertex_id].y = y;
        vertices_[vertex_id].z = z;
    }

    /** Set the edge weight type. */
    void set_edge_weight_type(std::string edge_weight_type)
    {
        edge_weight_type_ = edge_weight_type;
    }

    /** Set the node coord type. */
    void set_node_coord_type(std::string node_coord_type)
    {
        node_coord_type_ = node_coord_type;
    }

    /** Build an instance from a file. */
    Instance(
            std::string instance_path,
            std::string format = "");

    /*
     * Getters
     */

    /** Get the number of vertices. */
    inline VertexId number_of_vertices() const { return vertices_.size(); }

    /** Get a vertex. */
    inline const Vertex& vertex(VertexId vertex_id) const { return vertices_[vertex_id]; }

    /** Get the x-coordinate of a vertex. */
    inline double x(VertexId vertex_id) const { return vertices_[vertex_id].x; }

    /** Get the y-coordinate of a vertex. */
    inline double y(VertexId vertex_id) const { return vertices_[vertex_id].y; }

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

    /** Vertices. */
    std::vector<Vertex> vertices_;

    /** Distances between vertices. */
    std::vector<std::vector<Distance>> distances_;

    /*
     * Computed attributes
     */

    /** Maximum distance. */
    Distance distance_max_ = 0;

    /**
     * Edge weight type.
     *
     * See http://comopt.ifi.uni-heidelberg.de/software/TSPLIB95/tsp95.pdf
     */
    std::string edge_weight_type_ = "EXPLICIT";

    /**
     * Node coord type.
     *
     * See http://comopt.ifi.uni-heidelberg.de/software/TSPLIB95/tsp95.pdf
     */
    std::string node_coord_type_ = "TWOD_COORDS";

};

void init_display(
        const Instance& instance,
        optimizationtools::Info& info);

}

