#pragma once

#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <limits>
#include <cmath>

namespace travelingsalesmansolver
{

using VertexId = int64_t;
using Distance = int64_t;

/*
 * Structure for a vertex.
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
 * Class that stores distance information.
 */
class Distances
{

public:

    /*
     * Constructors and destructor
     */

    void compute_distances() const;

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
        if (!distances_.empty()) {
            if (vertex_id_1 > vertex_id_2) {
                return distances_[vertex_id_1][vertex_id_2];
            } else {
                return distances_[vertex_id_2][vertex_id_1];
            }
        } else if (edge_weight_type_ == "EUC_2D") {
            double xd = x(vertex_id_2) - x(vertex_id_1);
            double yd = y(vertex_id_2) - y(vertex_id_1);
            return std::round(std::sqrt(xd * xd + yd * yd));
        } else if (edge_weight_type_ == "CEIL_2D") {
            double xd = x(vertex_id_2) - x(vertex_id_1);
            double yd = y(vertex_id_2) - y(vertex_id_1);
            return std::ceil(std::sqrt(xd * xd + yd * yd));
        } else if (edge_weight_type_ == "GEO") {
            double rrr = 6378.388;
            double q1 = cos(longitudes_[vertex_id_1] - longitudes_[vertex_id_2]);
            double q2 = cos(latitudes_[vertex_id_1] - latitudes_[vertex_id_2]);
            double q3 = cos(latitudes_[vertex_id_1] + latitudes_[vertex_id_2]);
            return (Distance)(rrr * acos(0.5 * ((1.0 + q1) * q2 - (1.0 - q1) * q3)) + 1.0);
        } else if (edge_weight_type_ == "ATT") {
            double xd = x(vertex_id_1) - x(vertex_id_2);
            double yd = y(vertex_id_1) - y(vertex_id_2);
            double rij = sqrt((xd * xd + yd * yd) / 10.0);
            int tij = std::round(rij);
            return (tij < rij)? tij + 1: tij;
        } else {
            throw std::invalid_argument(
                "Unknown edge weight type \"" + edge_weight_type_ + "\".");
            return -1;
        }
    }

    /*
     * Export
     */

    /** Print the instance. */
    std::ostream& print(
            std::ostream& os,
            int verbose = 1) const;

    /** Write the instance to a file. */
    void write(std::ofstream& file) const;

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
    Distances() { }

    inline void init_distances() const
    {
        distances_ = std::vector<std::vector<Distance>>(number_of_vertices());
        for (VertexId vertex_id = 0;
                vertex_id < number_of_vertices();
                ++vertex_id) {
            distances_[vertex_id] = std::vector<Distance>(
                    vertex_id,
                    std::numeric_limits<Distance>::max());
        }
    }

    /** Set the distance between two vertices. */
    inline void set_distance(
            VertexId vertex_id_1,
            VertexId vertex_id_2,
            Distance distance) const
    {
        if (vertex_id_1 > vertex_id_2) {
            distances_[vertex_id_1][vertex_id_2] = distance;
        } else {
            distances_[vertex_id_2][vertex_id_1] = distance;
        }
    }

    /*
     * Private attributes
     */

    /** Vertices. */
    std::vector<Vertex> vertices_;

    /** Distances between vertices. */
    mutable std::vector<std::vector<Distance>> distances_;

    /*
     * Computed attributes
     */

    std::string edge_weight_format_ = "";

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

    /** Structure for GEO edge weight type. */
    mutable std::vector<double> latitudes_;

    /** Structure for GEO edge weight type. */
    mutable std::vector<double> longitudes_;

    friend class DistancesBuilder;

};

}

