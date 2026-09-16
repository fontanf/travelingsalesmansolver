#pragma once

#include "travelingsalesmansolver/distances/commons.hpp"

#include <vector>
#include <cmath>
#include <fstream>

namespace travelingsalesmansolver
{

/**
 * Manhattan (L1) distance between 2D coordinates ('MAN_2D' in the TSPLIB
 * format), stored as one pair of coordinates per vertex rather than an
 * explicit O(n^2) distance matrix -- unlike 'DistancesExplicit', this scales
 * to instances with tens of thousands of vertices or more, and lets LKH
 * build its candidate edge set from the coordinates (e.g.
 * 'CANDIDATE_SET_TYPE = POPMUSIC'), which an explicit matrix cannot offer
 * since it carries no geometric structure.
 *
 * Unlike 'DistancesEuc2D' and 'DistancesCeil2D', this involves no rounding:
 * the Manhattan distance between two points with integer coordinates is
 * always itself an integer, so this is an exact distance, not an
 * approximation of one.
 */
class DistancesMan2D
{

public:

    /** Get the distance between two vertices. */
    inline Distance distance(
            VertexId vertex_id_1,
            VertexId vertex_id_2) const
    {
        double xd = vertices_[vertex_id_2].x - vertices_[vertex_id_1].x;
        double yd = vertices_[vertex_id_2].y - vertices_[vertex_id_1].y;
        return (Distance)(std::abs(xd) + std::abs(yd));
    }

    /** Get the coordinates of a vertex. */
    inline const Coordinates2D& coordinates(VertexId vertex_id) const { return vertices_[vertex_id]; }

    /** Write to a file. */
    void write(std::ofstream& file) const
    {
        file << "EDGE_WEIGHT_TYPE : MAN_2D" << std::endl;
        file << "NODE_COORD_TYPE : TWOD_COORDS" << std::endl;
        file << "NODE_COORD_SECTION" << std::endl;
        for (VertexId vertex_id = 0;
                vertex_id < (VertexId)vertices_.size();
                ++vertex_id) {
            const Coordinates2D& vertex = this->vertices_[vertex_id];
            file << vertex_id + 1
                << " " << vertex.x
                << " " << vertex.y
                << std::endl;
        }
    }

private:

    /** Constructor. */
    DistancesMan2D() { }

    /** Vertices. */
    std::vector<Coordinates2D> vertices_;

    friend class DistancesMan2DBuilder;
};

class DistancesMan2DBuilder
{

public:

    /** Constructor. */
    DistancesMan2DBuilder() { }

    /** Set the number of vertices. */
    void set_number_of_vertices(
            VertexId number_of_vertices)
    {
        distances_.vertices_ = std::vector<Coordinates2D>(number_of_vertices);
    }

    /** Set the coordinates of a vertex. */
    void set_coordinates(
            VertexId vertex_id,
            double x,
            double y)
    {
        distances_.vertices_[vertex_id].x = x;
        distances_.vertices_[vertex_id].y = y;
    }

    /** Build. */
    DistancesMan2D build()
    {
        return std::move(distances_);
    }

private:

    /** Distances between vertices. */
    DistancesMan2D distances_;

};

}
