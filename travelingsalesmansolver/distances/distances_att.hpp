#pragma once

#include "travelingsalesmansolver/distances/commons.hpp"

#include <vector>
#include <cmath>
#include <fstream>

namespace travelingsalesmansolver
{

class DistancesAtt
{

public:

    /** Get the distance between two vertices. */
    inline Distance distance(
            VertexId vertex_id_1,
            VertexId vertex_id_2) const
    {
        double xd = vertices_[vertex_id_1].x - vertices_[vertex_id_2].x;
        double yd = vertices_[vertex_id_1].y - vertices_[vertex_id_2].y;
        double rij = sqrt((xd * xd + yd * yd) / 10.0);
        int tij = std::round(rij);
        return (tij < rij)? tij + 1: tij;
    }

    /** Write to a file. */
    void write(std::ofstream& file) const
    {
        file << "EDGE_WEIGHT_TYPE : EUC_2D" << std::endl;
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
    DistancesAtt() { }

    /** Vertices. */
    std::vector<Coordinates2D> vertices_;

    friend class DistancesAttBuilder;
};

class DistancesAttBuilder
{

public:

    /** Constructor. */
    DistancesAttBuilder() { }

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
    DistancesAtt build()
    {
        return std::move(distances_);
    }

private:

    /** Distances between vertices. */
    DistancesAtt distances_;

};

}
