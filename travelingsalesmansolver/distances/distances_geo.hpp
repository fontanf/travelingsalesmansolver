#pragma once

#include "travelingsalesmansolver/distances/commons.hpp"

#include <vector>
#include <cmath>
#include <fstream>

namespace travelingsalesmansolver
{

class DistancesGeo
{

public:

    struct Vertex
    {
        /** Coordinates. */
        Coordinates2D coordinates;

        /** Latitude. */
        double latitude;

        /** Longitude. */
        double longitude;
    };

    /** Get the distance between two vertices. */
    inline Distance distance(
            VertexId vertex_id_1,
            VertexId vertex_id_2) const
    {
        double rrr = 6378.388;
        double q1 = cos(vertices_[vertex_id_1].longitude - vertices_[vertex_id_2].longitude);
        double q2 = cos(vertices_[vertex_id_1].latitude - vertices_[vertex_id_2].latitude);
        double q3 = cos(vertices_[vertex_id_1].latitude + vertices_[vertex_id_2].latitude);
        return (Distance)(rrr * acos(0.5 * ((1.0 + q1) * q2 - (1.0 - q1) * q3)) + 1.0);
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
            const Vertex& vertex = this->vertices_[vertex_id];
            file << vertex_id + 1
                << " " << vertex.coordinates.x
                << " " << vertex.coordinates.y
                << std::endl;
        }
    }

private:

    /** Constructor. */
    DistancesGeo() { }

    /** Vertices. */
    std::vector<Vertex> vertices_;

    friend class DistancesGeoBuilder;
};

class DistancesGeoBuilder
{

public:

    /** Constructor. */
    DistancesGeoBuilder() { }

    /** Set the number of vertices. */
    void set_number_of_vertices(
            VertexId number_of_vertices)
    {
        distances_.vertices_ = std::vector<DistancesGeo::Vertex>(number_of_vertices);
    }

    /** Set the coordinates of a vertex. */
    void set_coordinates(
            VertexId vertex_id,
            double x,
            double y)
    {
        double pi = 3.141592;
        int deg_x = std::round(x);
        double min_x = x - deg_x;
        double latitude = pi * (deg_x + 5.0 * min_x / 3.0) / 180.0;
        int deg_y = std::round(y);
        double min_y = y - deg_y;
        double longitude = pi * (deg_y + 5.0 * min_y / 3.0) / 180.0;
        distances_.vertices_[vertex_id].coordinates.x = x;
        distances_.vertices_[vertex_id].coordinates.y = y;
        distances_.vertices_[vertex_id].latitude = latitude;
        distances_.vertices_[vertex_id].longitude = longitude;
    }

    /** Build. */
    DistancesGeo build()
    {
        return std::move(distances_);
    }

private:

    /** Distances between vertices. */
    DistancesGeo distances_;

};

}
