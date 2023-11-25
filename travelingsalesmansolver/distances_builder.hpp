#pragma once

#include "travelingsalesmansolver/distances.hpp"

namespace travelingsalesmansolver
{

class DistancesBuilder
{

public:

    /** Constructor. */
    DistancesBuilder() { }

    /** Read an distances from a file in 'tsplib' format. */
    bool read_tsplib(
            std::ifstream& file,
            const std::string& tmp,
            const std::vector<std::string>& line);

    /** Add vertices. */
    void add_vertices(VertexId number_of_vertices);

    /** Set the distance between two vertices. */
    inline void set_distance(
            VertexId vertex_id_1,
            VertexId vertex_id_2,
            Distance distance)
    {
        distances_.set_distance(vertex_id_1, vertex_id_2, distance);
    }

    /** Set the coordinates of a vertex. */
    void set_coordinates(
            VertexId vertex_id,
            double x,
            double y,
            double z = -1)
    {
        distances_.vertices_[vertex_id].x = x;
        distances_.vertices_[vertex_id].y = y;
        distances_.vertices_[vertex_id].z = z;
    }

    /** Set the edge weight type. */
    void set_edge_weight_type(std::string edge_weight_type)
    {
        distances_.edge_weight_type_ = edge_weight_type;
    }

    /** Set the node coord type. */
    void set_node_coord_type(std::string node_coord_type)
    {
        distances_.node_coord_type_ = node_coord_type;
    }

    /*
     * Build
     */

    /** Build. */
    Distances build();

private:

    /*
     * Private attributes
     */

    /** Distances. */
    Distances distances_;

};

}
