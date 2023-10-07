#pragma once

#include "travelingsalesmansolver/instance.hpp"

namespace travelingsalesmansolver
{

class InstanceBuilder
{

public:

    /** Constructor. */
    InstanceBuilder() { }

    /** Read an instance from a file. */
    InstanceBuilder(
            std::string instance_path,
            std::string format = "");

    /** Add vertices. */
    void add_vertices(VertexId number_of_vertices);

    /** Set the distance between two vertices. */
    inline void set_distance(
            VertexId vertex_id_1,
            VertexId vertex_id_2,
            Distance distance)
    {
        instance_.set_distance(vertex_id_1, vertex_id_2, distance);
    }

    /** Set the coordinates of a vertex. */
    void set_coordinates(
            VertexId vertex_id,
            double x,
            double y,
            double z = -1)
    {
        instance_.vertices_[vertex_id].x = x;
        instance_.vertices_[vertex_id].y = y;
        instance_.vertices_[vertex_id].z = z;
    }

    /** Set the edge weight type. */
    void set_edge_weight_type(std::string edge_weight_type)
    {
        instance_.edge_weight_type_ = edge_weight_type;
    }

    /** Set the node coord type. */
    void set_node_coord_type(std::string node_coord_type)
    {
        instance_.node_coord_type_ = node_coord_type;
    }

    /*
     * Build
     */

    /** Build. */
    Instance build();

private:

    /*
     * Private methods
     */

    /*
     * Read input file
     */

    /** Read an instance from a file in 'tsplib' format. */
    void read_tsplib(std::ifstream& file);

    /*
     * Private attributes
     */

    /** Instance. */
    Instance instance_;

};

}
