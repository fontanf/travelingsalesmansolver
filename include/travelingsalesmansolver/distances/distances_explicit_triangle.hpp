#pragma once

#include "travelingsalesmansolver/distances/commons.hpp"

#include <vector>
#include <limits>
#include <fstream>

namespace travelingsalesmansolver
{

class DistancesExplicitTriangle
{

public:

    /** Get the distance between two vertices. */
    inline Distance distance(
            VertexId vertex_id_1,
            VertexId vertex_id_2) const
    {
        if (vertex_id_1 > vertex_id_2) {
            return distances_[vertex_id_1][vertex_id_2];
        } else {
            return distances_[vertex_id_2][vertex_id_1];
        }
    }

    /** Write to a file. */
    void write(std::ofstream& file) const
    {
        file << "EDGE_WEIGHT_TYPE : EXPLICIT" << std::endl;
        file << "EDGE_WEIGHT_FORMAT : UPPER_DIAG_ROW" << std::endl;
        file << "EDGE_WEIGHT_SECTION" << std::endl;
        for (VertexId vertex_id_1 = 0;
                vertex_id_1 < (VertexId)distances_.size();
                ++vertex_id_1) {
            for (VertexId vertex_id_2 = vertex_id_1;
                    vertex_id_2 < (VertexId)distances_.size();
                    ++vertex_id_2) {
                file << distance(vertex_id_1, vertex_id_2) << " ";
            }
        }
        file << std::endl;
    }

private:

    /** Constructor. */
    DistancesExplicitTriangle() { }

    /** Distances between vertices. */
    std::vector<std::vector<Distance>> distances_;

    friend class DistancesExplicitTriangleBuilder;
};

class DistancesExplicitTriangleBuilder
{

public:

    /** Constructor. */
    DistancesExplicitTriangleBuilder() { }

    /** Set the number of vertices. */
    void set_number_of_vertices(VertexId number_of_vertices)
    {
        distances_.distances_ = std::vector<std::vector<Distance>>(number_of_vertices);
        for (VertexId vertex_id = 0;
                vertex_id < number_of_vertices;
                ++vertex_id) {
            distances_.distances_[vertex_id] = std::vector<Distance>(
                    vertex_id + 1,
                    std::numeric_limits<Distance>::max());
        }
    }

    /** Set the distance between two vertices. */
    inline void set_distance(
            VertexId vertex_id_1,
            VertexId vertex_id_2,
            Distance distance)
    {
        if (vertex_id_1 > vertex_id_2) {
            distances_.distances_[vertex_id_1][vertex_id_2] = distance;
        } else {
            distances_.distances_[vertex_id_2][vertex_id_1] = distance;
        }
    }

    /** Build. */
    DistancesExplicitTriangle build()
    {
        return std::move(distances_);
    }

private:

    /** Distances between vertices. */
    DistancesExplicitTriangle distances_;

};

}
