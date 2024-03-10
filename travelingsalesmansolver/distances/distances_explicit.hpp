#pragma once

#include "travelingsalesmansolver/distances/commons.hpp"

#include <vector>
#include <limits>
#include <fstream>

namespace travelingsalesmansolver
{

class DistancesExplicit
{

public:

    /** Get the distance between two vertices. */
    inline Distance distance(
            VertexId vertex_id_1,
            VertexId vertex_id_2) const
    {
        return distances_[vertex_id_1][vertex_id_2];
    }

    /** Write to a file. */
    void write(std::ofstream& file) const
    {
        file << "EDGE_WEIGHT_TYPE : EXPLICIT" << std::endl;
        file << "EDGE_WEIGHT_FORMAT : FULL_MATRIX" << std::endl;
        file << "EDGE_WEIGHT_SECTION" << std::endl;
        for (VertexId vertex_id_1 = 0;
                vertex_id_1 < (VertexId)distances_.size();
                ++vertex_id_1) {
            for (VertexId vertex_id_2 = 0;
                    vertex_id_2 < (VertexId)distances_.size();
                    ++vertex_id_2) {
                file << distance(vertex_id_1, vertex_id_2) << " ";
            }
        }
        file << std::endl;
    }

private:

    /** Constructor. */
    DistancesExplicit() { }

    /** Distances between vertices. */
    std::vector<std::vector<Distance>> distances_;

    friend class DistancesExplicitBuilder;
};

class DistancesExplicitBuilder
{

public:

    /** Constructor. */
    DistancesExplicitBuilder() { }

    /** Set the number of vertices. */
    void set_number_of_vertices(
            VertexId number_of_vertices)
    {
        distances_.distances_ = std::vector<std::vector<Distance>>(number_of_vertices);
        for (VertexId vertex_id = 0;
                vertex_id < number_of_vertices;
                ++vertex_id) {
            distances_.distances_[vertex_id] = std::vector<Distance>(
                    number_of_vertices,
                    std::numeric_limits<Distance>::max());
        }
    }

    /** Set the distance between two vertices. */
    inline void set_distance(
            VertexId vertex_id_1,
            VertexId vertex_id_2,
            Distance distance)
    {
        distances_.distances_[vertex_id_1][vertex_id_2] = distance;
        distances_.distances_[vertex_id_2][vertex_id_1] = distance;
    }

    /** Build. */
    DistancesExplicit build()
    {
        return std::move(distances_);
    }

private:

    /** Distances between vertices. */
    DistancesExplicit distances_;

};

}
