#pragma once

#include "travelingsalesmansolver/distances/commons.hpp"

#include "nlohmann/json.hpp"

#include <vector>
#include <cmath>
#include <fstream>

namespace travelingsalesmansolver
{

/**
 * Euclidean distances, without the TSPLIB 'EUC_2D' rounding-to-nearest-integer
 * convention.
 *
 * The distance between two vertices is a fixed-point representation of the
 * raw Euclidean distance, scaled by 'scale()' (and rounded to the nearest
 * integer only at that much finer resolution). This keeps the shared
 * 'Distance' ('int64_t') contract used throughout this library, while
 * matching literature/benchmarks that report raw (unrounded) Euclidean
 * distances. Callers that need the human-readable distance should divide
 * the result by 'scale()'.
 */
class DistancesEuc2DUnrounded
{

public:

    /** Scale factor used to represent the raw distance as an integer. */
    static constexpr Distance scale() { return 1000000; }

    /** Get the distance between two vertices, scaled by 'scale()'. */
    inline Distance distance(
            VertexId vertex_id_1,
            VertexId vertex_id_2) const
    {
        double xd = vertices_[vertex_id_2].x - vertices_[vertex_id_1].x;
        double yd = vertices_[vertex_id_2].y - vertices_[vertex_id_1].y;
        return (Distance)std::llround(std::sqrt(xd * xd + yd * yd) * (double)scale());
    }

    /** Get the coordinates of a vertex. */
    inline const Coordinates2D& coordinates(VertexId vertex_id) const { return vertices_[vertex_id]; }

    /** Write to a file. */
    void write(std::ofstream& file) const
    {
        file << "EDGE_WEIGHT_TYPE : EUC_2D_UNROUNDED" << std::endl;
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

    /** Export to a JSON structure. */
    nlohmann::json to_json() const
    {
        nlohmann::json j;
        j["type"] = "euc_2d_unrounded";
        j["number_of_vertices"] = (VertexId)vertices_.size();
        for (VertexId vertex_id = 0; vertex_id < (VertexId)vertices_.size(); ++vertex_id)
            j["coordinates"][vertex_id] = {vertices_[vertex_id].x, vertices_[vertex_id].y};
        return j;
    }

private:

    /** Constructor. */
    DistancesEuc2DUnrounded() { }

    /** Vertices. */
    std::vector<Coordinates2D> vertices_;

    friend class DistancesEuc2DUnroundedBuilder;
};

class DistancesEuc2DUnroundedBuilder
{

public:

    /** Constructor. */
    DistancesEuc2DUnroundedBuilder() { }

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
    DistancesEuc2DUnrounded build()
    {
        return std::move(distances_);
    }

private:

    /** Distances. */
    DistancesEuc2DUnrounded distances_;

};

}
