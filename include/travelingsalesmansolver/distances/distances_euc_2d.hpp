#pragma once

#include "travelingsalesmansolver/distances/commons.hpp"

#include "nlohmann/json.hpp"

#include <vector>
#include <cmath>
#include <fstream>

namespace travelingsalesmansolver
{

/**
 * Euclidean distances, multiplied by 'scale()' and rounded to the nearest
 * integer.
 *
 * With the default scale of 1, this is TSPLIB's 'EUC_2D' edge weight type.
 * A larger scale (LKH's 'SCALE' keyword) keeps more precision, e.g. a scale
 * of '1e6' matches literature/benchmarks that report raw (unrounded)
 * Euclidean distances, while keeping the shared 'Distance' ('int64_t')
 * contract used throughout this library. Callers that need the
 * human-readable distance should divide the result by 'scale()'.
 */
class DistancesEuc2D
{

public:

    /** Get the distance between two vertices, scaled by 'scale()'. */
    inline Distance distance(
            VertexId vertex_id_1,
            VertexId vertex_id_2) const
    {
        double xd = vertices_[vertex_id_2].x - vertices_[vertex_id_1].x;
        double yd = vertices_[vertex_id_2].y - vertices_[vertex_id_1].y;
        return (Distance)std::llround((double)scale_ * std::sqrt(xd * xd + yd * yd));
    }

    /** Get the scale factor applied to the Euclidean distances before rounding them. */
    inline Distance scale() const { return scale_; }

    /** Get the coordinates of a vertex. */
    inline const Coordinates2D& coordinates(VertexId vertex_id) const { return vertices_[vertex_id]; }

    /** Write to a file. */
    void write(std::ofstream& file) const
    {
        file << "EDGE_WEIGHT_TYPE : EUC_2D" << std::endl;
        if (scale_ != 1)
            file << "SCALE : " << scale_ << std::endl;
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
        j["type"] = "euc_2d";
        j["number_of_vertices"] = (VertexId)vertices_.size();
        j["scale"] = scale_;
        for (VertexId vertex_id = 0; vertex_id < (VertexId)vertices_.size(); ++vertex_id)
            j["coordinates"][vertex_id] = {vertices_[vertex_id].x, vertices_[vertex_id].y};
        return j;
    }

private:

    /** Constructor. */
    DistancesEuc2D() { }

    /** Vertices. */
    std::vector<Coordinates2D> vertices_;

    /** Scale factor applied to the Euclidean distances before rounding them. */
    Distance scale_ = 1;

    friend class DistancesEuc2DBuilder;
};

class DistancesEuc2DBuilder
{

public:

    /** Constructor. */
    DistancesEuc2DBuilder() { }

    /** Set the number of vertices. */
    void set_number_of_vertices(
            VertexId number_of_vertices)
    {
        distances_.vertices_ = std::vector<Coordinates2D>(number_of_vertices);
    }

    /** Set the scale factor applied to the Euclidean distances before rounding them. */
    void set_scale(Distance scale)
    {
        if (scale < 1) {
            throw std::invalid_argument(
                    "DistancesEuc2DBuilder::set_scale: "
                    "the scale must be positive; "
                    "scale: " + std::to_string(scale) + ".");
        }
        distances_.scale_ = scale;
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
    DistancesEuc2D build()
    {
        return std::move(distances_);
    }

private:

    /** Distances between vertices. */
    DistancesEuc2D distances_;

};

}
