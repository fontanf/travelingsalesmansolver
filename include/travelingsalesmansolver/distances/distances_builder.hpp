#pragma once

#include "travelingsalesmansolver/distances/distances.hpp"

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

    /** Set the number of vertices. */
    void set_number_of_vertices(
            VertexId number_of_vertices)
    {
        distances_.number_of_vertices_ = number_of_vertices;
    }

    /** Set explicit distances. */
    void set_distances_explicit(
            const DistancesExplicit& distances_explicit)
    {
        if (set_) {
            throw std::invalid_argument("Distances has already been set.");
        }
        distances_.distances_explicit_ = std::make_unique<const DistancesExplicit>(distances_explicit);
        set_ = true;
    }

    /** Set explicit triangle distances. */
    void set_distances_explicit_triangle(
            const DistancesExplicitTriangle& distances_explicit_triangle)
    {
        if (set_) {
            throw std::invalid_argument("Distances has already been set.");
        }
        distances_.distances_explicit_triangle_ = std::make_unique<const DistancesExplicitTriangle>(distances_explicit_triangle);
        set_ = true;
    }

    /** Set Euc2D distances. */
    void set_distances_euc_2d(
            const DistancesEuc2D& distances_euc_2d)
    {
        if (set_) {
            throw std::invalid_argument("Distances has already been set.");
        }
        distances_.distances_euc_2d_ = std::make_unique<const DistancesEuc2D>(distances_euc_2d);
        set_ = true;
    }

    /** Set Ceil2D distances. */
    void set_distances_ceil_2d(
            const DistancesCeil2D& distances_ceil_2d)
    {
        if (set_) {
            throw std::invalid_argument("Distances has already been set.");
        }
        distances_.distances_ceil_2d_ = std::make_unique<const DistancesCeil2D>(distances_ceil_2d);
        set_ = true;
    }

    /** Set Geo distances. */
    void set_distances_geo(
            const DistancesGeo& distances_geo)
    {
        if (set_) {
            throw std::invalid_argument("Distances has already been set.");
        }
        distances_.distances_geo_ = std::make_unique<const DistancesGeo>(distances_geo);
        set_ = true;
    }

    /** Set Att distances. */
    void set_distances_att(
            const DistancesAtt& distances_att)
    {
        if (set_) {
            throw std::invalid_argument("Distances has already been set.");
        }
        distances_.distances_att_ = std::make_unique<const DistancesAtt>(distances_att);
        set_ = true;
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

    bool set_ = false;

    std::string edge_weight_type_ = "";

    std::string edge_weight_format_ = "";

    std::string node_coord_type_ = "TWOD_COORDS";

};

}
