#pragma once

#include "travelingsalesmansolver/distances/distances_explicit.hpp"
#include "travelingsalesmansolver/distances/distances_explicit_triangle.hpp"
#include "travelingsalesmansolver/distances/distances_euc_2d.hpp"
#include "travelingsalesmansolver/distances/distances_ceil_2d.hpp"
#include "travelingsalesmansolver/distances/distances_geo.hpp"
#include "travelingsalesmansolver/distances/distances_att.hpp"

#include <fstream>
#include <iostream>
#include <memory>
#include <iomanip>

namespace travelingsalesmansolver
{

#define FUNCTION_WITH_DISTANCES_0(function, distances) \
    ((distances).distances_explicit() != nullptr)? \
    (function)(*(distances).distances_explicit()): \
    ((distances).distances_explicit_triangle() != nullptr)? \
    (function)(*(distances).distances_explicit_triangle()): \
    ((distances).distances_euc_2d() != nullptr)? \
    (function)(*(distances).distances_euc_2d()): \
    ((distances).distances_ceil_2d() != nullptr)? \
    (function)(*(distances).distances_ceil_2d()): \
    ((distances).distances_geo() != nullptr)? \
    (function)(*(distances).distances_geo()): \
    (function)(*(distances).distances_att());

#define FUNCTION_WITH_DISTANCES(function, distances, ...) \
    ((distances).distances_explicit() != nullptr)? \
    (function)(*(distances).distances_explicit(), __VA_ARGS__): \
    ((distances).distances_explicit_triangle() != nullptr)? \
    (function)(*(distances).distances_explicit_triangle(), __VA_ARGS__): \
    ((distances).distances_euc_2d() != nullptr)? \
    (function)(*(distances).distances_euc_2d(), __VA_ARGS__): \
    ((distances).distances_ceil_2d() != nullptr)? \
    (function)(*(distances).distances_ceil_2d(), __VA_ARGS__): \
    ((distances).distances_geo() != nullptr)? \
    (function)(*(distances).distances_geo(), __VA_ARGS__): \
    (function)(*(distances).distances_att(), __VA_ARGS__);

#define FUNCTION_WITH_DISTANCES_0_R(function, distances) \
    ((distances).distances_euc_2d() != nullptr)? \
    (function)(*(distances).distances_euc_2d()): \
    ((distances).distances_ceil_2d() != nullptr)? \
    (function)(*(distances).distances_ceil_2d()): \
    ((distances).distances_geo() != nullptr)? \
    (function)(*(distances).distances_geo()): \
    ((distances).distances_att() != nullptr)? \
    (function)(*(distances).distances_att()): \
    ((distances).distances_explicit_triangle() != nullptr)? \
    (function)(*(distances).distances_explicit_triangle()): \
    (function)(*(distances).distances_explicit());

#define FUNCTION_WITH_DISTANCES_R(function, distances, ...) \
    ((distances).distances_euc_2d() != nullptr)? \
    (function)(*(distances).distances_euc_2d(), __VA_ARGS__): \
    ((distances).distances_ceil_2d() != nullptr)? \
    (function)(*(distances).distances_ceil_2d(), __VA_ARGS__): \
    ((distances).distances_geo() != nullptr)? \
    (function)(*(distances).distances_geo(), __VA_ARGS__): \
    ((distances).distances_att() != nullptr)? \
    (function)(*(distances).distances_att(), __VA_ARGS__): \
    ((distances).distances_explicit_triangle() != nullptr)? \
    (function)(*(distances).distances_explicit_triangle(), __VA_ARGS__): \
    (function)(*(distances).distances_explicit(), __VA_ARGS__);

/**
 * Class that stores distance information.
 */
class Distances
{

public:

    /*
     * Constructors and destructor
     */

    void compute_distances_explicit() const;

    template <typename T>
    void compute_distances_explicit(
            const T& distances) const;

    void compute_distances_explicit_triangle() const;

    template <typename T>
    void compute_distances_explicit_triangle(
            const T& distances) const;

    /*
     * Getters
     */

    /** Get the number of vertices. */
    VertexId number_of_vertices() const { return number_of_vertices_; }

    const DistancesExplicit* distances_explicit() const { return distances_explicit_.get(); }

    const DistancesExplicitTriangle* distances_explicit_triangle() const { return distances_explicit_triangle_.get(); }

    const DistancesEuc2D* distances_euc_2d() const { return distances_euc_2d_.get(); }

    const DistancesCeil2D* distances_ceil_2d() const { return distances_ceil_2d_.get(); }

    const DistancesGeo* distances_geo() const { return distances_geo_.get(); }

    const DistancesAtt* distances_att() const { return distances_att_.get(); }

    /*
     * Export
     */

    /** Print the instance. */
    std::ostream& format(
            std::ostream& os,
            int verbosity_level = 1) const;

    /** Print the instance. */
    template <typename T>
    std::ostream& format(
            const T& distances,
            std::ostream& os,
            int verbosity_level = 1) const;

    /** Write the instance to a file. */
    void write(std::ofstream& file) const;

    /** Write the instance to a file. */
    template <typename T>
    void write(
            const T& distances,
            std::ofstream& file) const;

private:

    /*
     * Private methods
     */

    /** Create an instance manually. */
    Distances() { }

    /*
     * Private attributes
     */

    VertexId number_of_vertices_;

    mutable std::unique_ptr<const DistancesExplicit> distances_explicit_ = nullptr;

    mutable std::unique_ptr<const DistancesExplicitTriangle> distances_explicit_triangle_ = nullptr;

    std::unique_ptr<const DistancesEuc2D> distances_euc_2d_ = nullptr;

    std::unique_ptr<const DistancesCeil2D> distances_ceil_2d_ = nullptr;

    std::unique_ptr<const DistancesGeo> distances_geo_ = nullptr;

    std::unique_ptr<const DistancesAtt> distances_att_ = nullptr;

    friend class DistancesBuilder;

};

template <typename T>
void Distances::write(
        const T& distances,
        std::ofstream& file) const
{
    distances.write(file);
}

template <typename T>
void Distances::compute_distances_explicit(
        const T& distances) const
{
    DistancesExplicitBuilder distances_explicit_builder;
    distances_explicit_builder.set_number_of_vertices(number_of_vertices());
    for (VertexId vertex_id_1 = 0;
            vertex_id_1 < number_of_vertices();
            ++vertex_id_1) {
        for (VertexId vertex_id_2 = 0;
                vertex_id_2 < number_of_vertices();
                ++vertex_id_2) {
            distances_explicit_builder.set_distance(
                    vertex_id_1,
                    vertex_id_2,
                    distances.distance(vertex_id_1, vertex_id_2));
        }
    }
    distances_explicit_ = std::make_unique<const DistancesExplicit>(distances_explicit_builder.build());
}

template <typename T>
void Distances::compute_distances_explicit_triangle(
        const T& distances) const
{
    DistancesExplicitTriangleBuilder distances_explicit_triangle_builder;
    distances_explicit_triangle_builder.set_number_of_vertices(number_of_vertices());
    for (VertexId vertex_id_1 = 0;
            vertex_id_1 < number_of_vertices();
            ++vertex_id_1) {
        for (VertexId vertex_id_2 = 0;
                vertex_id_2 < vertex_id_1;
                ++vertex_id_2) {
            distances_explicit_triangle_builder.set_distance(
                    vertex_id_1,
                    vertex_id_2,
                    distances.distance(vertex_id_1, vertex_id_2));
        }
    }
    distances_explicit_triangle_ = std::make_unique<const DistancesExplicitTriangle>(distances_explicit_triangle_builder.build());
}

template <typename T>
std::ostream& Distances::format(
        const T& distances,
        std::ostream& os,
        int verbosity_level) const
{
    if (verbosity_level >= 2) {
        os << std::endl
            << std::setw(12) << "Vertex 1"
            << std::setw(12) << "Vertex 2"
            << std::setw(12) << "Distance"
            << std::endl
            << std::setw(12) << "------"
            << std::setw(12) << "------"
            << std::setw(12) << "--------"
            << std::endl;
        for (VertexId vertex_id_1 = 0;
                vertex_id_1 < number_of_vertices();
                ++vertex_id_1) {
            for (VertexId vertex_id_2 = vertex_id_1 + 1;
                    vertex_id_2 < number_of_vertices();
                    ++vertex_id_2) {
                os
                    << std::setw(12) << vertex_id_1
                    << std::setw(12) << vertex_id_2
                    << std::setw(12) << distances.distance(vertex_id_1, vertex_id_2)
                    << std::endl;
            }
        }
    }

    return os;
}

}
