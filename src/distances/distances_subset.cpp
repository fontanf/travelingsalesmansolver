#include "travelingsalesmansolver/distances/distances_subset.hpp"

#include "travelingsalesmansolver/distances/distances_builder.hpp"

using namespace travelingsalesmansolver;

namespace
{

Distances distances_subset_impl(
        const DistancesExplicit& distances,
        const std::vector<VertexId>& vertex_ids)
{
    VertexId number_of_vertices = (VertexId)vertex_ids.size();
    DistancesExplicitBuilder distances_explicit_builder;
    distances_explicit_builder.set_number_of_vertices(number_of_vertices);
    for (VertexId vertex_id_1 = 0; vertex_id_1 < number_of_vertices; ++vertex_id_1) {
        for (VertexId vertex_id_2 = 0; vertex_id_2 < number_of_vertices; ++vertex_id_2) {
            if (vertex_id_1 == vertex_id_2)
                continue;
            distances_explicit_builder.set_distance(
                    vertex_id_1,
                    vertex_id_2,
                    distances.distance(vertex_ids[vertex_id_1], vertex_ids[vertex_id_2]));
        }
    }
    DistancesBuilder distances_builder;
    distances_builder.set_number_of_vertices(number_of_vertices);
    distances_builder.set_distances_explicit(distances_explicit_builder.build());
    return distances_builder.build();
}

Distances distances_subset_impl(
        const DistancesExplicitTriangle& distances,
        const std::vector<VertexId>& vertex_ids)
{
    VertexId number_of_vertices = (VertexId)vertex_ids.size();
    DistancesExplicitTriangleBuilder distances_explicit_triangle_builder;
    distances_explicit_triangle_builder.set_number_of_vertices(number_of_vertices);
    for (VertexId vertex_id_1 = 0; vertex_id_1 < number_of_vertices; ++vertex_id_1) {
        for (VertexId vertex_id_2 = 0; vertex_id_2 < vertex_id_1; ++vertex_id_2) {
            distances_explicit_triangle_builder.set_distance(
                    vertex_id_1,
                    vertex_id_2,
                    distances.distance(vertex_ids[vertex_id_1], vertex_ids[vertex_id_2]));
        }
    }
    DistancesBuilder distances_builder;
    distances_builder.set_number_of_vertices(number_of_vertices);
    distances_builder.set_distances_explicit_triangle(distances_explicit_triangle_builder.build());
    return distances_builder.build();
}

Distances distances_subset_impl(
        const DistancesEuc2D& distances,
        const std::vector<VertexId>& vertex_ids)
{
    VertexId number_of_vertices = (VertexId)vertex_ids.size();
    DistancesEuc2DBuilder distances_euc_2d_builder;
    distances_euc_2d_builder.set_number_of_vertices(number_of_vertices);
    distances_euc_2d_builder.set_scale(distances.scale());
    for (VertexId vertex_id = 0; vertex_id < number_of_vertices; ++vertex_id) {
        const Coordinates2D& coordinates = distances.coordinates(vertex_ids[vertex_id]);
        distances_euc_2d_builder.set_coordinates(vertex_id, coordinates.x, coordinates.y);
    }
    DistancesBuilder distances_builder;
    distances_builder.set_number_of_vertices(number_of_vertices);
    distances_builder.set_distances_euc_2d(distances_euc_2d_builder.build());
    return distances_builder.build();
}

Distances distances_subset_impl(
        const DistancesCeil2D& distances,
        const std::vector<VertexId>& vertex_ids)
{
    VertexId number_of_vertices = (VertexId)vertex_ids.size();
    DistancesCeil2DBuilder distances_ceil_2d_builder;
    distances_ceil_2d_builder.set_number_of_vertices(number_of_vertices);
    for (VertexId vertex_id = 0; vertex_id < number_of_vertices; ++vertex_id) {
        const Coordinates2D& coordinates = distances.coordinates(vertex_ids[vertex_id]);
        distances_ceil_2d_builder.set_coordinates(vertex_id, coordinates.x, coordinates.y);
    }
    DistancesBuilder distances_builder;
    distances_builder.set_number_of_vertices(number_of_vertices);
    distances_builder.set_distances_ceil_2d(distances_ceil_2d_builder.build());
    return distances_builder.build();
}

Distances distances_subset_impl(
        const DistancesFloor2D& distances,
        const std::vector<VertexId>& vertex_ids)
{
    VertexId number_of_vertices = (VertexId)vertex_ids.size();
    DistancesFloor2DBuilder distances_floor_2d_builder;
    distances_floor_2d_builder.set_number_of_vertices(number_of_vertices);
    distances_floor_2d_builder.set_scale(distances.scale());
    for (VertexId vertex_id = 0; vertex_id < number_of_vertices; ++vertex_id) {
        const Coordinates2D& coordinates = distances.coordinates(vertex_ids[vertex_id]);
        distances_floor_2d_builder.set_coordinates(vertex_id, coordinates.x, coordinates.y);
    }
    DistancesBuilder distances_builder;
    distances_builder.set_number_of_vertices(number_of_vertices);
    distances_builder.set_distances_floor_2d(distances_floor_2d_builder.build());
    return distances_builder.build();
}

Distances distances_subset_impl(
        const DistancesMan2D& distances,
        const std::vector<VertexId>& vertex_ids)
{
    VertexId number_of_vertices = (VertexId)vertex_ids.size();
    DistancesMan2DBuilder distances_man_2d_builder;
    distances_man_2d_builder.set_number_of_vertices(number_of_vertices);
    for (VertexId vertex_id = 0; vertex_id < number_of_vertices; ++vertex_id) {
        const Coordinates2D& coordinates = distances.coordinates(vertex_ids[vertex_id]);
        distances_man_2d_builder.set_coordinates(vertex_id, coordinates.x, coordinates.y);
    }
    DistancesBuilder distances_builder;
    distances_builder.set_number_of_vertices(number_of_vertices);
    distances_builder.set_distances_man_2d(distances_man_2d_builder.build());
    return distances_builder.build();
}

Distances distances_subset_impl(
        const DistancesGeo& distances,
        const std::vector<VertexId>& vertex_ids)
{
    VertexId number_of_vertices = (VertexId)vertex_ids.size();
    DistancesGeoBuilder distances_geo_builder;
    distances_geo_builder.set_number_of_vertices(number_of_vertices);
    for (VertexId vertex_id = 0; vertex_id < number_of_vertices; ++vertex_id) {
        const Coordinates2D& coordinates = distances.coordinates(vertex_ids[vertex_id]);
        distances_geo_builder.set_coordinates(vertex_id, coordinates.x, coordinates.y);
    }
    DistancesBuilder distances_builder;
    distances_builder.set_number_of_vertices(number_of_vertices);
    distances_builder.set_distances_geo(distances_geo_builder.build());
    return distances_builder.build();
}

Distances distances_subset_impl(
        const DistancesAtt& distances,
        const std::vector<VertexId>& vertex_ids)
{
    VertexId number_of_vertices = (VertexId)vertex_ids.size();
    DistancesAttBuilder distances_att_builder;
    distances_att_builder.set_number_of_vertices(number_of_vertices);
    for (VertexId vertex_id = 0; vertex_id < number_of_vertices; ++vertex_id) {
        const Coordinates2D& coordinates = distances.coordinates(vertex_ids[vertex_id]);
        distances_att_builder.set_coordinates(vertex_id, coordinates.x, coordinates.y);
    }
    DistancesBuilder distances_builder;
    distances_builder.set_number_of_vertices(number_of_vertices);
    distances_builder.set_distances_att(distances_att_builder.build());
    return distances_builder.build();
}

}

Distances travelingsalesmansolver::distances_subset(
        const Distances& distances,
        const std::vector<VertexId>& vertex_ids)
{
    return FUNCTION_WITH_DISTANCES(
            distances_subset_impl,
            distances,
            vertex_ids);
}
