#include "travelingsalesmansolver/distances/distances_builder.hpp"

using namespace travelingsalesmansolver;

bool DistancesBuilder::read_tsplib(
        std::ifstream& file,
        const std::string& tmp,
        const std::vector<std::string>& line)
{
    std::unique_ptr<DistancesExplicitBuilder> distances_explicit_builder = nullptr;

    if (tmp.rfind("EDGE_WEIGHT_TYPE", 0) == 0) {
        edge_weight_type_ = line.back();

    } else if (tmp.rfind("EDGE_WEIGHT_FORMAT", 0) == 0) {
        edge_weight_format_ = line.back();

    } else if (tmp.rfind("NODE_COORD_TYPE", 0) == 0) {
        node_coord_type_ = line.back();

    } else if (tmp.rfind("EDGE_WEIGHT_SECTION", 0) == 0) {
        if (edge_weight_format_ == "UPPER_ROW") {
            DistancesExplicitTriangleBuilder distances_explicit_triangle_builder;
            distances_explicit_triangle_builder.set_number_of_vertices(distances_.number_of_vertices_);
            Distance distance;
            for (VertexId vertex_id_1 = 0;
                    vertex_id_1 < distances_.number_of_vertices_ - 1;
                    ++vertex_id_1) {
                for (VertexId vertex_id_2 = vertex_id_1 + 1;
                        vertex_id_2 < distances_.number_of_vertices_;
                        ++vertex_id_2) {
                    file >> distance;
                    distances_explicit_triangle_builder.set_distance(
                            vertex_id_1,
                            vertex_id_2,
                            distance);
                }
            }
            set_distances_explicit_triangle(distances_explicit_triangle_builder.build());
        } else if (edge_weight_format_ == "LOWER_ROW") {
            DistancesExplicitTriangleBuilder distances_explicit_triangle_builder;
            distances_explicit_triangle_builder.set_number_of_vertices(distances_.number_of_vertices_);
            Distance distance;
            for (VertexId vertex_id_1 = 1;
                    vertex_id_1 < distances_.number_of_vertices_;
                    ++vertex_id_1) {
                for (VertexId vertex_id_2 = 0;
                        vertex_id_2 < vertex_id_1;
                        ++vertex_id_2) {
                    file >> distance;
                    distances_explicit_triangle_builder.set_distance(
                            vertex_id_1,
                            vertex_id_2,
                            distance);
                }
            }
            set_distances_explicit_triangle(distances_explicit_triangle_builder.build());
        } else if (edge_weight_format_ == "UPPER_DIAG_ROW") {
            DistancesExplicitTriangleBuilder distances_explicit_triangle_builder;
            distances_explicit_triangle_builder.set_number_of_vertices(distances_.number_of_vertices_);
            Distance distance;
            for (VertexId vertex_id_1 = 0;
                    vertex_id_1 < distances_.number_of_vertices_;
                    ++vertex_id_1) {
                for (VertexId vertex_id_2 = vertex_id_1;
                        vertex_id_2 < distances_.number_of_vertices_;
                        ++vertex_id_2) {
                    file >> distance;
                    distances_explicit_triangle_builder.set_distance(
                            vertex_id_1,
                            vertex_id_2,
                            distance);
                }
            }
            set_distances_explicit_triangle(distances_explicit_triangle_builder.build());
        } else if (edge_weight_format_ == "LOWER_DIAG_ROW") {
            DistancesExplicitTriangleBuilder distances_explicit_triangle_builder;
            distances_explicit_triangle_builder.set_number_of_vertices(distances_.number_of_vertices_);
            Distance distance;
            for (VertexId vertex_id_1 = 0;
                    vertex_id_1 < distances_.number_of_vertices_;
                    ++vertex_id_1) {
                for (VertexId vertex_id_2 = 0;
                        vertex_id_2 <= vertex_id_1;
                        ++vertex_id_2) {
                    file >> distance;
                    distances_explicit_triangle_builder.set_distance(
                            vertex_id_1,
                            vertex_id_2,
                            distance);
                }
            }
            set_distances_explicit_triangle(distances_explicit_triangle_builder.build());
        } else if (edge_weight_format_ == "FULL_MATRIX") {
            DistancesExplicitBuilder distances_explicit_builder;
            distances_explicit_builder.set_number_of_vertices(distances_.number_of_vertices_);
            Distance distance;
            for (VertexId vertex_id_1 = 0;
                    vertex_id_1 < distances_.number_of_vertices_;
                    ++vertex_id_1) {
                for (VertexId vertex_id_2 = 0;
                        vertex_id_2 < distances_.number_of_vertices_;
                        ++vertex_id_2) {
                    file >> distance;
                    distances_explicit_builder.set_distance(
                            vertex_id_1,
                            vertex_id_2,
                            distance);
                }
            }
            set_distances_explicit(distances_explicit_builder.build());
        } else {
            throw std::invalid_argument(
                    "EDGE_WEIGHT_FORMAT \""
                    + edge_weight_format_
                    + "\" not implemented.");
        }

    } else if (tmp.rfind("NODE_COORD_SECTION", 0) == 0) {
        if (node_coord_type_ == "TWOD_COORDS") {
            std::vector<Coordinates2D> coordinates_2d;
            VertexId tmp;
            double x, y;
            for (VertexId vertex_id = 0;
                    vertex_id < distances_.number_of_vertices_;
                    ++vertex_id) {
                file >> tmp >> x >> y;
                coordinates_2d.push_back({x, y});
            }

            if (edge_weight_type_ == "EUC_2D") {
                DistancesEuc2DBuilder distances_euc_2d_builder;
                distances_euc_2d_builder.set_number_of_vertices(distances_.number_of_vertices_);
                for (VertexId vertex_id = 0;
                        vertex_id < distances_.number_of_vertices_;
                        ++vertex_id) {
                    distances_euc_2d_builder.set_coordinates(
                            vertex_id,
                            coordinates_2d[vertex_id].x,
                            coordinates_2d[vertex_id].y);
                }
                set_distances_euc_2d(distances_euc_2d_builder.build());

            } else if (edge_weight_type_ == "CEIL_2D") {
                DistancesCeil2DBuilder distances_ceil_2d_builder;
                distances_ceil_2d_builder.set_number_of_vertices(distances_.number_of_vertices_);
                for (VertexId vertex_id = 0;
                        vertex_id < distances_.number_of_vertices_;
                        ++vertex_id) {
                    distances_ceil_2d_builder.set_coordinates(
                            vertex_id,
                            coordinates_2d[vertex_id].x,
                            coordinates_2d[vertex_id].y);
                }
                set_distances_ceil_2d(distances_ceil_2d_builder.build());

            } else if (edge_weight_type_ == "GEO") {
                DistancesGeoBuilder distances_geo_builder;
                distances_geo_builder.set_number_of_vertices(distances_.number_of_vertices_);
                for (VertexId vertex_id = 0;
                        vertex_id < distances_.number_of_vertices_;
                        ++vertex_id) {
                    distances_geo_builder.set_coordinates(
                            vertex_id,
                            coordinates_2d[vertex_id].x,
                            coordinates_2d[vertex_id].y);
                }
                set_distances_geo(distances_geo_builder.build());

            } else if (edge_weight_type_ == "ATT") {
                DistancesAttBuilder distances_att_builder;
                distances_att_builder.set_number_of_vertices(distances_.number_of_vertices_);
                for (VertexId vertex_id = 0;
                        vertex_id < distances_.number_of_vertices_;
                        ++vertex_id) {
                    distances_att_builder.set_coordinates(
                            vertex_id,
                            coordinates_2d[vertex_id].x,
                            coordinates_2d[vertex_id].y);
                }
                set_distances_att(distances_att_builder.build());

            }
        } else if (node_coord_type_ == "THREED_COORDS") {
            std::vector<Coordinates3D> coordinates_3d;
            VertexId tmp;
            double x, y, z;
            for (VertexId vertex_id = 0;
                    vertex_id < distances_.number_of_vertices_;
                    ++vertex_id) {
                file >> tmp >> x >> y >> z;
                coordinates_3d.push_back({x, y, z});
            }
        }

    } else if (tmp.rfind("DISPLAY_DATA_TYPE", 0) == 0) {

    } else if (tmp.rfind("DISPLAY_DATA_SECTION", 0) == 0) {
        VertexId tmp;
        double x, y;
        for (VertexId vertex_id = 0;
                vertex_id < distances_.number_of_vertices_;
                ++vertex_id) {
            file >> tmp >> x >> y;
            // TODO
        }

    } else if (tmp.rfind("EOF", 0) == 0) {

    } else {
        return false;

    }

    return true;
}

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////// Build /////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

Distances DistancesBuilder::build()
{
    if (!set_) {
        throw std::invalid_argument("No distances has already been set.");
    }
    set_ = false;
    return std::move(distances_);
}
