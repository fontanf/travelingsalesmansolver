#include "travelingsalesmansolver/distances_builder.hpp"

#include "optimizationtools/utils/utils.hpp"

using namespace travelingsalesmansolver;

void DistancesBuilder::add_vertices(VertexId number_of_vertices)
{
    distances_.vertices_.insert(
            distances_.vertices_.end(),
            number_of_vertices,
            Vertex());
}

bool DistancesBuilder::read_tsplib(
        std::ifstream& file,
        const std::string& tmp,
        const std::vector<std::string>& line)
{
    if (tmp.rfind("EDGE_WEIGHT_TYPE", 0) == 0) {
        distances_.edge_weight_type_ = line.back();
    } else if (tmp.rfind("EDGE_WEIGHT_FORMAT", 0) == 0) {
        distances_.edge_weight_format_ = line.back();
    } else if (tmp.rfind("NODE_COORD_TYPE", 0) == 0) {
        distances_.node_coord_type_ = line.back();
    } else if (tmp.rfind("EDGE_WEIGHT_SECTION", 0) == 0) {
        distances_.init_distances();
        if (distances_.edge_weight_format_ == "UPPER_ROW") {
            Distance distance;
            for (VertexId vertex_id_1 = 0;
                    vertex_id_1 < distances_.number_of_vertices() - 1;
                    ++vertex_id_1) {
                for (VertexId vertex_id_2 = vertex_id_1 + 1;
                        vertex_id_2 < distances_.number_of_vertices();
                        ++vertex_id_2) {
                    file >> distance;
                    set_distance(vertex_id_1, vertex_id_2, distance);
                }
            }
        } else if (distances_.edge_weight_format_ == "LOWER_ROW") {
            Distance distance;
            for (VertexId vertex_id_1 = 1;
                    vertex_id_1 < distances_.number_of_vertices();
                    ++vertex_id_1) {
                for (VertexId vertex_id_2 = 0;
                        vertex_id_2 < vertex_id_1;
                        ++vertex_id_2) {
                    file >> distance;
                    set_distance(vertex_id_1, vertex_id_2, distance);
                }
            }
        } else if (distances_.edge_weight_format_ == "UPPER_DIAG_ROW") {
            Distance distance;
            for (VertexId vertex_id_1 = 0;
                    vertex_id_1 < distances_.number_of_vertices();
                    ++vertex_id_1) {
                for (VertexId vertex_id_2 = vertex_id_1;
                        vertex_id_2 < distances_.number_of_vertices();
                        ++vertex_id_2) {
                    file >> distance;
                    set_distance(vertex_id_1, vertex_id_2, distance);
                }
            }
        } else if (distances_.edge_weight_format_ == "LOWER_DIAG_ROW") {
            Distance distance;
            for (VertexId vertex_id_1 = 0;
                    vertex_id_1 < distances_.number_of_vertices();
                    ++vertex_id_1) {
                for (VertexId vertex_id_2 = 0;
                        vertex_id_2 <= vertex_id_1;
                        ++vertex_id_2) {
                    file >> distance;
                    set_distance(vertex_id_1, vertex_id_2, distance);
                }
            }
        } else if (distances_.edge_weight_format_ == "FULL_MATRIX") {
            Distance distance;
            for (VertexId vertex_id_1 = 0;
                    vertex_id_1 < distances_.number_of_vertices();
                    ++vertex_id_1) {
                for (VertexId vertex_id_2 = 0;
                        vertex_id_2 < distances_.number_of_vertices();
                        ++vertex_id_2) {
                    file >> distance;
                    set_distance(vertex_id_1, vertex_id_2, distance);
                }
            }
        } else {
            throw std::invalid_argument(
                    "EDGE_WEIGHT_FORMAT \""
                    + distances_.edge_weight_format_
                    + "\" not implemented.");
        }
    } else if (tmp.rfind("NODE_COORD_SECTION", 0) == 0) {
        if (distances_.node_coord_type_ == "TWOD_COORDS") {
            VertexId tmp;
            double x, y;
            for (VertexId vertex_id = 0;
                    vertex_id < distances_.number_of_vertices();
                    ++vertex_id) {
                file >> tmp >> x >> y;
                set_coordinates(vertex_id, x, y);
            }
        } else if (distances_.node_coord_type_ == "THREED_COORDS") {
            VertexId tmp;
            double x, y, z;
            for (VertexId vertex_id = 0;
                    vertex_id < distances_.number_of_vertices();
                    ++vertex_id) {
                file >> tmp >> x >> y >> z;
                set_coordinates(vertex_id, x, y, z);
            }
        }
    } else if (tmp.rfind("DISPLAY_DATA_TYPE", 0) == 0) {
    } else if (tmp.rfind("DISPLAY_DATA_SECTION", 0) == 0) {
        VertexId tmp;
        double x, y;
        for (VertexId vertex_id = 0;
                vertex_id < distances_.number_of_vertices();
                ++vertex_id) {
            file >> tmp >> x >> y;
            set_coordinates(vertex_id, x, y);
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
    return std::move(distances_);
}
