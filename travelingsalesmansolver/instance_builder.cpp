#include "travelingsalesmansolver/instance_builder.hpp"

using namespace travelingsalesmansolver;

void InstanceBuilder::add_vertices(VertexId number_of_vertices)
{
    instance_.vertices_.insert(
            instance_.vertices_.end(),
            number_of_vertices,
            Vertex());
}

InstanceBuilder::InstanceBuilder(
        std::string instance_path,
        std::string format)
{
    std::ifstream file(instance_path);
    if (!file.good()) {
        throw std::runtime_error(
                "Unable to open file \"" + instance_path + "\".");
    }

    if (format == "" || format == "tsplib") {
        read_tsplib(file);
    } else {
        throw std::invalid_argument(
                "Unknown instance format \"" + format + "\".");
    }
}

void InstanceBuilder::read_tsplib(std::ifstream& file)
{
    std::string tmp;
    std::vector<std::string> line;
    std::string edge_weight_format;
    while (getline(file, tmp)) {
        line = optimizationtools::split(tmp, ' ');
        if (line.size() == 0) {
        } else if (tmp.rfind("NAME", 0) == 0) {
        } else if (tmp.rfind("COMMENT", 0) == 0) {
        } else if (tmp.rfind("TYPE", 0) == 0) {
        } else if (tmp.rfind("DISPLAY_DATA_TYPE", 0) == 0) {
        } else if (tmp.rfind("DIMENSION", 0) == 0) {
            VertexId number_of_vertices = std::stol(line.back());
            add_vertices(number_of_vertices);
        } else if (tmp.rfind("EDGE_WEIGHT_TYPE", 0) == 0) {
            instance_.edge_weight_type_ = line.back();
        } else if (tmp.rfind("EDGE_WEIGHT_FORMAT", 0) == 0) {
            edge_weight_format = line.back();
        } else if (tmp.rfind("NODE_COORD_TYPE", 0) == 0) {
            instance_.node_coord_type_ = line.back();
        } else if (tmp.rfind("EDGE_WEIGHT_SECTION", 0) == 0) {
            instance_.init_distances();
            if (edge_weight_format == "UPPER_ROW") {
                Distance distance;
                for (VertexId vertex_id_1 = 0;
                        vertex_id_1 < instance_.number_of_vertices() - 1;
                        ++vertex_id_1) {
                    for (VertexId vertex_id_2 = vertex_id_1 + 1;
                            vertex_id_2 < instance_.number_of_vertices();
                            ++vertex_id_2) {
                        file >> distance;
                        set_distance(vertex_id_1, vertex_id_2, distance);
                    }
                }
            } else if (edge_weight_format == "LOWER_ROW") {
                Distance distance;
                for (VertexId vertex_id_1 = 1;
                        vertex_id_1 < instance_.number_of_vertices();
                        ++vertex_id_1) {
                    for (VertexId vertex_id_2 = 0;
                            vertex_id_2 < vertex_id_1;
                            ++vertex_id_2) {
                        file >> distance;
                        set_distance(vertex_id_1, vertex_id_2, distance);
                    }
                }
            } else if (edge_weight_format == "UPPER_DIAG_ROW") {
                Distance distance;
                for (VertexId vertex_id_1 = 0;
                        vertex_id_1 < instance_.number_of_vertices();
                        ++vertex_id_1) {
                    for (VertexId vertex_id_2 = vertex_id_1;
                            vertex_id_2 < instance_.number_of_vertices();
                            ++vertex_id_2) {
                        file >> distance;
                        set_distance(vertex_id_1, vertex_id_2, distance);
                    }
                }
            } else if (edge_weight_format == "LOWER_DIAG_ROW") {
                Distance distance;
                for (VertexId vertex_id_1 = 0;
                        vertex_id_1 < instance_.number_of_vertices();
                        ++vertex_id_1) {
                    for (VertexId vertex_id_2 = 0;
                            vertex_id_2 <= vertex_id_1;
                            ++vertex_id_2) {
                        file >> distance;
                        set_distance(vertex_id_1, vertex_id_2, distance);
                    }
                }
            } else if (edge_weight_format == "FULL_MATRIX") {
                Distance distance;
                for (VertexId vertex_id_1 = 0;
                        vertex_id_1 < instance_.number_of_vertices();
                        ++vertex_id_1) {
                    for (VertexId vertex_id_2 = 0;
                            vertex_id_2 < instance_.number_of_vertices();
                            ++vertex_id_2) {
                        file >> distance;
                        set_distance(vertex_id_1, vertex_id_2, distance);
                    }
                }
            } else {
                throw std::invalid_argument(
                        "EDGE_WEIGHT_FORMAT \""
                        + edge_weight_format
                        + "\" not implemented.");
            }
        } else if (tmp.rfind("NODE_COORD_SECTION", 0) == 0) {
            if (instance_.node_coord_type_ == "TWOD_COORDS") {
                VertexId tmp;
                double x, y;
                for (VertexId vertex_id = 0;
                        vertex_id < instance_.number_of_vertices();
                        ++vertex_id) {
                    file >> tmp >> x >> y;
                    set_coordinates(vertex_id, x, y);
                }
            } else if (instance_.node_coord_type_ == "THREED_COORDS") {
                VertexId tmp;
                double x, y, z;
                for (VertexId vertex_id = 0;
                        vertex_id < instance_.number_of_vertices();
                        ++vertex_id) {
                    file >> tmp >> x >> y >> z;
                    set_coordinates(vertex_id, x, y, z);
                }
            }
        } else if (tmp.rfind("DISPLAY_DATA_SECTION", 0) == 0) {
            VertexId tmp;
            double x, y;
            for (VertexId vertex_id = 0;
                    vertex_id < instance_.number_of_vertices();
                    ++vertex_id) {
                file >> tmp >> x >> y;
                set_coordinates(vertex_id, x, y);
            }
        } else if (tmp.rfind("EOF", 0) == 0) {
            break;
        } else {
            throw std::invalid_argument(
                    "ENTRY \""
                    + line[0]
                    + "\" not implemented.");
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////// Build /////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

Instance InstanceBuilder::build()
{
    return std::move(instance_);
}
