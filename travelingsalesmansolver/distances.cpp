#include "travelingsalesmansolver/distances.hpp"

#include "optimizationtools/utils/utils.hpp"

using namespace travelingsalesmansolver;

void Distances::compute_distances() const
{
    if (!distances_.empty())
        return;
    init_distances();

    if (edge_weight_type_ == "EUC_2D") {
        for (VertexId vertex_id_1 = 0;
                vertex_id_1 < number_of_vertices();
                ++vertex_id_1) {
            for (VertexId vertex_id_2 = vertex_id_1 + 1;
                    vertex_id_2 < number_of_vertices();
                    ++vertex_id_2) {
                double xd = x(vertex_id_2) - x(vertex_id_1);
                double yd = y(vertex_id_2) - y(vertex_id_1);
                Distance distance = std::round(std::sqrt(xd * xd + yd * yd));
                set_distance(vertex_id_1, vertex_id_2, distance);
            }
        }
    } else if (edge_weight_type_ == "CEIL_2D") {
        for (VertexId vertex_id_1 = 0;
                vertex_id_1 < number_of_vertices();
                ++vertex_id_1) {
            for (VertexId vertex_id_2 = vertex_id_1 + 1;
                    vertex_id_2 < number_of_vertices();
                    ++vertex_id_2) {
                double xd = x(vertex_id_2) - x(vertex_id_1);
                double yd = y(vertex_id_2) - y(vertex_id_1);
                Distance distance = std::ceil(std::sqrt(xd * xd + yd * yd));
                set_distance(vertex_id_1, vertex_id_2, distance);
            }
        }
    } else if (edge_weight_type_ == "GEO") {
        latitudes_ = std::vector<double>(number_of_vertices(), 0);
        longitudes_ = std::vector<double>(number_of_vertices(), 0);
        for (VertexId vertex_id = 0;
                vertex_id < number_of_vertices();
                ++vertex_id) {
            double pi = 3.141592;
            int deg_x = std::round(x(vertex_id));
            double min_x = x(vertex_id) - deg_x;
            latitudes_[vertex_id] = pi * (deg_x + 5.0 * min_x / 3.0) / 180.0;
            int deg_y = std::round(y(vertex_id));
            double min_y = y(vertex_id) - deg_y;
            longitudes_[vertex_id] = pi * (deg_y + 5.0 * min_y / 3.0) / 180.0;
        }
        double rrr = 6378.388;
        for (VertexId vertex_id_1 = 0;
                vertex_id_1 < number_of_vertices();
                ++vertex_id_1) {
            for (VertexId vertex_id_2 = vertex_id_1 + 1;
                    vertex_id_2 < number_of_vertices();
                    ++vertex_id_2) {
                double q1 = cos(longitudes_[vertex_id_1] - longitudes_[vertex_id_2]);
                double q2 = cos(latitudes_[vertex_id_1] - latitudes_[vertex_id_2]);
                double q3 = cos(latitudes_[vertex_id_1] + latitudes_[vertex_id_2]);
                Distance distance = (Distance)(rrr * acos(0.5 * ((1.0 + q1) * q2 - (1.0 - q1) * q3)) + 1.0);
                set_distance(vertex_id_1, vertex_id_2, distance);
            }
        }
    } else if (edge_weight_type_ == "ATT") {
        for (VertexId vertex_id_1 = 0;
                vertex_id_1 < number_of_vertices();
                ++vertex_id_1) {
            for (VertexId vertex_id_2 = vertex_id_1 + 1;
                    vertex_id_2 < number_of_vertices();
                    ++vertex_id_2) {
                double xd = x(vertex_id_1) - x(vertex_id_2);
                double yd = y(vertex_id_1) - y(vertex_id_2);
                double rij = sqrt((xd * xd + yd * yd) / 10.0);
                int tij = std::round(rij);
                Distance distance = (tij < rij)? tij + 1: tij;
                set_distance(vertex_id_1, vertex_id_2, distance);
            }
        }
    } else {
        throw std::invalid_argument(
                "EDGE_WEIGHT_TYPE \""
                + edge_weight_type_
                + "\" not implemented.");
    }
    for (VertexId vertex_id = 0;
            vertex_id < number_of_vertices();
            ++vertex_id) {
        distances_[vertex_id][vertex_id] = std::numeric_limits<Distance>::max();
    }
}

std::ostream& Distances::print(
        std::ostream& os,
        int verbose) const
{
    if (verbose >= 2) {
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
                    << std::setw(12) << distance(vertex_id_1, vertex_id_2)
                    << std::endl;
            }
        }
    }

    return os;
}

void Distances::write(std::ofstream& file) const
{
    file << "DIMENSION: " << number_of_vertices() << std::endl;

    file << "EDGE_WEIGHT_TYPE: " << edge_weight_type_ << std::endl;
    if (edge_weight_type_ == "EXPLICIT") {
        file << "EDGE_WEIGHT_FORMAT: UPPER_ROW" << std::endl;
        file << "EDGE_WEIGHT_SECTION" << std::endl;
        for (VertexId vertex_id_1 = 0;
                vertex_id_1 < number_of_vertices() - 1;
                ++vertex_id_1) {
            for (VertexId vertex_id_2 = vertex_id_1 + 1;
                    vertex_id_2 < number_of_vertices();
                    ++vertex_id_2) {
                file << distance(vertex_id_1, vertex_id_2) << " ";
            }
        }
        file << std::endl;
    }

    file << "NODE_COORD_TYPE " << node_coord_type_ << std::endl;
    if (node_coord_type_ == "TWOD_COORDS") {
        file << "NODE_COORD_SECTION" << std::endl;
        for (VertexId vertex_id = 0;
                vertex_id < number_of_vertices();
                ++vertex_id) {
            const Vertex& vertex = this->vertex(vertex_id);
            file << vertex_id + 1
                << " " << vertex.x << " "
                << vertex.y << std::endl;
        }
    } else if (node_coord_type_ == "THREED_COORDS") {
        file << "NODE_COORD_SECTION" << std::endl;
        for (VertexId vertex_id = 0;
                vertex_id < number_of_vertices();
                ++vertex_id) {
            const Vertex& vertex = this->vertex(vertex_id);
            file << vertex_id + 1 << " "
                << vertex.x << " "
                << vertex.y << " "
                << vertex.z << std::endl;
        }
    }
}
