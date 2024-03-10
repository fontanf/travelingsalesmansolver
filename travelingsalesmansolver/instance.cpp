#include "travelingsalesmansolver/instance.hpp"

#include "travelingsalesmansolver/distances/distances_builder.hpp"

#include "optimizationtools/utils/utils.hpp"
#include "optimizationtools/containers//indexed_set.hpp"

using namespace travelingsalesmansolver;

Instance::Instance(
        const std::string& instance_path,
        const std::string& format)
{
    std::ifstream file(instance_path);
    if (!file.good()) {
        throw std::runtime_error(
                "Unable to open file \"" + instance_path + "\".");
    }
    if (format == ""
            || format == "default"
            || format == "tsplib") {
        read_tsplib(file);
    } else {
        throw std::invalid_argument(
                "Unknown instance format \"" + format + "\".");
    }
}

void Instance::read_tsplib(std::ifstream& file)
{
    travelingsalesmansolver::DistancesBuilder distances_builder;

    std::string tmp;
    std::vector<std::string> line;
    while (getline(file, tmp)) {
        line = optimizationtools::split(tmp);
        if (line.size() == 0) {
        } else if (distances_builder.read_tsplib(file, tmp, line)) {
        } else if (tmp.rfind("NAME", 0) == 0) {
        } else if (tmp.rfind("PROBLEM NAME", 0) == 0) {
        } else if (tmp.rfind("COMMENT", 0) == 0) {
        } else if (tmp.rfind("TYPE", 0) == 0) {
        } else if (tmp.rfind("DIMENSION", 0) == 0) {
            VertexId number_of_vertices = std::stol(line.back());
            distances_builder.set_number_of_vertices(number_of_vertices);
        } else if (tmp.rfind("EOF", 0) == 0) {
            break;
        } else {
            throw std::invalid_argument(
                    "ENTRY \""
                    + line[0]
                    + "\" not implemented.");
        }
    }

    distances_ = std::shared_ptr<const travelingsalesmansolver::Distances>(
            new travelingsalesmansolver::Distances(distances_builder.build()));
}

std::ostream& Instance::format(
        std::ostream& os,
        int verbosity_level) const
{
    return FUNCTION_WITH_DISTANCES(
            (this->Instance::format),
            *distances_,
            os,
            verbosity_level);
}

void Instance::write(
        const std::string& instance_path) const
{
    if (instance_path.empty())
        return;
    std::ofstream file(instance_path);
    if (!file.good()) {
        throw std::runtime_error(
                "Unable to open file \"" + instance_path + "\".");
    }

    file << "NAME: XXX" << std::endl;
    file << "COMMENT: generated by fontanf/travelingsalesmansolver" << std::endl;
    file << "TYPE: TSP" << std::endl;
    file << "DIMENSION: " << number_of_vertices() << std::endl;
    distances().write(file);
    file << "EOF" << std::endl;
}


std::pair<bool, Distance> Instance::check(
        const std::string& certificate_path,
        std::ostream& os,
        int verbosity_level) const
{
    return FUNCTION_WITH_DISTANCES(
            (this->Instance::check),
            *distances_,
            certificate_path,
            os,
            verbosity_level);
}

template <typename Distances>
std::pair<bool, Distance> Instance::check(
        const Distances& distances,
        const std::string& certificate_path,
        std::ostream& os,
        int verbosity_level) const
{
    std::ifstream file(certificate_path);
    if (!file.good()) {
        throw std::runtime_error(
                "Unable to open file \"" + certificate_path + "\".");
    }

    if (verbosity_level >= 2) {
        os << std::endl
            << std::setw(12) << "Vertex"
            << std::setw(12) << "Distance"
            << std::endl
            << std::setw(12) << "------"
            << std::setw(12) << "--------"
            << std::endl;
    }

    VertexId vertex_id_pred = 0;
    VertexId vertex_id = -1;
    optimizationtools::IndexedSet vertices(number_of_vertices());
    vertices.add(0);
    VertexPos number_of_duplicates = 0;
    Distance total_distance = 0;
    while (file >> vertex_id) {

        // Check duplicates.
        if (vertices.contains(vertex_id)) {
            number_of_duplicates++;
            if (verbosity_level >= 2) {
                os << "Vertex " << vertex_id
                    << " has already been visited." << std::endl;
            }
        }
        vertices.add(vertex_id);

        total_distance += distances.distance(vertex_id_pred, vertex_id);

        if (verbosity_level >= 2) {
            os
                << std::setw(12) << vertex_id
                << std::setw(12) << total_distance
                << std::endl;
        }

        vertex_id_pred = vertex_id;
    }
    total_distance += distances.distance(vertex_id_pred, 0);

    bool feasible
        = (vertices.size() == number_of_vertices())
        && (number_of_duplicates == 0);

    if (verbosity_level >= 2)
        os << std::endl;
    if (verbosity_level >= 1) {
        os
            << "Number of vertices:     " << vertices.size() << " / " << number_of_vertices()  << std::endl
            << "Number of duplicates:   " << number_of_duplicates << std::endl
            << "Feasible:               " << feasible << std::endl
            << "Total distance:         " << total_distance << std::endl
            ;
    }
    return {feasible, total_distance};
}
