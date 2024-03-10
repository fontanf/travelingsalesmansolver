#include "travelingsalesmansolver/solution.hpp"

using namespace travelingsalesmansolver;

Solution::Solution(const Instance& instance):
    instance_(&instance),
    vertices_is_visited_(instance.number_of_vertices(), false)
{
    // Add initial vertex.
    vertex_ids_.push_back(0);
    vertices_is_visited_[0] = true;
}

Solution::Solution(
        const Instance& instance,
        const std::string& certificate_path):
    Solution(instance)
{
    if (certificate_path.empty())
        return;
    std::ifstream file(certificate_path);
    if (!file.good()) {
        throw std::runtime_error(
                "Unable to open file \"" + certificate_path + "\".");
    }

    // TODO
}

bool Solution::feasible() const
{
    return (number_of_vertices() == instance().number_of_vertices());
}

void Solution::format(
        std::ostream& os,
        int verbosity_level) const
{
    if (verbosity_level >= 1) {
        os
            << "Number of vertices:  " << optimizationtools::Ratio<VertexId>(number_of_vertices(), instance().number_of_vertices()) << std::endl
            << "Feasible:            " << feasible() << std::endl
            << "Distance:            " << distance() << std::endl
            ;
    }

    if (verbosity_level >= 2) {
        os << std::endl
            << std::setw(12) << "Vertex"
            << std::endl
            << std::setw(12) << "------"
            << std::endl;
        for (VertexId vertex_id: vertex_ids_) {
            os
                << std::setw(12) << vertex_id
                << std::endl;
        }
    }
}

nlohmann::json Solution::to_json() const
{
    return nlohmann::json {
        {"NumberOfVertices", number_of_vertices()},
        {"Feasible", feasible()},
        {"Distance", distance()}
    };
}

void Solution::write(
        const std::string& certificate_path) const
{
    if (certificate_path.empty())
        return;
    std::ofstream file(certificate_path);
    if (!file.good()) {
        throw std::runtime_error(
                "Unable to open file \"" + certificate_path + "\".");
    }

    for (VertexId vertex_id: vertex_ids_)
        file << " " << vertex_id + 1;
}
