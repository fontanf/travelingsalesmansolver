#include "travelingsalesmansolver/solution_builder.hpp"

#include <algorithm>

using namespace travelingsalesmansolver;

namespace
{

template <typename Distances>
Distance compute_distance(
        const Distances& distances,
        const std::vector<VertexId>& vertex_ids)
{
    if (vertex_ids.size() <= 1)
        return 0;
    Distance distance = 0;
    for (VertexPos pos = 0; pos < (VertexPos)vertex_ids.size() - 1; ++pos)
        distance += distances.distance(vertex_ids[pos], vertex_ids[pos + 1]);
    distance += distances.distance(vertex_ids.back(), vertex_ids.front());
    return distance;
}

}

SolutionBuilder::SolutionBuilder(const Instance& instance):
    instance_(&instance),
    vertices_is_added_(instance.number_of_vertices(), false)
{
}

void SolutionBuilder::add_vertex(VertexId vertex_id)
{
    if (vertex_id < 0 || vertex_id >= instance_->number_of_vertices()) {
        throw std::invalid_argument(
                "travelingsalesmansolver::SolutionBuilder::add_vertex: "
                "invalid vertex id " + std::to_string(vertex_id) + ".");
    }
    if (vertices_is_added_[vertex_id]) {
        throw std::invalid_argument(
                "travelingsalesmansolver::SolutionBuilder::add_vertex: "
                "vertex " + std::to_string(vertex_id)
                + " has already been added.");
    }
    vertices_is_added_[vertex_id] = true;
    vertex_ids_.push_back(vertex_id);
}

void SolutionBuilder::read(const std::string& certificate_path)
{
    std::ifstream file(certificate_path);
    if (!file.good()) {
        throw std::runtime_error(
                "Unable to open file \"" + certificate_path + "\".");
    }

    VertexId vertex_id = -1;
    while (file >> vertex_id) {
        try {
            add_vertex(vertex_id - 1);
        } catch (const std::invalid_argument&) {
            throw std::runtime_error(
                    "Invalid or duplicated vertex id "
                    + std::to_string(vertex_id)
                    + " in certificate \"" + certificate_path + "\".");
        }
    }
    if (!file.eof()) {
        throw std::runtime_error(
                "Unable to parse certificate \"" + certificate_path + "\".");
    }
}

Solution SolutionBuilder::build()
{
    // Start the tour at vertex 0.
    auto it = std::find(vertex_ids_.begin(), vertex_ids_.end(), 0);
    if (it != vertex_ids_.end()) {
        std::rotate(vertex_ids_.begin(), it, vertex_ids_.end());
    } else {
        vertex_ids_.insert(vertex_ids_.begin(), 0);
    }

    Solution solution(*instance_);
    solution.distance_ = FUNCTION_WITH_DISTANCES(
            compute_distance,
            instance_->distances(),
            vertex_ids_);
    solution.vertex_ids_ = std::move(vertex_ids_);
    return solution;
}
