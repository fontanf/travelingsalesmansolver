#include "travelingsalesmansolver/algorithms/lkh.hpp"

using namespace travelingsalesmansolver;

const LkhOutput travelingsalesmansolver::lkh(
        const Instance& instance,
        const LkhParameters& parameters)
{
    return FUNCTION_WITH_DISTANCES(
            lkh,
            instance.distances(),
            instance,
            parameters);
}

std::vector<LkhCandidate> travelingsalesmansolver::read_candidates(
        const std::string& candidate_file_content)
{
    std::stringstream ss(candidate_file_content);
    VertexId number_of_vertices = -1;
    ss >> number_of_vertices;

    std::vector<LkhCandidate> candidates(number_of_vertices);

    VertexId vertex_id_tmp = -1;
    VertexId parent_id = -1;
    VertexId vertex_number_of_candidates = -1;
    VertexId vertex_id_2 = -1;
    double alpha = 0;
    for (VertexId vertex_id = 0; vertex_id < number_of_vertices; ++vertex_id) {
        ss >> vertex_id_tmp >> parent_id >> vertex_number_of_candidates;
        parent_id--;
        candidates[vertex_id].parent_id = parent_id;
        for (VertexId pos = 0; pos < vertex_number_of_candidates; ++pos) {
            ss >> vertex_id_2 >> alpha;
            vertex_id_2--;

            LkhCandidateEdge edge;
            edge.vertex_id = vertex_id_2;
            edge.alpha = alpha;
            candidates[vertex_id].edges.push_back(edge);
        }

    }

    return candidates;
}
