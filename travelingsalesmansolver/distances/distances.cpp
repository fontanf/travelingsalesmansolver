#include "travelingsalesmansolver/distances/distances.hpp"

#include <iomanip>

using namespace travelingsalesmansolver;

void Distances::compute_distances_explicit() const
{
    //if (!distances_.empty())
    //    return;
    //init_distances(triangle);

    //    for (VertexId vertex_id_1 = 0;
    //            vertex_id_1 < number_of_vertices();
    //            ++vertex_id_1) {
    //        for (VertexId vertex_id_2 = 0;
    //                vertex_id_2 < vertex_id_1;
    //                ++vertex_id_2) {
    //            distances_[vertex_id_1][vertex_id_2] = compute_distance(vertex_id_1, vertex_id_2);
    //        }
    //    }
    //for (VertexId vertex_id = 0;
    //        vertex_id < number_of_vertices();
    //        ++vertex_id) {
    //    distances_[vertex_id][vertex_id] = std::numeric_limits<Distance>::max();
    //}
}

std::ostream& Distances::format(
        std::ostream& os,
        int verbosity_level) const
{
    return FUNCTION_WITH_DISTANCES(this->format, *this, os, verbosity_level);
}

void Distances::write(std::ofstream& file) const
{
    file << "DIMENSION: " << number_of_vertices() << std::endl;
    FUNCTION_WITH_DISTANCES(this->write, *this, file);
}
