#include "travelingsalesmansolver/distances/distances.hpp"

using namespace travelingsalesmansolver;

void Distances::compute_distances_explicit() const
{
    FUNCTION_WITH_DISTANCES_0(this->compute_distances_explicit, *this);
}

void Distances::compute_distances_explicit_triangle() const
{
    FUNCTION_WITH_DISTANCES_0(this->compute_distances_explicit_triangle, *this);
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
    FUNCTION_WITH_DISTANCES_R(this->write, *this, file);
}
