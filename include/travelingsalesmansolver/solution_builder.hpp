#pragma once

#include "travelingsalesmansolver/solution.hpp"

namespace travelingsalesmansolver
{

/**
 * Builder class for a 'travelingsalesman' solution.
 *
 * Vertices are only stored when added; the length of the tour is computed
 * once, in 'build'.
 */
class SolutionBuilder
{

public:

    /** Constructor. */
    SolutionBuilder(const Instance& instance);

    /** Append a vertex at the end of the tour. */
    void add_vertex(VertexId vertex_id);

    /**
     * Append the vertices of a certificate file at the end of the tour.
     *
     * Vertex ids are 1-based, as written by 'Solution::write'.
     */
    void read(const std::string& certificate_path);

    /*
     * Build
     */

    /**
     * Build the solution.
     *
     * The tour is rotated so that it starts at vertex 0. If vertex 0 has not
     * been added, it is inserted at the start of the tour.
     */
    Solution build();

private:

    /*
     * Private attributes
     */

    /** Instance. */
    const Instance* instance_;

    /** List of added vertices. */
    std::vector<VertexId> vertex_ids_;

    /**
     * Array indexed by vertices indicating whether of not they have been
     * added.
     */
    std::vector<uint8_t> vertices_is_added_;

};

}
