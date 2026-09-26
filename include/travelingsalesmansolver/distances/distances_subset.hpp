#pragma once

#include "travelingsalesmansolver/distances/distances.hpp"

namespace travelingsalesmansolver
{

/**
 * Build distances over a subset of the vertices of 'distances': local
 * vertex 'i' of the result corresponds to 'vertex_ids[i]' of the original
 * distances.
 *
 * The concrete representation is adapted to the input's: coordinate-based
 * types (Euc2D, Ceil2D, Man2D, Geo, Att) produce a distances
 * object of that same type, holding only the 'vertex_ids.size()' selected
 * vertices' coordinates (O(k), no distance ever computed or stored);
 * explicit distances (Explicit, ExplicitTriangle) have no coordinates to
 * carry over, so the subset's distances are computed and copied instead
 * (O(k^2), unavoidable for that representation).
 */
Distances distances_subset(
        const Distances& distances,
        const std::vector<VertexId>& vertex_ids);

}
