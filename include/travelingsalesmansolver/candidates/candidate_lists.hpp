#pragma once

#include "travelingsalesmansolver/distances/commons.hpp"

#include <vector>

namespace travelingsalesmansolver
{

/**
 * Candidate lists: for each vertex, the vertices its edges are tried with,
 * best first.
 */
using CandidateLists = std::vector<std::vector<VertexId>>;

}
