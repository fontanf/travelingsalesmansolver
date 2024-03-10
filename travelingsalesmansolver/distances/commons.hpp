#pragma once

#include <cstdint>

namespace travelingsalesmansolver
{

using VertexId = int64_t;
using Distance = int64_t;

struct Coordinates2D
{
    /** x-coordinate of the location. */
    double x;

    /** y-coordinate of the location. */
    double y;
};

struct Coordinates3D
{
    /** x-coordinate of the location. */
    double x;

    /** y-coordinate of the location. */
    double y;

    /** z-coordinate of the location. */
    double z;
};

}
