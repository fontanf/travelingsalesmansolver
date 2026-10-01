#pragma once

#include "travelingsalesmansolver/algorithm_formatter.hpp"

#include <random>

namespace travelingsalesmansolver
{

struct RandomPermutationParameters: Parameters
{
};

/** Random tour: a uniformly random permutation of the vertices. */
const Output random_permutation(
        const Instance& instance,
        std::mt19937_64& generator,
        const RandomPermutationParameters& parameters = {});

}
