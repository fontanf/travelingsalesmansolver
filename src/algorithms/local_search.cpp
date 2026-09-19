// The EAX genetic algorithm implemented in 'local_search.hpp' is adapted from
// Shujia Liu's C++ implementation (https://github.com/Sugia/GA-for-TSP),
// licensed under the Apache License, Version 2.0; see
// 'licenses/eax-ga/LICENSE' and 'licenses/eax-ga/NOTICE.md' for the license
// text and the full list of changes made to the original source.

#include "travelingsalesmansolver/algorithms/local_search.hpp"

#include <cmath>
#include <cstdlib>

using namespace travelingsalesmansolver;

namespace travelingsalesmansolver
{
namespace
{

bool Individual::operator==(
        const Individual& individual) const
{
    if (neighbors.size() != individual.neighbors.size())
        return false;
    if (length != individual.length)
        return false;

    VertexId current_vertex_id = 0;
    VertexId previous_vertex_id = -1;
    for (std::size_t i = 0; i < neighbors.size(); ++i) {
        VertexId next_vertex_id = (neighbors[current_vertex_id][0] == previous_vertex_id)? neighbors[current_vertex_id][1]: neighbors[current_vertex_id][0];
        if (individual.neighbors[current_vertex_id][0] != next_vertex_id && individual.neighbors[current_vertex_id][1] != next_vertex_id)
            return false;
        previous_vertex_id = current_vertex_id;
        current_vertex_id = next_vertex_id;
    }
    return true;
}

void seed_random(int seed)
{
    std::srand(seed);
}

int random_integer(int min, int max)
{
    return min + (std::rand() % (max - min + 1));
}

namespace
{

double random_double(double min, double max)
{
    return min + std::rand() % (int)(max - min);
}

}

double random_normal(double mu, double sigma)
{
    const double pi = 3.1415926;
    double u1;
    do {
        u1 = random_double(0.0, 1.0);
    } while (u1 == 0.0);
    double u2 = random_double(0.0, 1.0);
    double x = std::sqrt(-2.0 * std::log(u1)) * std::cos(2 * pi * u2);
    return mu + sigma * x;
}

void random_permutation(
        std::vector<int>& array,
        int number_of_elements,
        int number_of_samples)
{
    if (number_of_elements <= 0)
        return;
    std::vector<int> visited(number_of_elements, 0);
    for (int i = 0; i < number_of_samples; ++i) {
        int r = std::rand() % (number_of_elements - i);
        while (visited[r] == 1)
            r = (r + 1) % number_of_elements;
        array[i] = r;
        visited[r] = 1;
    }
}

void random_shuffle(
        std::vector<int>& array,
        int number_of_elements)
{
    std::vector<int> shuffled_positions(number_of_elements);
    random_permutation(shuffled_positions, number_of_elements, number_of_elements);
    std::vector<int> original(array.begin(), array.begin() + number_of_elements);
    for (int i = 0; i < number_of_elements; ++i)
        array[i] = original[shuffled_positions[i]];
}

namespace
{

void selection_sort(
        std::vector<int>& values,
        int l,
        int r)
{
    for (int i = l; i < r; ++i) {
        int id = i;
        for (int j = i + 1; j <= r; ++j)
            if (values[j] < values[id])
                id = j;
        std::swap(values[i], values[id]);
    }
}

int quick_sort_partition(
        std::vector<int>& values,
        int l,
        int r)
{
    int id = l + std::rand() % (r - l + 1);
    std::swap(values[l], values[id]);
    id = l;
    for (int i = l + 1; i <= r; ++i)
        if (values[i] < values[l])
            std::swap(values[++id], values[i]);
    std::swap(values[l], values[id]);
    return id;
}

void quick_sort(
        std::vector<int>& values,
        int l,
        int r)
{
    if (l < r) {
        if (r - l < 20) { // utilizes selection sort for small batch of data
            selection_sort(values, l, r);
            return;
        }
        int mid = quick_sort_partition(values, l, r);
        quick_sort(values, l, mid - 1);
        quick_sort(values, mid + 1, r);
    }
}

}

void sort_indices_ascending(
        const std::vector<int>& values,
        int number_of_values,
        std::vector<int>& sorted_indices,
        int number_of_indices)
{
    std::vector<int> checked(number_of_values, 0);
    for (int i = 0; i < number_of_indices; ++i) {
        int best_value = std::numeric_limits<int>::max();
        int best_index = 0;
        for (int j = 0; j < number_of_values; ++j) {
            if (values[j] < best_value && checked[j] == 0) {
                best_value = values[j];
                best_index = j;
            }
        }
        sorted_indices[i] = best_index;
        checked[best_index] = 1;
    }
}

void sort_indices_descending(
        const std::vector<int>& values,
        int number_of_values,
        std::vector<int>& sorted_indices,
        int number_of_indices)
{
    std::vector<int> checked(number_of_values, 0);
    for (int i = 0; i < number_of_indices; ++i) {
        int best_value = std::numeric_limits<int>::min();
        int best_index = 0;
        for (int j = 0; j < number_of_values; ++j) {
            if (values[j] > best_value && checked[j] == 0) {
                best_value = values[j];
                best_index = j;
            }
        }
        sorted_indices[i] = best_index;
        checked[best_index] = 1;
    }
}

void sort_ascending(
        std::vector<int>& values,
        int number_of_values)
{
    quick_sort(values, 0, number_of_values - 1);
}

}
}

const Output travelingsalesmansolver::local_search(
        const Instance& instance,
        const LocalSearchParameters& parameters)
{
    return FUNCTION_WITH_DISTANCES(
            local_search,
            instance.distances(),
            instance,
            parameters);
}
