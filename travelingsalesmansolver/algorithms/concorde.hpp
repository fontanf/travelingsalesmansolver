#pragma once

#include "travelingsalesmansolver/solution.hpp"

#include "travelingsalesmansolver/algorithm_formatter.hpp"

namespace travelingsalesmansolver
{

const Output concorde(
        const Instance& instance,
        const Parameters& parameters);

template <typename Distances>
const Output concorde(
        const Distances& distances,
        const Instance& instance,
        const Parameters& parameters = {});


template <typename Distances>
const Output concorde(
        const Distances& distances,
        const Instance& instance,
        const Parameters& parameters)
{
    Output output(instance);
    AlgorithmFormatter algorithm_formatter(parameters, output);
    algorithm_formatter.start("Concorde");
    algorithm_formatter.print_header();

    // Write instance file.
    char instance_path[L_tmpnam];
    tmpnam(instance_path);
    instance.write(instance_path);

    // Run.
    char output_path[L_tmpnam];
    tmpnam(output_path);
    char solution_path[L_tmpnam];
    tmpnam(solution_path);
    std::string command = (
            "concorde"
            " -x"  // delete files on completion (sav pul mas)
            " -o \"" + std::string(solution_path) + "\""
            + " \"" + std::string(instance_path) + "\""
            + " > \"" + std::string(output_path) + "\"");
    std::system(command.c_str());

    // Read solution.
    std::ifstream solution_file(solution_path);
    if (!solution_file.good()) {
        throw std::runtime_error(
                "Unable to open file \"" + std::string(solution_path) + "\".");
    }

    Solution solution(instance);
    VertexId number_of_vertices = -1;
    VertexId vertex_id = -1;
    solution_file >> number_of_vertices >> vertex_id;
    for (VertexId pos = 0; pos < number_of_vertices - 1; ++pos) {
        solution_file >> vertex_id;
        solution.add_vertex(distances, vertex_id);
    }

    // Remove temporary files.
    std::remove(instance_path);
    std::string sol_path = instance_path;
    sol_path = sol_path.substr(sol_path.find_last_of("/\\") + 1) + ".sol";
    std::remove(sol_path.c_str());
    std::remove(solution_path);
    std::remove(output_path);

    // Update output.
    algorithm_formatter.update_solution(solution, "final solution");
    algorithm_formatter.update_bound(solution.distance(), "final solution");

    algorithm_formatter.end();
    return output;
}

}
