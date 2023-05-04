#include "travelingsalesmansolver/algorithms/lkh.hpp"

using namespace travelingsalesmansolver;

Output travelingsalesmansolver::lkh(
        const Instance& instance,
        LkhOptionalParameters parameters)
{
    init_display(instance, parameters.info);
    parameters.info.os()
        << "Algorithm" << std::endl
        << "---------" << std::endl
        << "LKH" << std::endl
        << std::endl;

    Output output(instance, parameters.info);

    // Write instance file.
    char instance_path[L_tmpnam];
    tmpnam(instance_path);
    instance.write(instance_path);

    // Write parameters file.
    char parameters_path[L_tmpnam];
    tmpnam(parameters_path);

    std::ofstream parameters_file(parameters_path);
    if (!parameters_file.good()) {
        throw std::runtime_error(
                "Unable to open file \"" + std::string(parameters_path) + "\".");
    }

    parameters_file << "PROBLEM_FILE = " << instance_path << std::endl;
    char solution_path[L_tmpnam];
    tmpnam(solution_path);
    parameters_file << "OUTPUT_TOUR_FILE = " << solution_path << std::endl;
    if (parameters.info.time_limit != std::numeric_limits<double>::infinity())
        parameters_file << "TIME_LIMIT = " << parameters.info.remaining_time() << std::endl;
    if (!parameters.candidate_set_type.empty())
        parameters_file << "CANDIDATE_SET_TYPE = " << parameters.candidate_set_type << std::endl;
    if (!parameters.initial_period.empty())
        parameters_file << "INITIAL_PERIOD = " << parameters.initial_period << std::endl;
    if (!parameters.runs.empty())
        parameters_file << "RUNS = " << parameters.runs << std::endl;
    if (!parameters.max_trials.empty())
        parameters_file << "MAX_TRIALS = " << parameters.max_trials << std::endl;
    if (!parameters.seed.empty())
        parameters_file << "SEED = " << parameters.seed << std::endl;

    // Run.
    char output_path[L_tmpnam];
    tmpnam(output_path);
    std::string command = (
            "LKH"
            " \"" + std::string(parameters_path) + "\""
            + " > \"" + std::string(output_path) + "\"");
    std::system(command.c_str());

    // Read solution.
    std::ifstream solution_file(solution_path);
    if (!solution_file.good()) {
        throw std::runtime_error(
                "Unable to open file \"" + std::string(solution_path) + "\".");
    }

    Solution solution(instance);
    std::string tmp;
    std::vector<std::string> line;
    while (getline(solution_file, tmp)) {
        line = optimizationtools::split(tmp, ' ');
        if (line.size() == 0) {
        } else if (tmp.rfind("TOUR_SECTION", 0) == 0) {
            VertexId vertex_id = -1;
            solution_file >> vertex_id;
            for (;;) {
                solution_file >> vertex_id;
                if (vertex_id == -1)
                    break;
                solution.add_vertex(vertex_id - 1);
            }
        }
    }

    // Remove temporary files.
    std::remove(instance_path);
    std::remove(parameters_path);
    std::remove(solution_path);
    std::remove(output_path);

    // Update output.
    std::stringstream ss;
    ss << "final solution";
    output.update_solution(solution, ss, parameters.info);

    return output.algorithm_end(parameters.info);
}
