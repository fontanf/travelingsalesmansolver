#include "travelingsalesmansolver/algorithms/lkh.hpp"

#include "travelingsalesmansolver/algorithm_formatter.hpp"

using namespace travelingsalesmansolver;

const LkhOutput travelingsalesmansolver::lkh(
        const Instance& instance,
        const LkhParameters& parameters)
{
    LkhOutput output(instance);
    AlgorithmFormatter algorithm_formatter(parameters, output);
    algorithm_formatter.start("LKH");
    algorithm_formatter.print_header();

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
    if (parameters.timer.time_limit() != std::numeric_limits<double>::infinity())
        parameters_file << "TIME_LIMIT = " << parameters.timer.remaining_time() << std::endl;
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
    if (!parameters.max_candidates.empty())
        parameters_file << "MAX_CANDIDATES = " << parameters.max_candidates << std::endl;

    // Candidate file.
    char candidate_path[L_tmpnam];
    tmpnam(candidate_path);
    parameters_file << "CANDIDATE_FILE  = " << candidate_path << std::endl;
    if (!parameters.candidate_file_content.empty()) {
        std::ofstream candidate_file(candidate_path);
        if (!candidate_file.good()) {
            throw std::runtime_error(
                    "Unable to open file \"" + std::string(candidate_path) + "\".");
        }
        candidate_file << parameters.candidate_file_content;
    }

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

    // Retrieve candidate file content.
    if (parameters.candidate_file_content.empty()) {
        // https://stackoverflow.com/questions/2602013/read-whole-ascii-file-into-c-stdstring
        std::ifstream candidate_file(candidate_path);
        candidate_file.seekg(0, std::ios::end);
        size_t size = candidate_file.tellg();
        output.candidate_file_content = std::string(size, ' ');
        candidate_file.seekg(0);
        candidate_file.read(&output.candidate_file_content[0], size);
    } else {
        output.candidate_file_content = parameters.candidate_file_content;
    }

    // Remove temporary files.
    std::remove(instance_path);
    std::remove(parameters_path);
    std::remove(solution_path);
    std::remove(output_path);
    std::remove(candidate_path);

    // Update output.
    algorithm_formatter.update_solution(solution, "Final solution");

    algorithm_formatter.end();
    return output;
}

std::vector<LkhCandidate> travelingsalesmansolver::read_candidates(
        std::string candidate_file_content)
{
    std::stringstream ss(candidate_file_content);
    VertexId number_of_vertices = -1;
    ss >> number_of_vertices;

    std::vector<LkhCandidate> candidates(number_of_vertices);

    VertexId vertex_id_tmp = -1;
    VertexId parent_id = -1;
    VertexId vertex_number_of_candidates = -1;
    VertexId vertex_id_2 = -1;
    double alpha = 0;
    for (VertexId vertex_id = 0; vertex_id < number_of_vertices; ++vertex_id) {
        ss >> vertex_id_tmp >> parent_id >> vertex_number_of_candidates;
        parent_id--;
        candidates[vertex_id].parent_id = parent_id;
        for (VertexId pos = 0; pos < vertex_number_of_candidates; ++pos) {
            ss >> vertex_id_2 >> alpha;
            vertex_id_2--;

            LkhCandidateEdge edge;
            edge.vertex_id = vertex_id_2;
            edge.alpha = alpha;
            candidates[vertex_id].edges.push_back(edge);
        }

    }

    return candidates;
}
