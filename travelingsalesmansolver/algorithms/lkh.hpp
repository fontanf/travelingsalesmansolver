#pragma once

#include "travelingsalesmansolver/solution.hpp"
#include "travelingsalesmansolver/algorithm_formatter.hpp"

#include <iomanip>

namespace travelingsalesmansolver
{

struct LkhParameters: Parameters
{
    /** CANDIDATE_SET_TYPE */
    std::string candidate_set_type;

    /** INITIAL_PERIOD */
    std::string initial_period;

    /** RUNS */
    std::string runs;

    /** MAX_TRIALS */
    std::string max_trials;

    /** SEED */
    std::string seed;

    /** CANDIDATE_FILE content. */
    std::string candidate_file_content;

    /** MAX_CANDIDATES. */
    std::string max_candidates;


    virtual int format_width() const override { return 37; }

    virtual void format(std::ostream& os) const override
    {
        Parameters::format(os);
        int width = format_width();
        os
            << std::setw(width) << std::left << "Candidate set type: " << candidate_set_type << std::endl
            << std::setw(width) << std::left << "Initial period:  " << initial_period << std::endl
            << std::setw(width) << std::left << "Runs:  " << runs << std::endl
            << std::setw(width) << std::left << "Max trials:  " << max_trials << std::endl
            << std::setw(width) << std::left << "Seed:  " << seed << std::endl
            << std::setw(width) << std::left << "Has candidate flle content:  " << (!candidate_file_content.empty()) << std::endl
            << std::setw(width) << std::left << "Max candidates:  " << max_candidates << std::endl
            ;
    }

    virtual nlohmann::json to_json() const override
    {
        nlohmann::json json = Parameters::to_json();
        json.merge_patch({
                {"CandidateSetType", candidate_set_type},
                {"InitialPeriod", initial_period},
                {"Runs", runs},
                {"MaxTrials", max_trials},
                {"Seed", seed},
                {"HasCandidateFileContent", (!candidate_file_content.empty())},
                {"MaxCandidates", max_candidates},
                });
        return json;
    }
};

struct LkhOutput: Output
{
    LkhOutput(const Instance& instance):
        Output(instance) { }


    std::string candidate_file_content;
};

const LkhOutput lkh(
        const Instance& instance,
        const LkhParameters& parameters = {});


template <typename Distances>
const LkhOutput lkh(
        const Distances& distances,
        const Instance& instance,
        const LkhParameters& parameters = {});

struct LkhCandidateEdge
{
    /** Other end of the edge. */
    VertexId vertex_id;

    /** Alpha value. */
    double alpha;
};

struct LkhCandidate
{
    /**
     * Parent of the vertex in the minimum spanning tree.
     *
     * -1 if the vertex has no parent.
     */
    VertexId parent_id = -1;

    /** Candidate edges. */
    std::vector<LkhCandidateEdge> edges;
};

std::vector<LkhCandidate> read_candidates(
        const std::string& candidate_file_content);

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <typename Distances>
const LkhOutput lkh(
        const Distances& distances,
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
                solution.add_vertex(distances, vertex_id - 1);
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

}
