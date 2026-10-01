#include "travelingsalesmansolver/algorithms/lkh.hpp"
#include "travelingsalesmansolver/algorithms/concorde.hpp"
#include "travelingsalesmansolver/algorithms/greedy_edge.hpp"
#include "travelingsalesmansolver/algorithms/lin_kernighan.hpp"
#include "travelingsalesmansolver/candidates/alpha_nearness.hpp"
#include "travelingsalesmansolver/algorithms/eax.hpp"
#include "travelingsalesmansolver/algorithms/random_permutation.hpp"
#include "travelingsalesmansolver/algorithms/random_walk.hpp"
#include "travelingsalesmansolver/algorithms/two_opt.hpp"

#include <boost/program_options.hpp>

using namespace travelingsalesmansolver;

namespace po = boost::program_options;

void read_args(
        Parameters& parameters,
        const po::variables_map& vm)
{
    parameters.timer.set_sigint_handler();
    parameters.messages_to_stdout = true;
    if (vm.count("time-limit"))
        parameters.timer.set_time_limit(vm["time-limit"].as<double>());
    if (vm.count("verbosity-level"))
        parameters.verbosity_level = vm["verbosity-level"].as<int>();
    if (vm.count("log"))
        parameters.log_path = vm["log"].as<std::string>();
    parameters.log_to_stderr = vm.count("log-to-stderr");
    bool only_write_at_the_end = vm.count("only-write-at-the-end");
    if (!only_write_at_the_end) {
        std::string certificate_path = vm["certificate"].as<std::string>();
        std::string json_output_path = vm["output"].as<std::string>();
        parameters.new_solution_callback = [
            json_output_path,
            certificate_path](
                    const Output& output,
                    const std::string&)
        {
            output.write_json_output(json_output_path);
            output.solution.write(certificate_path);
        };
    }
}

/**
 * Compute the candidates given by the '--candidates' option ('nullptr'
 * output: the algorithm's default).
 */
std::unique_ptr<CandidateLists> compute_candidates(
        const Instance& instance,
        const po::variables_map& vm)
{
    if (!vm.count("candidates"))
        return nullptr;
    std::string type = vm["candidates"].as<std::string>();
    optimizationtools::Timer timer;
    std::unique_ptr<CandidateLists> candidates;
    if (type == "nearest-neighbor") {
        VertexId number_of_candidates = (vm.count("number-of-candidates"))?
            vm["number-of-candidates"].as<VertexId>(): 10;
        candidates = std::make_unique<CandidateLists>(
                nearest_neighbor_candidates(instance.distances(), number_of_candidates));
    } else if (type == "alpha-nearness") {
        AlphaNearnessParameters parameters;
        if (vm.count("number-of-candidates"))
            parameters.number_of_candidates = vm["number-of-candidates"].as<VertexId>();
        if (vm.count("ascent-graph-number-of-nearest-neighbors"))
            parameters.number_of_nearest_neighbors = vm["ascent-graph-number-of-nearest-neighbors"].as<VertexId>();
        if (vm.count("ascent-initial-period"))
            parameters.initial_period = vm["ascent-initial-period"].as<int64_t>();
        AlphaNearnessOutput output = alpha_nearness_candidates(instance.distances(), parameters);
        if (!vm.count("verbosity-level") || vm["verbosity-level"].as<int>() > 0) {
            std::cout
                << "Alpha-nearness candidates: "
                << "lower bound " << output.lower_bound
                << ", " << output.number_of_iterations << " iterations"
                << ", " << timer.elapsed_time() << " s" << std::endl;
        }
        candidates = std::make_unique<CandidateLists>(std::move(output.candidates));
    } else {
        throw std::invalid_argument(
                "Unknown candidates \"" + type + "\".");
    }
    return candidates;
}

Output run(
        const Instance& instance,
        const po::variables_map& vm)
{
    std::mt19937_64 generator(vm["seed"].as<Seed>());
    Solution solution(instance, vm["initial-solution"].as<std::string>());

    // Run algorithm.
    std::string algorithm = vm["algorithm"].as<std::string>();
    if (algorithm == "lkh") {
        LkhParameters parameters;
        if (vm.count("candidate-set-type"))
            parameters.candidate_set_type = vm["candidate-set-type"].as<std::string>();
        if (vm.count("initial-period"))
            parameters.initial_period = vm["initial-period"].as<std::string>();
        if (vm.count("runs"))
            parameters.runs = vm["runs"].as<std::string>();
        if (vm.count("max-trials"))
            parameters.max_trials = vm["max-trials"].as<std::string>();
        read_args(parameters, vm);
        return lkh(instance, generator, parameters);
    } else if (algorithm == "concorde") {
        Parameters parameters;
        read_args(parameters, vm);
        return concorde(instance, parameters);

    } else if (algorithm == "greedy-edge") {
        GreedyEdgeParameters parameters;
        read_args(parameters, vm);
        if (vm.count("number-of-candidates"))
            parameters.number_of_candidates = vm["number-of-candidates"].as<VertexId>();
        std::unique_ptr<CandidateLists> candidates = compute_candidates(instance, vm);
        return greedy_edge(instance, parameters, candidates.get());

    } else if (algorithm == "lin-kernighan") {
        LinKernighanParameters parameters;
        read_args(parameters, vm);
        if (vm.count("number-of-candidates"))
            parameters.number_of_candidates = vm["number-of-candidates"].as<VertexId>();
        if (vm.count("maximum-depth"))
            parameters.maximum_depth = vm["maximum-depth"].as<int>();
        if (vm.count("move-type"))
            parameters.move_type = vm["move-type"].as<int>();
        if (vm.count("restricted-search"))
            parameters.restricted_search = vm["restricted-search"].as<bool>();
        if (vm.count("non-sequential-moves"))
            parameters.non_sequential_moves = vm["non-sequential-moves"].as<bool>();
        if (vm.count("perturbation"))
            parameters.perturbation = vm["perturbation"].as<std::string>();
        if (vm.count("maximum-number-of-trials"))
            parameters.maximum_number_of_trials = vm["maximum-number-of-trials"].as<int64_t>();
        if (vm.count("candidates"))
            parameters.candidates = vm["candidates"].as<std::string>();
        if (vm.count("penalized-costs"))
            parameters.penalized_costs = vm["penalized-costs"].as<bool>();
        if (vm.count("ascent-initial-period"))
            parameters.ascent_initial_period = vm["ascent-initial-period"].as<int64_t>();
        return lin_kernighan(instance, generator, parameters);

    } else if (algorithm == "eax") {
        EaxParameters parameters;
        if (vm.count("population-size"))
            parameters.population_size = vm["population-size"].as<int>();
        if (vm.count("number-of-children"))
            parameters.number_of_children = vm["number-of-children"].as<int>();
        if (vm.count("initial-tours"))
            parameters.initial_tours = vm["initial-tours"].as<std::string>();
        if (vm.count("initial-local-search"))
            parameters.initial_local_search = vm["initial-local-search"].as<std::string>();
        if (vm.count("number-of-candidates"))
            parameters.number_of_candidates = vm["number-of-candidates"].as<VertexId>();
        read_args(parameters, vm);
        return eax(instance, generator, parameters);

    } else if (algorithm == "random-permutation") {
        RandomPermutationParameters parameters;
        read_args(parameters, vm);
        return random_permutation(instance, generator, parameters);

    } else if (algorithm == "random-walk") {
        RandomWalkParameters parameters;
        read_args(parameters, vm);
        if (vm.count("number-of-candidates"))
            parameters.number_of_candidates = vm["number-of-candidates"].as<VertexId>();
        if (vm.count("candidates"))
            parameters.candidates = vm["candidates"].as<std::string>();
        return random_walk(instance, generator, parameters);

    } else if (algorithm == "two-opt") {
        TwoOptParameters parameters;
        read_args(parameters, vm);
        if (vm.count("number-of-candidates"))
            parameters.number_of_candidates = vm["number-of-candidates"].as<VertexId>();
        std::unique_ptr<CandidateLists> candidates = compute_candidates(instance, vm);
        return two_opt(
                instance,
                generator,
                parameters,
                (solution.feasible())? &solution: nullptr,
                candidates.get());

    } else {
        throw std::invalid_argument(
                "Unknown algorithm \"" + algorithm + "\".");
    }
}

int main(int argc, char *argv[])
{
    // Parse program options
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "produce help message")
        ("algorithm,a", po::value<std::string>()->default_value("large-neighborhood-search"), "set algorithm")
        ("input,i", po::value<std::string>()->required(), "set input file (required)")
        ("format,f", po::value<std::string>()->default_value(""), "set input file format (default: standard)")
        ("output,o", po::value<std::string>()->default_value(""), "set JSON output file")
        ("initial-solution,", po::value<std::string>()->default_value(""), "")
        ("certificate,c", po::value<std::string>()->default_value(""), "set certificate file")
        ("goal,", po::value<Distance>(), "")
        ("seed,s", po::value<Seed>()->default_value(0), "set seed")
        ("time-limit,t", po::value<double>(), "set time limit in seconds")
        ("verbosity-level,v", po::value<int>(), "set verbosity level")
        ("only-write-at-the-end,e", "only write output and certificate files at the end")
        ("log,l", po::value<std::string>(), "set log file")
        ("log-to-stderr", "write log to stderr")

        ("candidate-set-type,", po::value<std::string>(), "set candidate set type")
        ("initial-period,", po::value<std::string>(), "set initial period")
        ("runs,", po::value<std::string>(), "set runs")
        ("max-trials,", po::value<std::string>(), "set max trials")

        ("candidates,", po::value<std::string>(), "set candidates: nearest-neighbor or alpha-nearness (greedy-edge, lin-kernighan, random-walk, two-opt)")
        ("number-of-candidates,", po::value<VertexId>(), "set number of candidates (greedy-edge, lin-kernighan, random-walk, two-opt, eax)")
        ("ascent-graph-number-of-nearest-neighbors,", po::value<VertexId>(), "set number of nearest neighbors of the graph of the ascent, -1 for the complete graph (alpha-nearness)")
        ("ascent-initial-period,", po::value<int64_t>(), "set initial period of the ascent (alpha-nearness)")
        ("penalized-costs,", po::value<bool>(), "guide the search with the costs penalized by the ascent of the alpha-nearness candidates (lin-kernighan)")
        ("maximum-number-of-trials,", po::value<int64_t>(), "set maximum number of trials or kicks (lin-kernighan)")
        ("move-type,", po::value<int>(), "set move type: 3 to 5 (lin-kernighan)")
        ("perturbation,", po::value<std::string>(), "set perturbation: walks, double-bridge or segment-swap (lin-kernighan)")
        ("restricted-search,", po::value<bool>(), "set restricted search (lin-kernighan)")
        ("non-sequential-moves,", po::value<bool>(), "set non-sequential moves (lin-kernighan)")
        ("maximum-depth,", po::value<int>(), "set maximum number of steps in a chain (lin-kernighan)")
        ("population-size,", po::value<int>(), "set population size (eax)")
        ("number-of-children,", po::value<int>(), "set number of children per generation (eax)")
        ("initial-tours,", po::value<std::string>(), "set initial tours of the population: random-permutation or random-walk (eax)")
        ("initial-local-search,", po::value<std::string>(), "set local search of the initial tours: 2-opt, lin-kernighan or none (eax)")
        ;
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    if (vm.count("help")) {
        std::cout << desc << std::endl;;
        return 1;
    }
    try {
        po::notify(vm);
    } catch (const po::required_option& e) {
        std::cout << desc << std::endl;;
        return 1;
    }

    // Build instance.
    const Instance instance(
            vm["input"].as<std::string>(),
            vm["format"].as<std::string>());

    // Run.
    Output output = run(instance, vm);

    // Write outputs.
    std::string certificate_path = vm["certificate"].as<std::string>();
    std::string json_output_path = vm["output"].as<std::string>();
    output.write_json_output(json_output_path);
    output.solution.write(certificate_path);

    return 0;
}
