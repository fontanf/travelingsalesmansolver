#include "travelingsalesmansolver/algorithms/algorithms.hpp"

#include <boost/program_options.hpp>

using namespace travelingsalesmansolver;
namespace po = boost::program_options;

Output travelingsalesmansolver::run(
        std::string algorithm,
        const Instance& instance,
        const Solution& initial_solution,
        std::mt19937_64& generator,
        optimizationtools::Info info)
{
    (void)initial_solution;
    (void)generator;
    std::vector<std::string> algorithm_args = po::split_unix(algorithm);
    std::vector<char*> algorithm_argv;
    for (Counter i = 0; i < (Counter)algorithm_args.size(); ++i)
        algorithm_argv.push_back(const_cast<char*>(algorithm_args[i].c_str()));

    if (algorithm.empty() || algorithm_args[0].empty()) {
        throw std::invalid_argument("Missing algorithm.");

    } else if (algorithm_args[0] == "concorde") {
        return concorde(instance, info);

    } else if (algorithm_args[0] == "lkh") {
        LkhOptionalParameters parameters;
        parameters.info = info;
        return lkh(instance, parameters);

    } else {
        throw std::invalid_argument(
                "Unknown algorithm \"" + algorithm_args[0] + "\".");
    }
}

