#include "travelingsalesmansolver/solution.hpp"

using namespace travelingsalesmansolver;

Solution::Solution(const Instance& instance):
    instance_(&instance),
    vertices_is_visited_(instance.number_of_vertices(), false)
{
    // Add initial vertex.
    vertex_ids_.push_back(0);
    vertices_is_visited_[0] = true;
}

Solution::Solution(
        const Instance& instance,
        std::string certificate_path):
    Solution(instance)
{
    if (certificate_path.empty())
        return;
    std::ifstream file(certificate_path);
    if (!file.good()) {
        throw std::runtime_error(
                "Unable to open file \"" + certificate_path + "\".");
    }

    // TODO
}

bool Solution::feasible() const
{
    return (number_of_vertices() == instance().number_of_vertices());
}

void Solution::add_vertex(VertexId vertex_id)
{
    // Check that the vertex has not already been visited.
    if (vertices_is_visited_[vertex_id]) {
        throw std::runtime_error("");  // TODO
    }

    VertexId vertex_id_prev = vertex_ids_.back();
    vertex_ids_.push_back(vertex_id);
    vertices_is_visited_[vertex_id] = true;
    distance_cur_ += instance().distance(vertex_id_prev, vertex_id);
    distance_ = distance_cur_ + instance().distance(vertex_id, 0);
}

std::ostream& Solution::print(
        std::ostream& os,
        int verbose) const
{
    if (verbose >= 1) {
        os
            << "Number of vertices:  " << optimizationtools::Ratio<VertexId>(number_of_vertices(), instance().number_of_vertices()) << std::endl
            << "Feasible:            " << feasible() << std::endl
            << "Distance:            " << distance() << std::endl
            ;
    }

    if (verbose >= 2) {
        os << std::endl
            << std::setw(12) << "Vertex"
            << std::endl
            << std::setw(12) << "------"
            << std::endl;
        for (VertexId vertex_id: vertex_ids_) {
            os
                << std::setw(12) << vertex_id
                << std::endl;
        }
    }

    return os;
}

void Solution::write(std::string certificate_path) const
{
    if (certificate_path.empty())
        return;
    std::ofstream file(certificate_path);
    if (!file.good()) {
        throw std::runtime_error(
                "Unable to open file \"" + certificate_path + "\".");
    }

    for (VertexId vertex_id: vertex_ids_)
        file << " " << vertex_id + 1;
}

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////// Output ////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

Output::Output(
        const Instance& instance,
        optimizationtools::Info& info):
    solution(instance)
{
    info.os()
        << std::setw(12) << "T (s)"
        << std::setw(12) << "UB"
        << std::setw(12) << "LB"
        << std::setw(12) << "GAP"
        << std::setw(12) << "GAP (%)"
        << std::setw(24) << "Comment"
        << std::endl
        << std::setw(12) << "-----"
        << std::setw(12) << "--"
        << std::setw(12) << "--"
        << std::setw(12) << "---"
        << std::setw(12) << "-------"
        << std::setw(24) << "-------"
        << std::endl;
    print(info, std::stringstream(""));
}

void Output::print(
        optimizationtools::Info& info,
        const std::stringstream& s) const
{
    std::string solution_value = optimizationtools::solution_value(
            optimizationtools::ObjectiveDirection::Minimize,
            solution.feasible(),
            solution.distance());
    double absolute_optimality_gap = optimizationtools::absolute_optimality_gap(
            optimizationtools::ObjectiveDirection::Minimize,
            solution.feasible(),
            solution.distance(),
            bound);
    double relative_optimality_gap = optimizationtools::relative_optimality_gap(
            optimizationtools::ObjectiveDirection::Minimize,
            solution.feasible(),
            solution.distance(),
            bound);
    double t = info.elapsed_time();
    std::streamsize precision = std::cout.precision();

    info.os()
        << std::setw(12) << std::fixed << std::setprecision(3) << t << std::defaultfloat << std::setprecision(precision)
        << std::setw(12) << solution_value
        << std::setw(12) << bound
        << std::setw(12) << absolute_optimality_gap
        << std::setw(12) << std::fixed << std::setprecision(2) << relative_optimality_gap * 100 << std::defaultfloat << std::setprecision(precision)
        << std::setw(24) << s.str() << std::endl;

    if (!info.output->only_write_at_the_end)
        info.write_json_output();
}

void Output::update_solution(
        const Solution& solution_new,
        const std::stringstream& s,
        optimizationtools::Info& info)
{
    info.lock();

    if (solution_new.feasible()
            && (!solution.feasible()
                || solution.distance() > solution_new.distance())) {
        solution = solution_new;
        print(info, s);

        std::string solution_value = optimizationtools::solution_value(
                optimizationtools::ObjectiveDirection::Minimize,
                solution.feasible(),
                solution.distance());
        double t = info.elapsed_time();

        info.output->number_of_solutions++;
        std::string sol_str = "Solution" + std::to_string(info.output->number_of_solutions);
        info.add_to_json(sol_str, "Value", solution_value);
        info.add_to_json(sol_str, "Time", t);
        info.add_to_json(sol_str, "String", s.str());
        if (!info.output->only_write_at_the_end) {
            info.write_json_output();
            solution.write(info.output->certificate_path);
        }
    }

    info.unlock();
}

void Output::update_bound(
        Distance bound_new,
        const std::stringstream& s,
        optimizationtools::Info& info)
{
    if (bound >= bound_new)
        return;

    info.lock();

    if (bound < bound_new) {
        bound = bound_new;
        print(info, s);

        double t = info.elapsed_time();
        info.output->number_of_bounds++;
        std::string sol_str = "Bound" + std::to_string(info.output->number_of_bounds);
        info.add_to_json(sol_str, "Bound", bound);
        info.add_to_json(sol_str, "Time", t);
        info.add_to_json(sol_str, "String", s.str());
        if (!info.output->only_write_at_the_end)
            solution.write(info.output->certificate_path);
    }

    info.unlock();
}

Output& Output::algorithm_end(optimizationtools::Info& info)
{
    std::string solution_value = optimizationtools::solution_value(
            optimizationtools::ObjectiveDirection::Minimize,
            solution.feasible(),
            solution.distance());
    double absolute_optimality_gap = optimizationtools::absolute_optimality_gap(
            optimizationtools::ObjectiveDirection::Minimize,
            solution.feasible(),
            solution.distance(),
            bound);
    double relative_optimality_gap = optimizationtools::relative_optimality_gap(
            optimizationtools::ObjectiveDirection::Minimize,
            solution.feasible(),
            solution.distance(),
            bound);
    time = info.elapsed_time();

    info.add_to_json("Solution", "Value", solution_value);
    info.add_to_json("Bound", "Value", bound);
    info.add_to_json("Solution", "Time", time);
    info.add_to_json("Bound", "Time", time);
    info.os()
        << std::endl
        << "Final statistics" << std::endl
        << "----------------" << std::endl
        << "Value:                         " << solution_value << std::endl
        << "Bound:                         " << bound << std::endl
        << "Absolute optimality gap:       " << absolute_optimality_gap << std::endl
        << "Relative optimality gap (%):   " << relative_optimality_gap * 100 << std::endl
        << "Time (s):                      " << time << std::endl
        ;
    print_statistics(info);
    info.os() << std::endl
        << "Solution" << std::endl
        << "--------" << std::endl ;
    solution.print(info.os(), info.verbosity_level());

    info.write_json_output();
    solution.write(info.output->certificate_path);
    return *this;
}
