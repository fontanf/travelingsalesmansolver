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
        << std::setw(12) << "Solution"
        << std::setw(12) << "Bound"
        << std::setw(24) << "Comment"
        << std::endl
        << std::setw(12) << "-----"
        << std::setw(12) << "--------"
        << std::setw(12) << "-----"
        << std::setw(24) << "-------"
        << std::endl;
    print(info, std::stringstream(""));
}

void Output::print(
        optimizationtools::Info& info,
        const std::stringstream& s) const
{
    double t = info.elapsed_time();
    std::streamsize precision = std::cout.precision();

    info.os()
        << std::setw(12) << std::fixed << std::setprecision(3) << t << std::defaultfloat << std::setprecision(precision)
        << std::setw(12) << solution.distance()
        << std::setw(12) << bound
        << std::setw(24) << s.str()
        << std::endl;

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
        // Update solution
        solution = solution_new;
        print(info, s);

        info.output->number_of_solutions++;
        double t = round(info.elapsed_time() * 10000) / 10000;
        std::string sol_str = "Solution" + std::to_string(info.output->number_of_solutions);
        info.add_to_json(sol_str, "Value", solution.distance());
        info.add_to_json(sol_str, "Time", t);
        info.add_to_json(sol_str, "String", s.str());
        if (!info.output->only_write_at_the_end) {
            info.write_json_output();
            solution.write(info.output->certificate_path);
        }
    }

    info.unlock();
}

Output& Output::algorithm_end(optimizationtools::Info& info)
{
    const Instance& instance = solution.instance();
    double t = round(info.elapsed_time() * 10000) / 10000;
    info.add_to_json("Solution", "Value", solution.distance());
    info.add_to_json("Solution", "Time", t);
    info.os()
        << std::endl
        << "Final statistics" << std::endl
        << "----------------" << std::endl
        << "# vertices:            " << solution.number_of_vertices() << " / " << instance.number_of_vertices() << " (" << (double)solution.number_of_vertices() / instance.number_of_vertices() * 100 << "%)" << std::endl
        << "Distance:              " << solution.distance() << std::endl
        << "Time (s):              " << t << std::endl;

    info.write_json_output();
    solution.write(info.output->certificate_path);
    return *this;
}
