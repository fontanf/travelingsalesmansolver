#include "travelingsalesmansolver/distances/ball_tree.hpp"

#include <numeric>
#include <algorithm>

using namespace travelingsalesmansolver;

BallTree::BallTree(const Distances& distances):
    distances_(distances),
    neighbors_tmp_(distances.number_of_vertices())
{
    build_tree();
}

void BallTree::build_tree()
{
    return FUNCTION_WITH_DISTANCES_0(
            (this->BallTree::build_tree),
            distances_);
}

template <typename T>
void BallTree::build_tree(
        const T& distances)
{
    Node root;
    tree_.push_back(root);

    std::vector<StackElement> stack;
    StackElement stack_initial_element;
    stack_initial_element.node_id = 0;
    stack_initial_element.vertex_ids = std::vector<VertexId>(distances_.number_of_vertices());
    std::iota(stack_initial_element.vertex_ids.begin(), stack_initial_element.vertex_ids.end(), 0);
    stack.push_back(stack_initial_element);

    std::vector<Distance> distances_from_center(distances_.number_of_vertices(), 0);;
    while (!stack.empty()) {
        StackElement stack_element = stack.back();
        stack.pop_back();

        NodeId node_id = stack_element.node_id;
        Node& node = tree_[node_id];

        // Choose a center.
        //std::cout
        //    << "node_id " << node_id << std::endl
        //    << " stack_element.vertex_ids.size() " << stack_element.vertex_ids.size() << std::endl;
        node.center_vertex_id = stack_element.vertex_ids[0];
        //std::cout << " center_vertex_id " << node.center_vertex_id << std::endl;

        // Compute the distances to the center.
        distances_from_center.clear();
        for (VertexId pos = 1;
                pos < stack_element.vertex_ids.size();
                ++pos) {
            VertexId vertex_id = stack_element.vertex_ids[pos];
            Distance distance = distances.distance(
                    node.center_vertex_id,
                    vertex_id);
            //std::cout << "pos " << pos
            //    << " vertex_id " << vertex_id
            //    << " distance " << distance
            //    << std::endl;
            distances_from_center.push_back(distance);
            if (node.radius < distance)
                node.radius = distance;
        }
        //std::cout << " radius " << node.radius << std::endl;

        if (node.radius == 0) {
            node.vertex_ids = stack_element.vertex_ids;
            continue;
        }

        // Find the median.
        VertexId middle = distances_from_center.size() / 2;
        std::nth_element(
                distances_from_center.begin(),
                distances_from_center.begin() + middle,
                distances_from_center.end());
        //std::cout << " middle " << middle << std::endl;
        node.median_distance = distances_from_center[middle];
        //std::cout << " node.median_distance " << node.median_distance << std::endl;

        // Compute left/right/bottom/top shapes.
        node.lesser_child_id = tree_.size();
        node.greater_child_id = tree_.size() + 1;
        //std::cout << " node.lesser_child_id " << node.lesser_child_id << std::endl;
        //std::cout << " node.greater_child_id " << node.greater_child_id << std::endl;
        Node child_lesser;
        Node child_greater;
        StackElement stack_element_lesser;
        StackElement stack_element_greater;
        stack_element_lesser.vertex_ids = {node.center_vertex_id};
        stack_element_lesser.node_id = node.lesser_child_id;
        stack_element_greater.node_id = node.greater_child_id;
        for (VertexId pos = 1;
                pos < stack_element.vertex_ids.size();
                ++pos) {
            VertexId vertex_id = stack_element.vertex_ids[pos];
            Distance distance = distances.distance(
                        node.center_vertex_id,
                        vertex_id);
            if (node.median_distance == 0) {
                if (distance == 0) {
                    stack_element_lesser.vertex_ids.push_back(vertex_id);
                } else {
                    stack_element_greater.vertex_ids.push_back(vertex_id);
                }
            } else {
                if (distance < node.median_distance) {
                    stack_element_lesser.vertex_ids.push_back(vertex_id);
                } else {
                    stack_element_greater.vertex_ids.push_back(vertex_id);
                }
            }
        }

        tree_.push_back(child_lesser);
        tree_.push_back(child_greater);
        stack.push_back(stack_element_lesser);
        stack.push_back(stack_element_greater);
    }
}

std::vector<VertexId> BallTree::nearest_neighbors(
        VertexId vertex_id,
        VertexId number_of_neighbors)
{
    return FUNCTION_WITH_DISTANCES(
            (this->BallTree::nearest_neighbors),
            distances_,
            vertex_id,
            number_of_neighbors);
}

template <typename T>
std::vector<VertexId> BallTree::nearest_neighbors(
        const T& distances,
        VertexId vertex_id,
        VertexId number_of_neighbors)
{
    //std::cout << "nearest_neighbors"
    //    << " vertex_id " << vertex_id
    //    << " number_of_neighbors " << number_of_neighbors
    //    << std::endl;
    neighbors_tmp_.clear();

    if (distances_.number_of_vertices() == 0)
        return {};

    std::vector<NodeId> stack = {0};

    while (!stack.empty()) {

        NodeId node_id = stack.back();
        stack.pop_back();
        const Node& node = tree_[node_id];
        //std::cout << "node_id " << node_id << std::endl;

        // Compute the minimum distance between vertex 'vertex_id' and any point
        // in the current node.
        if (neighbors_tmp_.size() == number_of_neighbors) {
            Distance distance_min = distances.distance(
                    vertex_id,
                    node.center_vertex_id)
                - node.radius;
            //std::cout << " distance_min " << distance_min << std::endl;
            if (distance_min >= -neighbors_tmp_.top().second)
                continue;
        }

        if (!node.vertex_ids.empty()) {
            //std::cout << " leaf " << std::endl;
            // The node is a leaf.
            Distance distance = distances.distance(
                    vertex_id,
                    node.center_vertex_id);
            for (VertexId neighbor_id: node.vertex_ids) {
                if (neighbor_id == vertex_id)
                    continue;
                //std::cout << " add neighbor_id " << neighbor_id
                //    << " distance " << distance
                //    << std::endl;
                neighbors_tmp_.update_key(neighbor_id, -distance);
                if (neighbors_tmp_.size() > number_of_neighbors)
                    neighbors_tmp_.pop();
            }
        } else {
            // Recursion.
            stack.push_back(node.lesser_child_id);
            stack.push_back(node.greater_child_id);
        }
    }

    // Convert heap to vector.
    std::vector<VertexId> neighbor_ids;
    for (VertexId pos = neighbors_tmp_.size() - 1;
            pos >= 0;
            --pos) {
        VertexId vertex_id = neighbors_tmp_.top(pos).first;
        neighbor_ids.push_back(vertex_id);
    }

    return neighbor_ids;
}
