#ifndef GRAPH_SPARSE_GRAPH_HPP
#define GRAPH_SPARSE_GRAPH_HPP

#include <algorithm>   // std::find_if, std::swap
#include <concepts>    // concepts
#include <cstddef>     // std::size_t
#include <fstream>     // std::ofstream
#include <functional>  // std::function
#include <map>         // std::map
#include <queue>       // std::queue
#include <stack>       // std::stack
#include <stdexcept>   // std::invalid_argument
#include <string>      // std::string
#include <vector>      // std::vector

#include <iostream>

#include "graph.hpp"

namespace graph {
template <typename NodeType, typename EdgeType>
requires HasId<NodeType>&& HasSourceAndTarget<EdgeType> class SparseGraph : public Graph<NodeType, EdgeType> {
 public:
  SparseGraph() = default;

  SparseGraph(const SparseGraph&) = default;
  SparseGraph(SparseGraph&&) noexcept = default;

  SparseGraph& operator=(const SparseGraph&) = default;
  SparseGraph& operator=(SparseGraph&&) noexcept = default;

  ~SparseGraph() = default;

  size_t addNode() override {
    nodes_.emplace_back(NodeType{.id = nodes_.size()});
    adjacency_list_.emplace_back();

    return nodes_.size() - 1;
  }

  void addEdge(size_t node1_id, size_t node2_id) override {
    EdgeType edge{};
    if (node1_id < node2_id) {
      edge.source = node1_id;
      edge.target = node2_id;
    }
    else {
      edge.source = node2_id;
      edge.target = node1_id;
    }

    adjacency_list_[node1_id].emplace_back(edges_.size());
    adjacency_list_[node2_id].emplace_back(edges_.size());
    edges_.emplace_back(edge);
  }

  [[nodiscard]] bool edgeExists(size_t node1_id, size_t node2_id) const override {
    const size_t source = node1_id < node2_id ? node1_id : node2_id;
    const size_t target = node1_id < node2_id ? node2_id : node1_id;

    return std::ranges::any_of(adjacency_list_[source],
                               [&](size_t index) { return edges_[index].source == source && edges_[index].target == target; });
  }

  [[nodiscard]] NodeType& getNode(size_t node_id) override { return nodes_[node_id]; }

  [[nodiscard]] NodeType& getLastAddedNode() override { return nodes_.back(); }

  [[nodiscard]] std::vector<NodeType> getNodes() const override { return nodes_; }

  [[nodiscard]] EdgeType& getEdge(size_t source, size_t target) override {
    auto iter = std::find_if(adjacency_list_[source].begin(), adjacency_list_[source].end(), [&](size_t index) {
      return (edges_[index].source == source && edges_[index].target == target) ||
             (edges_[index].source == target && edges_[index].target == source);
    });

    if (iter != adjacency_list_[source].end()) {
      return edges_[*iter];
    }

    throw std::invalid_argument("The edge (" + std::to_string(source) + "," + std::to_string(target) + ") does not exist!");
  }

  [[nodiscard]] EdgeType& getLastAddedEdge() override { return edges_.back(); }

  [[nodiscard]] std::vector<EdgeType> getEdges() const override { return edges_; }

  [[nodiscard]] size_t getNodesNumber() const override { return nodes_.size(); }

  [[nodiscard]] size_t getEdgesNumber() const override { return edges_.size(); }

  [[nodiscard]] float getDensity() const override {
    return static_cast<float>(2 * getEdgesNumber()) / (getNodesNumber() * (getNodesNumber() - 1));
  }

  [[nodiscard]] std::vector<NodeType> getNeighbours(size_t node_id) const override {
    std::vector<NodeType> neighbours;
    for (const size_t index : adjacency_list_[node_id]) {
      if (edges_[index].source == node_id) {
        neighbours.emplace_back(nodes_[edges_[index].target]);
      }
      else {
        neighbours.emplace_back(nodes_[edges_[index].source]);
      }
    }

    return neighbours;
  }

  [[nodiscard]] std::vector<EdgeType> getAdjacentEdges(size_t node_id) const override {
    std::vector<EdgeType> edges;
    for (const EdgeType& edge : edges_) {
      if (edge.source == node_id || edge.target == node_id) {
        edges.emplace_back(edge);
      }
    }

    return edges;
  }

  [[nodiscard]] size_t getDegree(size_t node_id) const override { return adjacency_list_[node_id].size(); }

  [[nodiscard]] std::map<size_t, size_t> getDegreesHistogram() const override {
    std::map<size_t, size_t> histogram;
    for (const auto& node_edges : adjacency_list_) {
      if (histogram.find(node_edges.size()) == histogram.end()) {
        histogram[node_edges.size()] = 1;
      }
      else {
        histogram[node_edges.size()] += 1;
      }
    }

    return histogram;
  }

  void dfs(size_t start_node_id, std::function<void(const NodeType&)> callback) const override {
    std::vector<bool> visited(getNodesNumber(), false);
    std::stack<size_t> stack;

    stack.push(start_node_id);
    while (!stack.empty()) {
      const size_t node_id = stack.top();
      stack.pop();

      if (!visited[node_id]) {
        callback(nodes_[node_id]);
        visited[node_id] = true;

        for (const size_t neighbor_index : adjacency_list_.at(node_id)) {
          const EdgeType& neighbor = edges_[neighbor_index];
          if (neighbor_index == neighbor.source) {
            if (!visited[neighbor.target]) {
              stack.push(neighbor.target);
            }
          }
          else {
            if (!visited[neighbor.source]) {
              stack.push(neighbor.source);
            }
          }
        }
      }
    }
  }

  void bfs(size_t start_node_id, std::function<void(const NodeType&)> callback) const override {
    std::vector<bool> visited(getNodesNumber(), false);
    std::queue<size_t> queue;

    queue.push(start_node_id);
    while (!queue.empty()) {
      const size_t node_id = queue.front();
      queue.pop();

      if (!visited[node_id]) {
        callback(nodes_[node_id]);
        visited[node_id] = true;

        for (const size_t neighbor_index : adjacency_list_.at(node_id)) {
          const EdgeType& neighbor = edges_[neighbor_index];
          if (neighbor_index == neighbor.source) {
            if (!visited[neighbor.target]) {
              queue.push(neighbor.target);
            }
          }
          else {
            if (!visited[neighbor.source]) {
              queue.push(neighbor.source);
            }
          }
        }
      }
    }
  }

  void saveToFile(const std::string& filename) const override {
    std::ofstream file{filename};

    file << "{" << std::endl << "  \"nodes\": [" << std::endl;
    for (size_t i = 0; i < getNodesNumber(); ++i) {
      file << "    {\"id\": " << nodes_[i].id << "}," << std::endl;
    }
    file << "    { }" << std::endl;

    file << "  ]," << std::endl << "  \"edges\": [" << std::endl;
    for (size_t i = 0; i < getEdgesNumber(); ++i) {
      file << "    {\"source\": " << edges_[i].source << ", \"target\": " << edges_[i].target << "}," << std::endl;
    }
    file << "    { }" << std::endl;

    file << "  ]" << std::endl << "}" << std::endl;
  }

  void saveToFile(const std::string& filename, std::function<void(std::ofstream& file, const NodeType&)> vertexCallback,
                  std::function<void(std::ofstream& file, const EdgeType&)> edgeCallback) const override {
    std::ofstream file{filename};

    file << "{" << std::endl << "  \"nodes\": [" << std::endl;
    for (size_t i = 0; i < getNodesNumber(); ++i) {
      file << "    {\"id\": " << nodes_[i].id;
      vertexCallback(file, nodes_[i]);
      file << "}," << std::endl;
    }
    file << "    { }" << std::endl;

    file << "  ]," << std::endl << "  \"edges\": [" << std::endl;
    for (size_t i = 0; i < getEdgesNumber(); ++i) {
      file << "    {\"source\": " << edges_[i].source << ", \"target\": " << edges_[i].target;
      edgeCallback(file, edges_[i]);
      file << "}," << std::endl;
    }
    file << "    { }" << std::endl;

    file << "  ]" << std::endl << "}" << std::endl;
  }

  int loadFromFile(const std::string& filename) override {
    try {
      std::ifstream input_file(filename);
      json graph_data;
      input_file >> graph_data;

      for (const auto& node : graph_data["nodes"]) {
        addNode();
      }

      for (const auto& edge : graph_data["edges"]) {
        addEdge(edge["source"], edge["target"]);
      }
    } catch (const std::exception& e) {
      std::cerr << e.what() << std::endl;
      return 1;
    }

    return 0;
  }

 private:
  std::vector<NodeType> nodes_;
  std::vector<EdgeType> edges_;
  std::vector<std::vector<size_t>> adjacency_list_;
};
}  // namespace graph

#endif  // GRAPH_SPARSE_GRAPH_HPP
