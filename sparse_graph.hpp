#ifndef GRAPH_SPARSE_GRAPH_HPP
#define GRAPH_SPARSE_GRAPH_HPP

#include <algorithm>    // std::find
#include <concepts>     // concepts
#include <cstddef>      // std::size_t
#include <functional>   // std::function
#include <map>          // std::map
#include <queue>        // std::queue
#include <stack>        // std::stack
#include <vector>       // std::vector

#include <iostream>

#include "graph.hpp"

namespace graph {
template<typename NodeType, typename EdgeType> requires HasId<NodeType> && HasSourceAndTarget<EdgeType>
class SparseGraph : public Graph<NodeType, EdgeType> {
 public:
  SparseGraph() = default;

  SparseGraph(const SparseGraph &) = default;
  SparseGraph(SparseGraph &&) noexcept = default;

  SparseGraph &operator=(const SparseGraph &) = default;
  SparseGraph &operator=(SparseGraph &&) noexcept = default;

  ~SparseGraph() = default;

  void addNode(NodeType &node) override {
	node.id = nodes_.size();
	nodes_.emplace_back(node);
	adjacency_list_.emplace_back(std::vector<EdgeType>());
  }

  bool addEdge(EdgeType &edge) override {
	if (edge.source > edge.target) {
	  std::swap(edge.source, edge.target);
	}

	if (std::find_if(adjacency_list_[edge.source].begin(), adjacency_list_[edge.source].end(),
					 [&edge](const EdgeType &e) { return e.target == edge.target && e.source == edge.source; })
		!= adjacency_list_[edge.source].end()) {
	  return false;
	}

	adjacency_list_[edge.source].emplace_back(edge);
	adjacency_list_[edge.target].emplace_back(edge);

	return true;
  }

  [[nodiscard]]
  NodeType getNode(size_t node_id) const override {
	return nodes_[node_id];
  }

  [[nodiscard]]
  EdgeType getEdge(size_t source, size_t target) const override {
	return adjacency_list_[source][target];
  }

  [[nodiscard]]
  size_t getNodesNumber() const override {
	return nodes_.size();
  }

  [[nodiscard]]
  size_t getEdgesNumber() const override {
	size_t adjacency_list_number = 0;
	for (const auto &node_edges : adjacency_list_) {
	  adjacency_list_number += node_edges.size();
	}

	return adjacency_list_number / 2;
  }

  [[nodiscard]]
  float getDensity() const override {
	return static_cast<float>(2 * getEdgesNumber()) / (getNodesNumber() * (getNodesNumber() - 1));
  }

  [[nodiscard]]
  std::vector<NodeType> getNeighbours(size_t node_id) const override {
	std::vector<NodeType> neighbours;
	for (const auto &edge : adjacency_list_[node_id]) {
	  neighbours.emplace_back(nodes_[edge.target]);
	}

	return neighbours;
  }

  [[nodiscard]]
  std::vector<EdgeType> getAdjacentEdges(size_t node_id) const override {
	return adjacency_list_[node_id];
  }

  [[nodiscard]]
  size_t getDegree(size_t node_id) const override {
	return adjacency_list_[node_id].size();
  }

  [[nodiscard]]
  std::map<size_t, size_t> getDegreesHistogram() const override {
	std::map<size_t, size_t> histogram;
	for (const auto &node_edges : adjacency_list_) {
	  if (histogram.find(node_edges.size()) == histogram.end()) {
		histogram[node_edges.size()] = 1;
	  } else {
		histogram[node_edges.size()] += 1;
	  }
	}

	return histogram;
  }

  void dfs(size_t start_node_id, std::function<void(const NodeType &)> callback) const override {
	std::vector<bool> visited(getNodesNumber(), false);
	std::stack<size_t> stack;

	stack.push(start_node_id);
	while (!stack.empty()) {
	  size_t node_id = stack.top();
	  stack.pop();

	  if (!visited[node_id]) {
		callback(nodes_[node_id]);
		visited[node_id] = true;

		for (const EdgeType &neighbor : adjacency_list_.at(node_id)) {
		  if (!visited[neighbor.target]) {
			stack.push(neighbor.target);
		  }
		}
	  }
	}
  }

  void bfs(size_t start_node_id, std::function<void(const NodeType &)> callback) const override {
	std::vector<bool> visited(getNodesNumber(), false);
	std::queue<size_t> queue;

	queue.push(start_node_id);
	while (!queue.empty()) {
	  size_t node_id = queue.front();
	  queue.pop();

	  if (!visited[node_id]) {
		callback(nodes_[node_id]);
		visited[node_id] = true;

		for (const EdgeType &neighbor : adjacency_list_.at(node_id)) {
		  if (!visited[neighbor.target]) {
			queue.push(neighbor.target);
		  }
		}
	  }
	}
  }

 private:
  std::vector<NodeType> nodes_;
  std::vector<std::vector<EdgeType>> adjacency_list_;
};
} // namespace graph

#endif // GRAPH_SPARSE_GRAPH_HPP
