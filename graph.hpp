#ifndef GRAPH_GRAPH_HPP
#define GRAPH_GRAPH_HPP

#include <concepts>      // concepts
#include <cstddef>       // std::size_t
#include <functional>    // std::function
#include <map>           // std::map
#include <vector>        // std::vector

namespace graph {
template<typename T>
concept HasId = requires(T t) {
  { t.id } -> std::convertible_to<std::size_t>;
};

template<typename T>
concept HasSourceAndTarget = requires(T t) {
  { t.source } -> std::convertible_to<std::size_t>;
  { t.target } -> std::convertible_to<std::size_t>;
};

template<typename NodeType, typename EdgeType> requires HasId<NodeType> && HasSourceAndTarget<EdgeType>
class Graph {
 public:
  Graph() = default;

  Graph(const Graph &) = default;
  Graph(Graph &&) noexcept = default;

  Graph &operator=(const Graph &) = default;
  Graph &operator=(Graph &&) noexcept = default;

  virtual ~Graph() = default;

  virtual void addNode(NodeType &node) = 0;
  virtual bool addEdge(EdgeType &edge) = 0;

  [[nodiscard]]
  virtual NodeType getNode(size_t node_id) const = 0;
  [[nodiscard]]
  virtual EdgeType getEdge(size_t source, size_t target) const = 0;

  [[nodiscard]]
  virtual size_t getNodesNumber() const = 0;
  [[nodiscard]]
  virtual size_t getEdgesNumber() const = 0;
  [[nodiscard]]
  virtual float getDensity() const = 0;

  [[nodiscard]]
  virtual std::vector<NodeType> getNeighbours(size_t node_id) const = 0;
  [[nodiscard]]
  virtual std::vector<EdgeType> getAdjacentEdges(size_t node_id) const = 0;

  [[nodiscard]]
  virtual size_t getDegree(size_t node_id) const = 0;
  [[nodiscard]]
  virtual std::map<size_t, size_t> getDegreesHistogram() const = 0;

  virtual void dfs(size_t start_node_id, std::function<void(const NodeType &)> callback) const = 0;
  virtual void bfs(size_t start_node_id, std::function<void(const NodeType &)> callback) const = 0;
};
} // namespace graph

#endif // GRAPH_GRAPH_HPP
