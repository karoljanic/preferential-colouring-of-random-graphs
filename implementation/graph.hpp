#ifndef GRAPH_GRAPH_HPP
#define GRAPH_GRAPH_HPP

#include <concepts>    // concept, requires, std::convertible_to
#include <cstddef>     // std::size_t
#include <functional>  // std::function
#include <map>         // std::map
#include <string>      // std::string
#include <vector>      // std::vector

#include <nlohmann/json.hpp>
using json = nlohmann::json;

namespace graph {
template <typename T>
concept HasId = requires(T t) {
  { t.id } -> std::convertible_to<std::size_t>;
};

template <typename T>
concept HasSourceAndTarget = requires(T t) {
  { t.source } -> std::convertible_to<std::size_t>;
  { t.target } -> std::convertible_to<std::size_t>;
};

template <typename NodeType, typename EdgeType>
requires HasId<NodeType>&& HasSourceAndTarget<EdgeType> class Graph {
 public:
  Graph() = default;

  Graph(const Graph&) = default;
  Graph(Graph&&) noexcept = default;

  Graph& operator=(const Graph&) = default;
  Graph& operator=(Graph&&) noexcept = default;

  virtual ~Graph() = default;

  virtual size_t addNode() = 0;
  virtual void addEdge(size_t node1_id, size_t node2_id) = 0;
  [[nodiscard]] virtual bool edgeExists(size_t node1_id, size_t node2_id) const = 0;

  [[nodiscard]] virtual NodeType& getNode(size_t node_id) = 0;
  [[nodiscard]] virtual NodeType& getLastAddedNode() = 0;
  [[nodiscard]] virtual std::vector<NodeType> getNodes() const = 0;

  [[nodiscard]] virtual EdgeType& getEdge(size_t source, size_t target) = 0;
  [[nodiscard]] virtual EdgeType& getLastAddedEdge() = 0;
  [[nodiscard]] virtual std::vector<EdgeType> getEdges() const = 0;

  [[nodiscard]] virtual size_t getNodesNumber() const = 0;
  [[nodiscard]] virtual size_t getEdgesNumber() const = 0;
  [[nodiscard]] virtual float getDensity() const = 0;

  [[nodiscard]] virtual std::vector<NodeType> getNeighbours(size_t node_id) const = 0;
  [[nodiscard]] virtual std::vector<EdgeType> getAdjacentEdges(size_t node_id) const = 0;

  [[nodiscard]] virtual size_t getDegree(size_t node_id) const = 0;
  [[nodiscard]] virtual std::map<size_t, size_t> getDegreesHistogram() const = 0;

  virtual void dfs(size_t start_node_id, std::function<void(const NodeType&)> callback) const = 0;
  virtual void bfs(size_t start_node_id, std::function<void(const NodeType&)> callback) const = 0;

  virtual void saveToFile(const std::string& filename) const = 0;
  virtual void saveToFile(const std::string& filename, std::function<void(std::ofstream& file, const NodeType&)> vertexCallback,
                          std::function<void(std::ofstream& file, const EdgeType&)> edgeCallback) const = 0;

  virtual int loadFromFile(const std::string& filename) = 0;
};
}  // namespace graph

#endif  // GRAPH_GRAPH_HPP
