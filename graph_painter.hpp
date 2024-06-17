#ifndef GRAPH_RANDOM_GRAPH_PAINTER_HPP
#define GRAPH_RANDOM_GRAPH_PAINTER_HPP

#include <vector>    // std::vector

#include "ba_graph.hpp"    // graph::random::BAGraph, graph::random::BANode, graph::random::BAEdge

namespace graph::random {
class GraphPainter {
 public:
  GraphPainter() = default;
  GraphPainter(std::vector<std::string> nodes_colors, std::vector<std::string> edges_colors)
	  : nodes_colors_(std::move(nodes_colors)), edges_colors_(std::move(edges_colors)) {}

  GraphPainter(const GraphPainter &) = default;
  GraphPainter(GraphPainter &&) = default;

  GraphPainter &operator=(const GraphPainter &) = default;
  GraphPainter &operator=(GraphPainter &&) = default;

  virtual ~GraphPainter() = default;

  [[nodiscard]]
  inline size_t getNodesColorsNumber() const { return nodes_colors_.size(); }
  [[nodiscard]]
  inline size_t getEdgesColorsNumber() const { return edges_colors_.size(); }

  virtual void paintNode(BAGraph &graph, BANode &node) = 0;
  virtual void paintEdge(BAGraph &graph, BAEdge &edge) = 0;

  virtual void reset() = 0;

 protected:
  std::vector<std::string> nodes_colors_;
  std::vector<std::string> edges_colors_;
};
} // namespace graph::random

#endif // GRAPH_RANDOM_GRAPH_PAINTER_HPP
