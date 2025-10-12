#ifndef GRAPH_RANDOM_PREFERENTIAL_COLORING_HPP
#define GRAPH_RANDOM_PREFERENTIAL_COLORING_HPP

#include <map>
#include <vector>

#include "ba_graph.hpp"      // graph::random::BAGraph, graph::random::BANode, graph::random::BAEdge
#include "disjoint_set.hpp"  // DisjointSet

namespace graph::random {
class PreferentialColoring {
 public:
  PreferentialColoring() = default;

  PreferentialColoring(const PreferentialColoring&) = default;
  PreferentialColoring(PreferentialColoring&&) = default;

  PreferentialColoring& operator=(const PreferentialColoring&) = default;
  PreferentialColoring& operator=(PreferentialColoring&&) = default;

  ~PreferentialColoring() = default;

  static std::vector<BAGraph> color(const BAGraph& graph, size_t colors_count);
};
}  // namespace graph::random

#endif  // GRAPH_RANDOM_PREFERENTIAL_COLORING_HPP