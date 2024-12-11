#ifndef GRAPH_RANDOM_MAX_PLANAR_SUBGRAPH_HPP
#define GRAPH_RANDOM_MAX_PLANAR_SUBGRAPH_HPP

#include <functional>
#include <queue>
#include <stack>
#include <utility>
#include <vector>

#include "ba_graph.hpp"        // graph::random::BAGraph, graph::random::BANode, graph::random::BAEdge
#include "disjoint_set.hpp"    // DisjointSet
#include "planarity_test.hpp"  // PlanarityTest

namespace graph::random {
class MaxPlanarSubgraph {
 public:
  MaxPlanarSubgraph() = default;

  MaxPlanarSubgraph(const MaxPlanarSubgraph&) = default;
  MaxPlanarSubgraph(MaxPlanarSubgraph&&) = default;

  MaxPlanarSubgraph& operator=(const MaxPlanarSubgraph&) = default;
  MaxPlanarSubgraph& operator=(MaxPlanarSubgraph&&) = default;

  ~MaxPlanarSubgraph() = default;

  static void mstBased(const BAGraph& graph, BAGraph& max_planar_subgraph);
  static void weightedMstBased(const BAGraph& graph, BAGraph& max_planar_subgraph);

  static void cactusBased(const BAGraph& graph, BAGraph& max_planar_subgraph);
  static void weightedCactusBased(const BAGraph& graph, BAGraph& max_planar_subgraph);

  static void maximizeSubgraph(const BAGraph& graph, BAGraph& max_planar_subgraph);
  static void weightedMaximizeSubgraph(const BAGraph& graph, BAGraph& max_planar_subgraph);

  static size_t crossingEdges(const BAGraph& graph);

 private:
  static size_t edgeWeight(const BAGraph& graph, size_t source, size_t target);

  // copy vertices from the original graph to the max planar subgraph
  static void copyVertices(const BAGraph& graph, BAGraph& max_planar_subgraph);

  // add edge to the max planar subgraph and mark it as used
  static void addEdge(size_t source, size_t target, BAGraph& max_planar_subgraph,
                      std::map<size_t, std::map<size_t, bool>>& usedEdges);

  // merge components with oldIds to the newId
  static void mergeComponents(DisjointSet& components, std::vector<size_t> ids);

  struct Compare {
    explicit Compare(const BAGraph& graph) : graph_{graph} {}

    bool operator()(const std::pair<size_t, size_t>& lhs, const std::pair<size_t, size_t>& rhs) const {
      return edgeWeight(graph_, lhs.first, lhs.second) < edgeWeight(graph_, rhs.first, rhs.second);
    }

   private:
    const BAGraph& graph_;
  };
};
}  // namespace graph::random

#endif  // GRAPH_RANDOM_MAX_PLANAR_SUBGRAPH_HPP