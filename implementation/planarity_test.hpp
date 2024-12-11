#ifndef PLANARITY_TEST_HPP
#define PLANARITY_TEST_HPP

#include <vector>

#include <ogdf/basic/Graph.h>
#include <ogdf/planarity/MaximalPlanarSubgraphSimple.h>
#include <ogdf/planarity/PlanarSubgraphBoyerMyrvold.h>

#include "ba_graph.hpp"

class PlanarityTest {
 public:
  PlanarityTest() = default;

  PlanarityTest(const PlanarityTest&) = default;
  PlanarityTest(PlanarityTest&&) = default;

  PlanarityTest& operator=(const PlanarityTest&) = default;
  PlanarityTest& operator=(PlanarityTest&&) = default;

  ~PlanarityTest() = default;

  static bool isPlanar(const graph::random::BAGraph& graph) {
    ogdf::Graph ogdf_graph;
    toOgdf(graph, ogdf_graph);

    return ogdf::isPlanar(ogdf_graph);
  }

  static void boyerMyrvoldPlanarSubgraph(const graph::random::BAGraph& graph, graph::random::BAGraph& subgraph) {
    ogdf::Graph ogdf_graph;
    toOgdf(graph, ogdf_graph);

    ogdf::PlanarSubgraphBoyerMyrvold psbm;
    ogdf::List<ogdf::edge> delEdges;
    psbm.call(ogdf_graph, delEdges);

    for (auto* edge : delEdges) {
      ogdf_graph.delEdge(edge);
    }

    subgraph = graph::random::BAGraph();
    fromOgdf(ogdf_graph, subgraph);
  }

 private:
  static void toOgdf(const graph::random::BAGraph& graph, ogdf::Graph& ogdf_graph) {
    std::vector<ogdf::node> nodes(graph.getNodesNumber());
    for (size_t i = 0; i < graph.getNodesNumber(); ++i) {
      nodes[i] = ogdf_graph.newNode();
    }

    for (const auto& edge : graph.getEdges()) {
      ogdf_graph.newEdge(nodes[edge.source], nodes[edge.target]);
    }
  }

  static void fromOgdf(const ogdf::Graph& ogdf_graph, graph::random::BAGraph& graph) {
    for (const auto& node : ogdf_graph.nodes) {
      graph.addNode();
    }

    for (const auto& edge : ogdf_graph.edges) {
      graph.addEdge(edge->source()->index(), edge->target()->index());
    }
  }
};

#endif  // PLANARITY_TEST_HPP