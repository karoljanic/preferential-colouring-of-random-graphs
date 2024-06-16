#include <iostream>

#include "sparse_graph.hpp"
#include "random_graph_factory.hpp"

void printHistogram(const std::map<size_t, size_t> &histogram) {
  for (const auto &[degree, count] : histogram) {
	std::cout << "Degree: " << degree << " Count: " << count << std::endl;
  }
  std::cout << std::endl;
}

void saveWithColors(const graph::Graph<graph::random::BANode, graph::random::BAEdge> &graph,
					const std::string &filename, bool node_coloring, bool edge_coloring) {

  auto node_color = node_coloring ?
					[](std::ofstream &file, const graph::random::BANode &node) {
					  file << R"(, "color": ")" << node.color << "\"";
					} : [](std::ofstream &file, const graph::random::BANode &node) {};

  auto edge_color = edge_coloring ?
					[](std::ofstream &file, const graph::random::BAEdge &edge) {
					  file << R"(, "color": ")" << edge.color << "\"";
					} : [](std::ofstream &file, const graph::random::BAEdge &edge) {};

  graph.saveToFile(filename, node_color, edge_color);
}

int main() {
  const graph::random::BAGraph
	  graph1 = graph::random::RandomGraphFactory::createBarabasiAlbertWithPreferentialAttachment(10, 50, 2, 1.2F);
  std::cout << "Graph 1 density: " << graph1.getDensity() << std::endl;
  std::cout << "Graph 1 degrees histogram:" << std::endl;
  printHistogram(graph1.getDegreesHistogram());
  graph1.saveToFile("graph1.txt");
  saveWithColors(graph1, "graph1-nc.txt", true, false);

  const graph::random::BAGraph
	  graph2 = graph::random::RandomGraphFactory::createBarabasiAlbertWithLinkSelection(10, 50, 2);
  std::cout << "Graph 2 density: " << graph2.getDensity() << std::endl;
  std::cout << "Graph 2 degrees histogram:" << std::endl;
  printHistogram(graph2.getDegreesHistogram());
  graph2.saveToFile("graph2.txt");
  saveWithColors(graph2, "graph2-ec.txt", false, true);

  const graph::random::BAGraph
	  graph3 = graph::random::RandomGraphFactory::createBarabasiAlbertWithCopyingModel(10, 50, 2, 0.5F);

  std::cout << "Graph 3 density: " << graph3.getDensity() << std::endl;
  std::cout << "Graph 3 degrees histogram:" << std::endl;
  printHistogram(graph3.getDegreesHistogram());
  graph3.saveToFile("graph3.txt");
  saveWithColors(graph3, "graph3-nc-ec.txt", true, true);

  return 0;
}
