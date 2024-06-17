#include <iostream>

#include "ba_graph.hpp"                    // graph::random::BAGraph
#include "random_graph_factory.hpp"        // graph::random::RandomGraphFactory
#include "colors_balance_globally_painter.hpp"    // graph::random::ColorsBalanceLocallyPainter
#include "colors_balance_locally_painter.hpp"    // graph::random::ColorsBalanceLocallyPainter
#include "avoiding_kuratowski_graphs_painter.hpp"    // graph::random::AvoidingKuratowskiGraphsPainter

using namespace graph::random;

constexpr std::string Blue{"#1f77b4"};
constexpr std::string Orange{"#ff7f0e"};
constexpr std::string Green{"#2ca02c"};
constexpr std::string Red{"#d62728"};
constexpr std::string Purple{"#9467bd"};
constexpr std::string Brown{"#8c564b"};
constexpr std::string Pink{"#e377c2"};
constexpr std::string Gray{"#7f7f7f"};
constexpr std::string Olive{"#bcbd22"};
constexpr std::string Cyan{"#17becf"};

std::vector<std::string> colors{Blue, Orange, Green};

void printDegreesHistogram(const std::map<size_t, size_t> &histogram) {
  for (const auto &[degree, count] : histogram) {
	std::cout << "Degree: " << degree << " Count: " << count << std::endl;
  }
  std::cout << std::endl;
}

void printCustomHistogram(const ColorsBalanceGloballyPainter &painter) {
  for (const auto &[color, count] : painter.getColorsHistogram()) {
	std::cout << "Color: " << color << " Count: " << count << std::endl;
  }
  std::cout << std::endl;
}

void printCustomHistogram(const ColorsBalanceLocallyPainter &painter) {
  for (const auto &[color, count] : painter.getColorsHistogram()) {
	std::cout << "Color: " << color << " Count: " << count << std::endl;
  }
  std::cout << std::endl;
}

void printCustomHistogram(const AvoidingKuratowskiGraphsPainter &painter) {

}

void saveWithColors(const graph::Graph<BANode, BAEdge> &graph,
					const std::string &filename, bool node_coloring, bool edge_coloring) {

  auto node_color = node_coloring ?
					[](std::ofstream &file, const BANode &node) {
					  file << R"(, "color": ")" << node.color << "\"";
					} : [](std::ofstream &file, const BANode &node) {};

  auto edge_color = edge_coloring ?
					[](std::ofstream &file, const BAEdge &edge) {
					  file << R"(, "color": ")" << edge.color << "\"";
					} : [](std::ofstream &file, const BAEdge &edge) {};

  graph.saveToFile(filename, node_color, edge_color);
}

int main() {
//  ColorsBalanceGloballyPainter painter{colors};
//  ColorsBalanceLocallyPainter painter{colors};
  AvoidingKuratowskiGraphsPainter painter{colors};

  const BAGraph graph1 =
	  RandomGraphFactory::createBarabasiAlbertWithPreferentialAttachment(10,
																		 50,
																		 2,
																		 1.2F,
																		 &painter);

  std::cout << "Graph 1 density: " << graph1.getDensity() << std::endl;
  std::cout << "Graph 1 degrees histogram:" << std::endl;
  printDegreesHistogram(graph1.getDegreesHistogram());
  printCustomHistogram(painter);
  graph1.saveToFile("graph1.txt");
  saveWithColors(graph1, "graph1.txt", true, true);

  painter.reset();
  const BAGraph graph2 = RandomGraphFactory::createBarabasiAlbertWithLinkSelection(10,
																				   50,
																				   2,
																				   &painter);
  std::cout << "Graph 2 density: " << graph2.getDensity() << std::endl;
  std::cout << "Graph 2 degrees histogram:" << std::endl;
  printDegreesHistogram(graph2.getDegreesHistogram());
  printCustomHistogram(painter);
  graph2.saveToFile("graph2.txt");
  saveWithColors(graph2, "graph2.txt", true, true);

  painter.reset();
  const BAGraph graph3 = RandomGraphFactory::createBarabasiAlbertWithCopyingModel(10,
																				  50,
																				  2,
																				  0.5F,
																				  &painter);

  std::cout << "Graph 3 density: " << graph3.getDensity() << std::endl;
  std::cout << "Graph 3 degrees histogram:" << std::endl;
  printDegreesHistogram(graph3.getDegreesHistogram());
  printCustomHistogram(painter);
  graph3.saveToFile("graph3.txt");
  saveWithColors(graph3, "graph3.txt", true, true);

  return 0;
}
