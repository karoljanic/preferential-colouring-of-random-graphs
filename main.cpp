#include <iomanip>
#include <iostream>
#include <set>

#include "avoiding_kuratowski_graphs_painter.hpp"  // graph::random::AvoidingKuratowskiGraphsPainter
#include "ba_graph.hpp"                            // graph::random::BAGraph
#include "colors_balance_globally_painter.hpp"     // graph::random::ColorsBalanceLocallyPainter
#include "colors_balance_locally_painter.hpp"      // graph::random::ColorsBalanceLocallyPainter
#include "random_graph_factory.hpp"                // graph::random::RandomGraphFactory

using namespace graph::random;

constexpr ColorType Blue{"#1f77b4"};
constexpr ColorType Orange{"#ff7f0e"};
constexpr ColorType Green{"#2ca02c"};
constexpr ColorType Red{"#d62728"};
constexpr ColorType Purple{"#9467bd"};
constexpr ColorType Brown{"#8c564b"};
constexpr ColorType Pink{"#e377c2"};
constexpr ColorType Gray{"#7f7f7f"};
constexpr ColorType Olive{"#bcbd22"};
constexpr ColorType Cyan{"#17becf"};

std::vector<ColorType> colors{Blue, Orange, Green};

void printDegreesHistogram(const std::map<size_t, size_t>& histogram) {
  for (const auto& [degree, count] : histogram) {
    std::cout << "Degree: " << degree << " Count: " << count << std::endl;
  }
  std::cout << std::endl;
}

void saveWithColors(const graph::Graph<BANode, BAEdge>& graph, const std::string& filename, bool node_coloring,
                    bool edge_coloring) {

  auto node_color = node_coloring ?
					[](std::ofstream &file, const BANode &node) {
					  file << R"(, "color": ")" << node.color.hex << "\"";
					} : [](std::ofstream &/*file*/, const BANode &/*node*/) {};

  auto edge_color = edge_coloring ?
					[](std::ofstream &file, const BAEdge &edge) {
					  file << R"(, "color": ")" << edge.color.hex << "\"";
					} : [](std::ofstream &/*file*/, const BAEdge &/*edge*/) {};

  graph.saveToFile(filename, node_color, edge_color);
}

void printMetrics(const AvoidingKuratowskiGraphsPainter::MetricsMap& metrics) {
  std::set<size_t> columnKeys;
  for (const auto& outerPair : metrics) {
    for (const auto& innerPair : outerPair.second) {
      columnKeys.insert(innerPair.first);
    }
  }

  for (size_t col : columnKeys) {
    std::cout << std::setw(15) << col;
    std::cout << " ";
  }
  std::cout << std::endl;

  for (const auto& outerPair : metrics) {
    size_t outerKey = outerPair.first;
    std::cout << std::setw(5) << outerKey;
    for (size_t col : columnKeys) {
      auto it = outerPair.second.find(col);
      if (it != outerPair.second.end()) {
        std::cout << std::setw(5) << "(" << std::fixed << std::setprecision(2) << it->second.k33 << ", " << it->second.k5 << ")";
      }
      else {
        std::cout << std::setw(5) << "(" << std::fixed << std::setprecision(2) << 0.0F << ", " << 0.0F << ")";
      }
    }
    std::cout << std::endl;
  }
}

void saveTimes(const AvoidingKuratowskiGraphsPainter& painter, const std::string& file_name) {
  std::ofstream file{file_name};
  for (const auto& time_pair : painter.coloring_times_vector_) {
    file << time_pair.first << " " << time_pair.second << std::endl;
  }
}

int main() {
  //  ColorsBalanceGloballyPainter painter{colors};
  //  ColorsBalanceLocallyPainter painter{colors};
  AvoidingKuratowskiGraphsPainter painter{colors};

  RandomGraphFactory random_graph_factory;

  const BAGraph graph1 = random_graph_factory.createBarabasiAlbertWithPreferentialAttachmentBatageljBrandes(10, 40, 5, &painter);
  std::cout << "Graph 1 density: " << graph1.getDensity() << std::endl;
  std::cout << "Graph 1 degrees histogram:" << std::endl;
  printDegreesHistogram(graph1.getDegreesHistogram());
  saveWithColors(graph1, "graph1.txt", true, true);

  saveTimes(painter, "times1.txt");

  for (const auto& embedding : painter.getEmbeddings()) {
    std::cout << "Embedding density: " << embedding.second.getDensity() << std::endl;
    std::cout << "Embedding degrees histogram:" << std::endl;
    printDegreesHistogram(embedding.second.getDegreesHistogram());
    saveWithColors(embedding.second, "embedding_" + std::string{embedding.first.hex} + ".txt", true, true);
  }

  std::cout << "Total edges: " << graph1.getEdgesNumber() << std::endl;
  for (const auto& embedding : painter.getEmbeddings()) {
    std::cout << "Embedding " << embedding.first.hex << " edges: " << embedding.second.getEdgesNumber() << std::endl;
  }

  //  painter.reset();
  //
  //  const BAGraph graph2 = random_graph_factory.createBarabasiAlbertWithLinkSelection(5,
  //																				   2000,
  //																				   3,
  //																				   &painter);
  //  std::cout << "Graph 2 density: " << graph2.getDensity() << std::endl;
  //  std::cout << "Graph 2 degrees histogram:" << std::endl;
  //  printDegreesHistogram(graph2.getDegreesHistogram());
  //  saveWithColors(graph2, "graph2.txt", true, true);
  //
  //  painter.reset();
  //  const BAGraph graph3 = random_graph_factory.createBarabasiAlbertWithCopyingModel(5,
  //																				  2000,
  //																				  3,
  //																				  0.5F,
  //																				  &painter);
  //  std::cout << "Graph 3 density: " << graph3.getDensity() << std::endl;
  //  std::cout << "Graph 3 degrees histogram:" << std::endl;
  //  printDegreesHistogram(graph3.getDegreesHistogram());
  //  saveWithColors(graph3, "graph3.txt", true, true);

  return 0;
}
