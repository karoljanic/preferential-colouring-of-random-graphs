#include <iomanip>
#include <iostream>
#include <set>

#include "avoiding_kuratowski_graphs_painter.hpp"  // graph::random::AvoidingKuratowskiGraphsPainter
#include "ba_graph.hpp"                            // graph::random::BAGraph
#include "colors_balance_globally_painter.hpp"     // graph::random::ColorsBalanceLocallyPainter
#include "colors_balance_locally_painter.hpp"      // graph::random::ColorsBalanceLocallyPainter
#include "random_graph_factory.hpp"                // graph::random::RandomGraphFactory

using ColorType = graph::random::ColorType;
using BAGraph = graph::random::BAGraph;
using BANode = graph::random::BANode;
using BAEdge = graph::random::BAEdge;
using RandomGraphFactory = graph::random::RandomGraphFactory;
using GraphPainter = graph::random::GraphPainter;
using AvoidingKuratowskiGraphsPainter = graph::random::AvoidingKuratowskiGraphsPainter;
using ColorsBalanceGloballyPainter = graph::random::ColorsBalanceGloballyPainter;
using ColorsBalanceLocallyPainter = graph::random::ColorsBalanceLocallyPainter;

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

const std::vector<ColorType> allColors{Blue, Orange, Green, Red, Purple, Brown, Pink, Gray, Olive, Cyan};

void saveGraph(const graph::Graph<BANode, BAEdge>& graph, const std::string& filename, bool node_coloring, bool edge_coloring) {

  auto node_color = node_coloring ?
                                  [](std::ofstream &file, const BANode &node) {
                                    file << R"(, "color": ")" << node.color.hex << "\"";
                                  } : [](std::ofstream &/*file*/, const BANode &/*node*/) {};

  auto edge_color = edge_coloring ?
                                  [](std::ofstream &file, const BAEdge &edge) {
                                    file << R"(, "color": ")" << edge.color.hex << "\"";
                                  } : [](std::ofstream &/*file*/, const BAEdge &/*edge*/) {};

  if (node_coloring || edge_coloring) {
    graph.saveToFile(filename, node_color, edge_color);
  }
  else {
    graph.saveToFile(filename);
  }
}

void saveMetrics(const AvoidingKuratowskiGraphsPainter::MetricsMap& metrics, std::ostream& out) {
  constexpr int column_width = 15;
  constexpr int number_width = 5;

  std::set<size_t> column_keys;
  for (const auto& outer_pair : metrics) {
    for (const auto& inner_pair : outer_pair.second) {
      column_keys.insert(inner_pair.first);
    }
  }

  for (const size_t col : column_keys) {
    out << std::setw(column_width) << col;
    out << " ";
  }
  out << std::endl;

  for (const auto& outer_pair : metrics) {
    const size_t outer_key = outer_pair.first;
    out << std::setw(number_width) << outer_key;
    for (const size_t col : column_keys) {
      const auto iter = outer_pair.second.find(col);
      if (iter != outer_pair.second.end()) {
        out << std::setw(number_width) << "(" << std::fixed << std::setprecision(2) << iter->second.k33 << ", " << iter->second.k5
            << ")";
      }
      else {
        out << std::setw(number_width) << "(" << std::fixed << std::setprecision(2) << 0.0F << ", " << 0.0F << ")";
      }
    }
    out << std::endl;
  }
}

void saveTimes(const AvoidingKuratowskiGraphsPainter& painter, const std::string& file_name) {
  std::ofstream file{file_name};
  for (const auto& time_pair : painter.getColoringTimes()) {
    file << time_pair.first << " " << time_pair.second << std::endl;
  }
}

int generateGraph(char* argv[]) {
  const size_t initial_nodes = std::stoi(argv[1]);
  const size_t edges_per_vertex = std::stoi(argv[2]);
  const size_t final_nodes = std::stoi(argv[3]);

  graph::random::RandomGraphFactory random_graph_factory;
  const graph::random::BAGraph graph = random_graph_factory.createBarabasiAlbertWithPreferentialAttachmentBatageljBrandes(
      initial_nodes, final_nodes, edges_per_vertex);

  saveGraph(graph, argv[4], false, false);

  return 0;
}

int generateAndColorGraph(char* argv[]) {
  const size_t initial_nodes = std::stoi(argv[1]);
  const size_t edges_per_vertex = std::stoi(argv[2]);
  const size_t final_nodes = std::stoi(argv[3]);
  const int colors_numbers = std::stoi(argv[4]);

  const std::vector<graph::random::ColorType> colors(allColors.begin(), allColors.begin() + colors_numbers);
  graph::random::AvoidingKuratowskiGraphsPainter painter{colors};

  graph::random::RandomGraphFactory random_graph_factory;
  const graph::random::BAGraph graph = random_graph_factory.createBarabasiAlbertWithPreferentialAttachmentBatageljBrandes(
      initial_nodes, final_nodes, edges_per_vertex, &painter);

  saveGraph(graph, argv[5] + std::string("/graph.json"), true, true);
  for (const auto& embedding : painter.getEmbeddings()) {
    saveGraph(embedding.second, argv[5] + std::string("/embedding_") + std::string{embedding.first.hex} + ".txt", true, true);

    std::ofstream file{argv[5] + std::string("/metrics_") + std::string{embedding.first.hex} + ".txt"};
    saveMetrics(painter.getMetricsMap().at(embedding.first), file);
  }

  saveTimes(painter, argv[5] + std::string("/times.txt"));

  return 0;
}

int main(int argc, char* argv[]) {
  constexpr int generate_graph_arguments_number = 5;
  constexpr int color_graph_arguments_number = 6;

  if (argc == generate_graph_arguments_number) {
    const int res = generateGraph(argv);
    return res;
  }

  if (argc == color_graph_arguments_number) {
    const int res = generateAndColorGraph(argv);
    return res;
  }

  std::cerr << "Usage: " << argv[0] << " <initial_nodes> <edges_per_vertex> <final_nodes> <output_file>" << std::endl;
  std::cerr << "or" << std::endl;
  std::cerr << "Usage: " << argv[0] << " <initial_nodes> <edges_per_vertex> <final_nodes> <colors_number> <output_dir>"
            << std::endl;

  return 0;
}
