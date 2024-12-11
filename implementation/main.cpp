#include <iomanip>
#include <iostream>
#include <set>

#include "avoiding_kuratowski_graphs_painter.hpp"  // graph::random::AvoidingKuratowskiGraphsPainter
#include "ba_graph.hpp"                            // graph::random::BAGraph
#include "colors_balance_globally_painter.hpp"     // graph::random::ColorsBalanceLocallyPainter
#include "colors_balance_locally_painter.hpp"      // graph::random::ColorsBalanceLocallyPainter
#include "max_planar_subgraph.hpp"                 // graph::random::MaxPlanarSubgraph
#include "metric_calculator.hpp"                   // graph::random::StaticPainter
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
using MetricCalculator = graph::random::MetricCalculator;
using MetricsMap = graph::random::MetricCalculator::MetricsMap;

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

void saveGraph(const graph::Graph<BANode, BAEdge>& graph,
               const std::map<size_t, std::map<size_t, MetricCalculator::Metric>>& metrics, const std::string& filename) {

  auto edge_metric = [&](std::ofstream& file, const BAEdge& edge) {
    if (edge.source <= edge.target) {
      file << R"(, "metrics": {)";
      file << R"("k33": )" << metrics.at(edge.source).at(edge.target).k33 << ", ";
      file << R"("k5": )" << metrics.at(edge.source).at(edge.target).k5 << "}";
    }
    else {
      file << R"(, "metrics": {)";
      file << R"("k33": )" << metrics.at(edge.target).at(edge.source).k33 << ", ";
      file << R"("k5": )" << metrics.at(edge.target).at(edge.source).k5 << "}";
    }
  };

  graph.saveToFile(
      filename, [](std::ofstream& /*file*/, const BANode& /*node*/) {}, edge_metric);
}

void saveMetrics(const MetricsMap& metrics, std::ostream& out) {
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
  for (const auto& time : painter.getColoringTimes()) {
    file << time.edges_number << " " << time.total_time << " " << time.metrics_calculation_time << " " << time.metrics_merge_time
         << std::endl;
  }
}

int generateGraph(char* argv[]) {
  const size_t initial_nodes = std::stoi(argv[2]);
  const size_t edges_per_vertex = std::stoi(argv[3]);
  const size_t final_nodes = std::stoi(argv[4]);

  RandomGraphFactory random_graph_factory;
  const BAGraph graph = random_graph_factory.createBarabasiAlbertWithPreferentialAttachmentBatageljBrandes(
      initial_nodes, final_nodes, edges_per_vertex);

  saveGraph(graph, argv[5], false, false);

  return 0;
}

int generateAndColorGraph(char* argv[]) {
  const size_t initial_nodes = std::stoi(argv[2]);
  const size_t edges_per_vertex = std::stoi(argv[3]);
  const size_t final_nodes = std::stoi(argv[4]);
  const int colors_numbers = std::stoi(argv[5]);

  const std::vector<ColorType> colors(allColors.begin(), allColors.begin() + colors_numbers);
  AvoidingKuratowskiGraphsPainter painter{colors, MetricCalculator::MetricType::EXPECTED};

  RandomGraphFactory random_graph_factory;
  const BAGraph graph = random_graph_factory.createBarabasiAlbertWithPreferentialAttachmentBatageljBrandes(
      initial_nodes, final_nodes, edges_per_vertex, &painter);

  saveGraph(graph, argv[6] + std::string("/graph.json"), true, true);
  for (const auto& embedding : painter.getEmbeddings()) {
    saveGraph(embedding.second, argv[6] + std::string("/embedding_") + std::string{embedding.first.hex} + ".txt", true, true);

    std::ofstream file{argv[6] + std::string("/metrics_") + std::string{embedding.first.hex} + ".txt"};
    saveMetrics(painter.getMetricsMap().at(embedding.first), file);
  }

  saveTimes(painter, argv[6] + std::string("/times.txt"));

  return 0;
}

int generateAndCalculateMetricForGraph(char* argv[]) {
  const size_t initial_nodes = std::stoi(argv[2]);
  const size_t edges_per_vertex = std::stoi(argv[3]);
  const size_t final_nodes = std::stoi(argv[4]);

  RandomGraphFactory random_graph_factory;
  const BAGraph graph = random_graph_factory.createBarabasiAlbertWithPreferentialAttachmentBatageljBrandes(
      initial_nodes, final_nodes, edges_per_vertex);

  MetricsMap metrics_map;
  MetricCalculator::calculate(MetricCalculator::MetricType::EXPECTED, graph, metrics_map);

  saveGraph(graph, metrics_map, argv[5]);

  return 0;
}

int setMetricToGraph(char* argv[]) {
  const std::string input_graph_filename = argv[2];
  const std::string output_graph_filename = argv[3];

  BAGraph graph;
  graph.loadFromFile(input_graph_filename);

  MetricsMap metrics_map;
  MetricCalculator::calculate(MetricCalculator::MetricType::EXPECTED, graph, metrics_map);

  saveGraph(graph, metrics_map, output_graph_filename);

  return 0;
}

int main(int argc, char* argv[]) {
  if (argc < 2) {
    return 0;
  }

  constexpr int generate_graph_arguments_number = 6;
  constexpr int color_graph_arguments_number = 7;
  constexpr int calculate_metric_arguments_number = 6;
  constexpr int set_metric_to_graph_arguments_number = 4;

  const char mode = argv[1][0];

  if (mode == 'G' && argc == generate_graph_arguments_number) {
    const int res = generateGraph(argv);
    return res;
  }

  if (mode == 'C' && argc == color_graph_arguments_number) {
    const int res = generateAndColorGraph(argv);
    return res;
  }

  if (mode == 'M' && argc == calculate_metric_arguments_number) {
    const int res = generateAndCalculateMetricForGraph(argv);
    return res;
  }

  if (mode == 'S' && argc == set_metric_to_graph_arguments_number) {
    const int res = setMetricToGraph(argv);
    return res;
  }

  std::cerr << "Usage: " << argv[0] << " <mode> <initial_nodes> <edges_per_vertex> <final_nodes> <output_file>" << std::endl;
  std::cerr << "or" << std::endl;
  std::cerr << "Usage: " << argv[0] << " <mode> <initial_nodes> <edges_per_vertex> <final_nodes> <colors_number> <output_dir>"
            << std::endl;
  std::cerr << "or" << std::endl;
  std::cerr << "Usage: " << argv[0] << " <mode> <input_file> <output_file>" << std::endl;

  return 0;
}
