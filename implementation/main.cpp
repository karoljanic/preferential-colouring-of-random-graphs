#include <iomanip>
#include <iostream>
#include <set>

#include "avoiding_kuratowski_graphs_painter.hpp"  // graph::random::AvoidingKuratowskiGraphsPainter
#include "ba_graph.hpp"                            // graph::random::BAGraph
#include "colors_balance_globally_painter.hpp"     // graph::random::ColorsBalanceLocallyPainter
#include "colors_balance_locally_painter.hpp"      // graph::random::ColorsBalanceLocallyPainter
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
using Metric = graph::random::Metric;

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
  MetricCalculator::calculate(MetricCalculator::MetricType::MAX_FASTER, graph, metrics_map);

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

//  saveMetrics(metrics_map, std::cout);

  saveGraph(graph, metrics_map, output_graph_filename);

  return 0;
}

void debug() {
  const size_t initial_nodes = 3;
  const size_t edges_per_vertex = 2;
  const size_t final_nodes = 12;

  const int colors_numbers = 2;
  const std::vector<ColorType> colors(allColors.begin(), allColors.begin() + colors_numbers);

  RandomGraphFactory random_graph_factory;
  const BAGraph graph = random_graph_factory.createBarabasiAlbertWithPreferentialAttachmentBatageljBrandes(
      initial_nodes, final_nodes, edges_per_vertex);

  MetricsMap max_metrics;
  MetricCalculator::calculate(MetricCalculator::MetricType::ALL_PATHS, graph, max_metrics);

  MetricsMap max_faster_metrics;
  MetricCalculator::calculate(MetricCalculator::MetricType::EXPECTED, graph, max_faster_metrics);

  for (size_t node = 0; node < graph.getNodesNumber(); ++node) {
    std::cout << "Node " << node << "(" << graph.getDegree(node) << "):";
    for (const auto& neighbor : graph.getNeighbours(node)) {
      std::cout << "  " << neighbor.id << " ";
    }
    std::cout << std::endl;
  }

  std::ofstream file1{"all_metrics.txt"};
  saveMetrics(max_metrics, file1);

  std::ofstream file2{"expected_metrics.txt"};
  saveMetrics(max_faster_metrics, file2);

  float k33_all_sum = 0.0F;
  float k5_all_sum = 0.0F;
  for (const auto& outer_pair : max_metrics) {
    for (const auto& inner_pair : outer_pair.second) {
      k33_all_sum += inner_pair.second.k33;
      k5_all_sum += inner_pair.second.k5;
    }
  }

  float k33_expected_sum = 0.0F;
  float k5_expected_sum = 0.0F;
  for (const auto& outer_pair : max_faster_metrics) {
    for (const auto& inner_pair : outer_pair.second) {
      k33_expected_sum += inner_pair.second.k33;
      k5_expected_sum += inner_pair.second.k5;
    }
  }

  std::cout << "All sum: " << k33_all_sum << " " << k5_all_sum << std::endl;
  std::cout << "Expected sum: " << k33_expected_sum << " " << k5_expected_sum << std::endl;
}

void debug2() {
  constexpr size_t initial_nodes = 3;
  constexpr size_t edges_per_vertex = 2;
  constexpr size_t final_nodes = 12;
  constexpr size_t tries = 1000;

  float k33_all_sum = 0.0F;
  float k5_all_sum = 0.0F;
  float k33_expected_sum = 0.0F;
  float k5_expected_sum = 0.0F;
  for (size_t i = 0; i < tries; ++i) {
    RandomGraphFactory random_graph_factory;
    const BAGraph graph = random_graph_factory.createBarabasiAlbertWithPreferentialAttachmentBatageljBrandes(
        initial_nodes, final_nodes, edges_per_vertex);

    MetricsMap max_metrics;
    MetricCalculator::calculate(MetricCalculator::MetricType::ALL_PATHS, graph, max_metrics);
    for (const auto& outer_pair : max_metrics) {
      for (const auto& inner_pair : outer_pair.second) {
        k33_all_sum += inner_pair.second.k33;
        k5_all_sum += inner_pair.second.k5;
      }
    }

    MetricsMap max_faster_metrics;
    MetricCalculator::calculate(MetricCalculator::MetricType::EXPECTED, graph, max_faster_metrics);
    for (const auto& outer_pair : max_faster_metrics) {
      for (const auto& inner_pair : outer_pair.second) {
        k33_expected_sum += inner_pair.second.k33;
        k5_expected_sum += inner_pair.second.k5;
      }
    }
  }

  std::cout << "All sum: " << k33_all_sum / static_cast<float>(tries) << " " << k5_all_sum / static_cast<float>(tries)
            << std::endl;
  std::cout << "Expected sum: " << k33_expected_sum / static_cast<float>(tries) << " "
            << k5_expected_sum / static_cast<float>(tries) << std::endl;
}

void debug3() {
  const size_t initial_nodes = 3;
  const size_t edges_per_vertex = 2;
  const size_t final_nodes = 12;

  RandomGraphFactory random_graph_factory;
  const BAGraph graph = random_graph_factory.createBarabasiAlbertWithPreferentialAttachmentBatageljBrandes(
      initial_nodes, final_nodes, edges_per_vertex);

  MetricsMap all_metrics;
  MetricCalculator::calculate(MetricCalculator::MetricType::ALL_PATHS, graph, all_metrics);

  std::ofstream file{"v1.txt"};
  saveMetrics(all_metrics, file);

  graph.saveToFile(
      "g.json", [](std::ofstream& /*file*/, const BANode& /*node*/) {}, [](std::ofstream& /*file*/, const BAEdge& /*node*/) {});
}

void debug4() {
  BAGraph  graph;
  graph.loadFromFile("g.json");

  MetricsMap all_metrics;
  MetricCalculator::calculate(MetricCalculator::MetricType::ALL_PATHS, graph, all_metrics);

  std::ofstream file{"v1.txt"};
  saveMetrics(all_metrics, file);
}

int main(int argc, char* argv[]) {
//  debug4();
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
