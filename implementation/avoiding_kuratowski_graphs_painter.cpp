#include "avoiding_kuratowski_graphs_painter.hpp"

#include <chrono>  // std::chrono
#include <limits>  // std::numeric_limits
#include <stack>   // std::stack
#include <vector>  // std::vector

namespace graph::random {
AvoidingKuratowskiGraphsPainter::AvoidingKuratowskiGraphsPainter(std::vector<ColorType> edges_colors, MetricType metric_type)
    : GraphPainter({}, std::move(edges_colors)), metric_type_(metric_type) {

  for (const auto& edges_color : edges_colors_) {
    embeddings_[edges_color] = BAGraph{};
    metrics_map_[edges_color] = MetricsMap{};
  }
}

void AvoidingKuratowskiGraphsPainter::paintNode(BAGraph& /*graph*/, BANode& node) {
  for (const auto& edges_color : edges_colors_) {
    embeddings_[edges_color].addNode();

    for (auto& iter : metrics_map_[edges_color]) {
      iter.second[node.id] = Metric{};
    }

    metrics_map_[edges_color][node.id] = std::map<size_t, Metric>{};
    for (size_t node_v = 0; node_v < node.id; ++node_v) {
      metrics_map_[edges_color][node.id][node_v] = Metric{};
    }
  }
}

void AvoidingKuratowskiGraphsPainter::paintEdge(BAGraph& graph, BAEdge& edge) {
  auto start_coloring = std::chrono::high_resolution_clock::now();

  std::map<ColorType, double> metrics;
  std::map<ColorType, MetricsMap> metrics_map_copy;
  for (const auto& edges_color : edges_colors_) {
    metrics_map_copy[edges_color] = MetricsMap{metrics_map_[edges_color]};

    BAGraph embedding_copy{embeddings_[edges_color]};
    embedding_copy.addEdge(edge.source, edge.target);

    MetricCalculator::calculate(metric_type_, embedding_copy, metrics_map_copy[edges_color]);
  }

  auto end_metrics_calculation = std::chrono::high_resolution_clock::now();

  for (const auto& edges_color : edges_colors_) {
    const Metric metric = sumMetrics(metrics_map_copy[edges_color]);

    metrics[edges_color] = scaleFinalMetric(metric);
  }

  auto end_metrics_merge = std::chrono::high_resolution_clock::now();

  //  double lowest_metric_value{std::numeric_limits<double>::max()};
  //  for (const auto& edges_color : edges_colors_) {
  //    if (lowest_metric_value > metrics[edges_color]) {
  //      lowest_metric_value = metrics[edges_color];
  //    }
  //  }

  double lowest_metric_value{std::numeric_limits<double>::lowest()};
    for (const auto& edges_color : edges_colors_) {
      if (lowest_metric_value < metrics[edges_color]) {
        lowest_metric_value = metrics[edges_color];
      }
    }

  std::vector<ColorType> low_colors;
  for (const auto& edges_color : edges_colors_) {
    if (std::fabs(lowest_metric_value - metrics[edges_color]) < EPSILON) {
      low_colors.push_back(edges_color);
    }
  }

  std::uniform_int_distribution<size_t> distribution(0, low_colors.size() - 1);
  const size_t index = distribution(generator_);
  embeddings_[low_colors[index]].addEdge(edge.source, edge.target);
  edge.color = low_colors[index];
  embeddings_[low_colors[index]].getEdge(edge.source, edge.target).color = low_colors[index];
  metrics_map_[low_colors[index]] = metrics_map_copy[low_colors[index]];

  auto end_coloring = std::chrono::high_resolution_clock::now();

  auto total_duration = std::chrono::duration_cast<std::chrono::microseconds>(end_coloring - start_coloring);
  auto metrics_calculation_duration =
      std::chrono::duration_cast<std::chrono::microseconds>(end_metrics_calculation - start_coloring);
  auto metrics_merge_duration =
      std::chrono::duration_cast<std::chrono::microseconds>(end_metrics_merge - end_metrics_calculation);

  coloring_times_vector_.push_back({.edges_number = graph.getEdgesNumber(),
                                    .total_time = static_cast<double>(total_duration.count()),
                                    .metrics_calculation_time = static_cast<double>(metrics_calculation_duration.count()),
                                    .metrics_merge_time = static_cast<double>(metrics_merge_duration.count())});

//  std::cout << "Painting edge: " << edge.source << " " << edge.target << " to " << low_colors[index].hex << std::endl;
}

void AvoidingKuratowskiGraphsPainter::reset() {
  *this = AvoidingKuratowskiGraphsPainter(edges_colors_, metric_type_);
}

Metric AvoidingKuratowskiGraphsPainter::sumMetrics(const MetricsMap& metrics_map) {
  double k33_sum{0.0F};
  double k5_sum{0.0F};

  for (const auto& outerPair : metrics_map) {
    for (const auto& innerPair : outerPair.second) {
      k33_sum += innerPair.second.k33;
      k5_sum += innerPair.second.k5;
    }
  }

  return {.k33 = k33_sum, .k5 = k5_sum};
}

double AvoidingKuratowskiGraphsPainter::scaleFinalMetric(Metric metric) {
  return (499.0F * metric.k33 + 1.0F * metric.k5) / 500.0F;
}
}  // namespace graph::random
