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

  std::map<ColorType, float> metrics;
  std::map<ColorType, MetricsMap> metrics_map_copy;
  for (const auto& edges_color : edges_colors_) {
    metrics_map_copy[edges_color] = MetricsMap{metrics_map_[edges_color]};

    BAGraph embedding_copy{embeddings_[edges_color]};
    embedding_copy.addEdge(edge.source, edge.target);

    MetricCalculator::calculate(metric_type_, embedding_copy, metrics_map_copy[edges_color]);
  }

  auto end_metrics_calculation = std::chrono::high_resolution_clock::now();

  std::cout << std::endl;
  for (const auto& edges_color : edges_colors_) {
    const Metric metric = sumMetrics(metrics_map_copy[edges_color]);
//    const Metric metric = productSubgraphMetrics(metrics_map_copy[edges_color]);

    metrics[edges_color] = metric.k33;
//    metrics[edges_color] = metric.k33 + metric.k5;
    std::cout << "Color: " << edges_color.hex << " metric: " << metrics[edges_color] << std::endl;
  }

  auto end_metrics_merge = std::chrono::high_resolution_clock::now();

//  float lowest_metric_value{std::numeric_limits<float>::max()};
//  for (const auto& edges_color : edges_colors_) {
//    if (lowest_metric_value > metrics[edges_color]) {
//      lowest_metric_value = metrics[edges_color];
//    }
//  }

  float lowest_metric_value{std::numeric_limits<float>::lowest()};
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
                                    .total_time = static_cast<float>(total_duration.count()),
                                    .metrics_calculation_time = static_cast<float>(metrics_calculation_duration.count()),
                                    .metrics_merge_time = static_cast<float>(metrics_merge_duration.count())});

  std::cout << "Painting edge: " << edge.source << " " << edge.target << " to " << low_colors[index].hex << std::endl;
}

void AvoidingKuratowskiGraphsPainter::reset() {
  *this = AvoidingKuratowskiGraphsPainter(edges_colors_, metric_type_);
}

Metric AvoidingKuratowskiGraphsPainter::sumMetrics(const MetricsMap& metrics_map) {
  float k33_sum{0.0F};
  float k5_sum{0.0F};

  for (const auto& outerPair : metrics_map) {
    for (const auto& innerPair : outerPair.second) {
      k33_sum += innerPair.second.k33;
      k5_sum += innerPair.second.k5;
    }
  }

  return {.k33 = k33_sum, .k5 = k5_sum};
}

Metric AvoidingKuratowskiGraphsPainter::productSubgraphMetrics(const MetricsMap& metrics_map) {
  constexpr float THRESHOLD{0.1F};

  float k33_sum{0.0F};
  float k5_sum{0.0F};

  const size_t nodes_number{metrics_map.size()};

  // check for K3,3
  for (size_t node1 = 0; node1 < nodes_number; ++node1) {
    for (size_t node2 = node1 + 1; node2 < nodes_number; ++node2) {
      for (size_t node3 = node2 + 1; node3 < nodes_number; ++node3) {
        for (size_t node4 = 0; node4 < nodes_number; ++node4) {
          if (node1 == node4 || node2 == node4 || node3 == node4) {
            continue;
          }
          if (metrics_map.at(node1).at(node4).k33 < THRESHOLD || metrics_map.at(node2).at(node4).k33 < THRESHOLD ||
              metrics_map.at(node3).at(node4).k33 < THRESHOLD) {
            continue;
          }

          for (size_t node5 = node4 + 1; node5 < nodes_number; ++node5) {
            if (node1 == node5 || node2 == node5 || node3 == node5) {
              continue;
            }
            if (metrics_map.at(node1).at(node5).k33 < THRESHOLD || metrics_map.at(node2).at(node5).k33 < THRESHOLD ||
                metrics_map.at(node3).at(node5).k33 < THRESHOLD) {
              continue;
            }

            for (size_t node6 = node5 + 1; node6 < nodes_number; ++node6) {
              if (node1 == node6 || node2 == node6 || node3 == node6) {
                continue;
              }
              if (metrics_map.at(node1).at(node6).k33 < THRESHOLD || metrics_map.at(node2).at(node6).k33 < THRESHOLD ||
                  metrics_map.at(node3).at(node6).k33 < THRESHOLD) {
                continue;
              }

              k33_sum += 0.5F * metrics_map.at(node1).at(node4).k33 * metrics_map.at(node1).at((node5)).k33 *
                         metrics_map.at(node1).at(node6).k33 * metrics_map.at(node2).at(node4).k33 *
                         metrics_map.at(node2).at((node5)).k33 * metrics_map.at(node2).at(node6).k33 *
                         metrics_map.at(node3).at(node4).k33 * metrics_map.at(node3).at((node5)).k33 *
                         metrics_map.at(node3).at(node6).k33;
            }
          }
        }
      }
    }
  }

  // check for K5
  for (size_t node1 = 0; node1 < nodes_number; ++node1) {
    for (size_t node2 = node1 + 1; node2 < nodes_number; ++node2) {
      if (metrics_map.at(node1).at(node2).k5 < THRESHOLD) {
        continue;
      }

      for (size_t node3 = node2 + 1; node3 < nodes_number; ++node3) {
        if (metrics_map.at(node1).at(node3).k5 < THRESHOLD || metrics_map.at(node2).at(node3).k5 < THRESHOLD) {
          continue;
        }

        for (size_t node4 = node3 + 1; node4 < nodes_number; ++node4) {
          if (metrics_map.at(node1).at(node4).k5 < THRESHOLD || metrics_map.at(node2).at(node4).k5 < THRESHOLD ||
              metrics_map.at(node3).at(node4).k5 < THRESHOLD) {
            continue;
          }

          for (size_t node5 = node4 + 1; node5 < nodes_number; ++node5) {
            if (metrics_map.at(node1).at(node5).k5 < THRESHOLD || metrics_map.at(node2).at(node5).k5 < THRESHOLD ||
                metrics_map.at(node3).at(node5).k5 < THRESHOLD || metrics_map.at(node4).at(node5).k5 < THRESHOLD) {
              continue;
            }

            k5_sum += 0.5F * metrics_map.at(node1).at(node2).k5 * metrics_map.at(node1).at(node3).k5 *
                      metrics_map.at(node1).at(node4).k5 * metrics_map.at(node1).at(node5).k5 *
                      metrics_map.at(node2).at(node3).k5 * metrics_map.at(node2).at(node4).k5 *
                      metrics_map.at(node2).at(node5).k5 * metrics_map.at(node3).at(node4).k5 *
                      metrics_map.at(node3).at(node5).k5 * metrics_map.at(node4).at(node5).k5;
          }
        }
      }
    }
  }

  return {.k33 = k33_sum, .k5 = k5_sum};
}
}  // namespace graph::random
