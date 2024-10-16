#ifndef GRAPH_RANDOM_AVOIDING_KURATOWSKI_GRAPHS_PAINTER_HPP
#define GRAPH_RANDOM_AVOIDING_KURATOWSKI_GRAPHS_PAINTER_HPP

#include <functional>  // std::function
#include <iostream>    // std::cout
#include <random>      // std::mt19937, std::random_device, std::uniform_real_distribution, std::uniform_int_distribution
#include <utility>     // std::pair

#include "graph_painter.hpp"      // graph::random::GraphPainter
#include "metric_calculator.hpp"  // graph::random::MetricCalculator

namespace graph::random {
using Metric = MetricCalculator::Metric;
using MetricType = MetricCalculator::MetricType;
using MetricsMap = graph::random::MetricCalculator::MetricsMap;

class AvoidingKuratowskiGraphsPainter : public GraphPainter {
 public:
  struct TimeInfo {
    size_t edges_number;
    double total_time;
    double metrics_calculation_time;
    double metrics_merge_time;
  };

  AvoidingKuratowskiGraphsPainter() = default;
  explicit AvoidingKuratowskiGraphsPainter(std::vector<ColorType> edges_colors, MetricType metric_type);

  AvoidingKuratowskiGraphsPainter(const AvoidingKuratowskiGraphsPainter&) = default;
  AvoidingKuratowskiGraphsPainter(AvoidingKuratowskiGraphsPainter&&) = default;

  AvoidingKuratowskiGraphsPainter& operator=(const AvoidingKuratowskiGraphsPainter&) = default;
  AvoidingKuratowskiGraphsPainter& operator=(AvoidingKuratowskiGraphsPainter&&) = default;

  ~AvoidingKuratowskiGraphsPainter() override = default;

  void paintNode(BAGraph& graph, BANode& node) override;
  void paintEdge(BAGraph& graph, BAEdge& edge) override;

  [[nodiscard]] const std::map<ColorType, BAGraph>& getEmbeddings() const { return embeddings_; }
  [[nodiscard]] const std::map<ColorType, MetricsMap>& getMetricsMap() const { return metrics_map_; }
  [[nodiscard]] const std::vector<TimeInfo>& getColoringTimes() const { return coloring_times_vector_; }

  void reset() override;

 private:
  [[nodiscard]] static MetricCalculator::Metric sumMetrics(const MetricsMap& metrics_map);
  [[nodiscard]] static double scaleFinalMetric(Metric metric);

  std::mt19937 generator_{std::random_device{}()};
  MetricType metric_type_;
  std::map<ColorType, BAGraph> embeddings_;
  std::map<ColorType, MetricsMap> metrics_map_;
  std::vector<TimeInfo> coloring_times_vector_;

  static constexpr double EPSILON{1e-6};
};
}  // namespace graph::random

#endif  // GRAPH_RANDOM_AVOIDING_KURATOWSKI_GRAPHS_PAINTER_HPP
