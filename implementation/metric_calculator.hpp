#ifndef GRAPH_RANDOM_METRIC_CALCUlATOR_HPP
#define GRAPH_RANDOM_METRIC_CALCUlATOR_HPP

#include <algorithm>
#include <map>
#include <set>
#include <stack>
#include <vector>

#include "ba_graph.hpp"  // graph::random::BAGraph, graph::random::BANode, graph::random::BAEdge

namespace graph::random {
class MetricCalculator {
 public:
  enum class MetricType { ALL_PATHS, GIVEN_SIZE_PATHS, MAX, GIVEN_SIZE_MAX, MAX_FASTER, EXPECTED, SUM_OF_NEIGHBOURS_DEGREES };

  struct Metric {
    double k33{0.0F};
    double k5{0.0F};
  };

  typedef std::map<size_t, std::map<size_t, MetricCalculator::Metric>> MetricsMap;

  MetricCalculator() = default;

  MetricCalculator(const MetricCalculator&) = default;
  MetricCalculator(MetricCalculator&&) = default;

  MetricCalculator& operator=(const MetricCalculator&) = default;
  MetricCalculator& operator=(MetricCalculator&&) = default;

  ~MetricCalculator() = default;

  static void calculate(MetricType metric_type, const BAGraph& graph, MetricsMap& metrics_map);

 private:
  static void allPathsMetric(const BAGraph& graph, MetricsMap& metrics_map);
  static void givenSizePathsMetric(const BAGraph& graph, MetricsMap& metrics_map);
  static void maxMetric(const BAGraph& graph, MetricsMap& metrics_map);
  static void givenSizeMaxMetric(const BAGraph& graph, MetricsMap& metrics_map);
  static void maxMetricFaster(const BAGraph& graph, MetricsMap& metrics_map);
  static void expectedMetric(const BAGraph& graph, MetricsMap& metrics_map);
  static void sumOfNeighboursDegreesMetric(const BAGraph& graph, MetricsMap& metrics_map);

  static void getAllPathsUtil(const BAGraph& graph, size_t node, std::vector<bool>& visited, std::vector<size_t>& path,
                              const std::function<void(const std::vector<size_t>&)>& callback);

  static void getAllPathsWithMaxLengthUtil(const BAGraph& graph, size_t node, std::vector<bool>& visited,
                                           std::vector<size_t>& path, size_t max_length,
                                           const std::function<void(const std::vector<size_t>&)>& callback);

  [[nodiscard]] static std::map<size_t, std::vector<size_t>> enumerateComponents(const BAGraph& graph);

  [[nodiscard]] static size_t getNearestNoTwoDegreeNode(const BAGraph& graph, size_t node, std::set<size_t>& banned_nodes);

  template <size_t MAX_DEG>
  [[nodiscard]] static double calculatePathImpact(const BAGraph& graph, const std::vector<size_t>& path);

  [[nodiscard]] static double phi(size_t node_degree);

  template <size_t MAX_DEG>
  [[nodiscard]] static double psi(size_t node_degree);

  static constexpr size_t K33_MAX_DEG{3};
  static constexpr size_t K5_MAX_DEG{4};
};
}  // namespace graph::random

#endif  // GRAPH_RANDOM_METRIC_CALCUlATOR_HPP