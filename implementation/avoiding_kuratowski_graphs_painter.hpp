#ifndef GRAPH_RANDOM_AVOIDING_KURATOWSKI_GRAPHS_PAINTER_HPP
#define GRAPH_RANDOM_AVOIDING_KURATOWSKI_GRAPHS_PAINTER_HPP

#include <functional>  // std::function
#include <iostream>    // std::cout
#include <random>      // std::mt19937, std::random_device, std::uniform_real_distribution, std::uniform_int_distribution
#include <utility>     // std::pair

#include "graph_painter.hpp"  // graph::random::GraphPainter

namespace graph::random {
class AvoidingKuratowskiGraphsPainter : public GraphPainter {
 public:
  struct Metric {
    float k33{0.0F};
    float k5{0.0F};
  };

  struct TimeInfo {
    size_t edges_number;
    float total_time;
    float metrics_calculation_time;
    float metrics_merge_time;
  };

  typedef std::map<size_t, std::map<size_t, Metric>> MetricsMap;

  AvoidingKuratowskiGraphsPainter() = default;
  explicit AvoidingKuratowskiGraphsPainter(std::vector<ColorType> edges_colors);

  AvoidingKuratowskiGraphsPainter(const AvoidingKuratowskiGraphsPainter&) = default;
  AvoidingKuratowskiGraphsPainter(AvoidingKuratowskiGraphsPainter&&) = default;

  AvoidingKuratowskiGraphsPainter& operator=(const AvoidingKuratowskiGraphsPainter&) = default;
  AvoidingKuratowskiGraphsPainter& operator=(AvoidingKuratowskiGraphsPainter&&) = default;

  ~AvoidingKuratowskiGraphsPainter() override = default;

  void paintNode(BAGraph& graph, BANode& node) override;
  void paintEdge(BAGraph& graph, BAEdge& edge) override;

  [[nodiscard]] const std::map<ColorType, BAGraph>& getEmbeddings() const { return embeddings_; }
  [[nodiscard]] const std::vector<TimeInfo>& getColoringTimes() const { return coloring_times_vector_; }
  [[nodiscard]] const std::map<ColorType, MetricsMap>& getMetricsMap() const { return metrics_map_; }

  void reset() override;

 private:
  static constexpr size_t K33_MAX_DEG{3};
  static constexpr size_t K5_MAX_DEG{4};
  static constexpr float EPSILON{1e-6F};

  std::mt19937 generator_{std::random_device{}()};
  std::map<ColorType, BAGraph> embeddings_;
  std::map<ColorType, MetricsMap> metrics_map_;

  std::vector<TimeInfo> coloring_times_vector_;

  static void updateAllPathsMetric(BAGraph& graph, MetricsMap& metrics_map);
  static void updateAllPathsMetricSmart(BAGraph& graph, MetricsMap& metrics_map);
  static void updateShortestPathsMetric(BAGraph& graph, MetricsMap& metrics_map);

  static void getAllPathsUtil(BAGraph& graph, size_t node, std::vector<bool>& visited, std::vector<size_t>& path,
                              const std::function<void(const std::vector<size_t>&)>& callback);

  static void getAllPathsWithMaxLengthUtil(BAGraph& graph, size_t node, std::vector<bool>& visited, std::vector<size_t>& path,
                                           size_t max_length, const std::function<void(const std::vector<size_t>&)>& callback);

  [[nodiscard]] static Metric sumMetric(const MetricsMap& metrics_map);
  [[nodiscard]] static Metric productSubgraphMetric(const MetricsMap& metrics_map);

  template <size_t MAX_DEG>
  [[nodiscard]] static float calculatePathImpact(BAGraph& graph, const std::vector<size_t>& path);

  [[nodiscard]] static float phi(size_t node_degree);

  template <size_t MAX_DEG>
  [[nodiscard]] static float psi(size_t node_degree);

  static size_t binomial(size_t n_value, size_t k_value);
  static size_t factorial(size_t n_value);

  static void iterateOverSubsets(const std::vector<size_t>& elements, size_t target_size,
                                 const std::function<void(const std::vector<size_t>&)>& callback);
};
}  // namespace graph::random

#endif  // GRAPH_RANDOM_AVOIDING_KURATOWSKI_GRAPHS_PAINTER_HPP
