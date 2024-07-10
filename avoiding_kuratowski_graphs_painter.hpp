#ifndef GRAPH_RANDOM_AVOIDING_KURATOWSKI_GRAPHS_PAINTER_HPP
#define GRAPH_RANDOM_AVOIDING_KURATOWSKI_GRAPHS_PAINTER_HPP

#include <iostream>                // std::cout
#include <random>                // std::mt19937, std::random_device, std::uniform_real_distribution, std::uniform_int_distribution
#include <utility>                // std::pair

#include "graph_painter.hpp"    // graph::random::GraphPainter

namespace graph::random {
class AvoidingKuratowskiGraphsPainter : public GraphPainter {
 public:
  AvoidingKuratowskiGraphsPainter() = default;
  explicit AvoidingKuratowskiGraphsPainter(std::vector<ColorType> edges_colors);

  AvoidingKuratowskiGraphsPainter(const AvoidingKuratowskiGraphsPainter &) = default;
  AvoidingKuratowskiGraphsPainter(AvoidingKuratowskiGraphsPainter &&) = default;

  AvoidingKuratowskiGraphsPainter &operator=(const AvoidingKuratowskiGraphsPainter &) = default;
  AvoidingKuratowskiGraphsPainter &operator=(AvoidingKuratowskiGraphsPainter &&) = default;

  ~AvoidingKuratowskiGraphsPainter() override = default;

  void paintNode(BAGraph &graph, BANode &node) override;
  void paintEdge(BAGraph &graph, BAEdge &edge) override;

  [[nodiscard]]
  const std::map<ColorType, BAGraph> &getEmbeddings() const { return embeddings_; }

  void reset() override;

// private:
  struct Metric {
	float k33{0.0F};
	float k5{0.0F};
  };

  typedef std::map<size_t, std::map<size_t, Metric>> MetricsMap;

  std::mt19937 generator_{std::random_device{}()};
  std::map<ColorType, BAGraph> embeddings_;
  std::map<ColorType, MetricsMap> metrics_map_;

  std::vector<std::pair<size_t, size_t>> coloring_times_vector_;

  static void updateAllPathsMetric(BAGraph &graph, MetricsMap &metrics_map);
  static void updateAllPathsMetricSmart(BAGraph &graph, MetricsMap &metrics_map);
  static void updateShortestPathsMetric(BAGraph &graph, MetricsMap &metrics_map);

  static void normalizeMetric(std::vector<Metric> &metric);

  [[nodiscard]]
  static Metric sumMetric(MetricsMap &metrics_map);

  [[nodiscard]]
  static float calculatePathImpact(BAGraph &graph, std::vector<size_t> &path, int max_degree);
};
} // namespace graph::random

#endif // GRAPH_RANDOM_AVOIDING_KURATOWSKI_GRAPHS_PAINTER_HPP
