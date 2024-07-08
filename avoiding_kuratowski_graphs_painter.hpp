#ifndef GRAPH_RANDOM_AVOIDING_KURATOWSKI_GRAPHS_PAINTER_HPP
#define GRAPH_RANDOM_AVOIDING_KURATOWSKI_GRAPHS_PAINTER_HPP

#include <iostream> // std::cout
#include <random>    // std::mt19937, std::random_device, std::uniform_real_distribution, std::uniform_int_distribution
#include <utility>    // std::pair

#include "graph_painter.hpp"    // graph::random::GraphPainter

namespace graph::random {
class AvoidingKuratowskiGraphsPainter : public GraphPainter {
 public:
  struct NodeProperties {
	float k_5_probability = 0.0F;
	float k_3_3_probability = 0.0F;
  };

  AvoidingKuratowskiGraphsPainter() = default;
  explicit AvoidingKuratowskiGraphsPainter(std::vector<std::string> edges_colors);

  AvoidingKuratowskiGraphsPainter(const AvoidingKuratowskiGraphsPainter &) = default;
  AvoidingKuratowskiGraphsPainter(AvoidingKuratowskiGraphsPainter &&) = default;

  AvoidingKuratowskiGraphsPainter &operator=(const AvoidingKuratowskiGraphsPainter &) = default;
  AvoidingKuratowskiGraphsPainter &operator=(AvoidingKuratowskiGraphsPainter &&) = default;

  ~AvoidingKuratowskiGraphsPainter() override = default;

  void paintNode(BAGraph &graph, BANode &node) override;
  void paintEdge(BAGraph &graph, BAEdge &edge) override;

  void reset() override;

// private:
  typedef std::map<size_t, std::map<size_t, std::pair<float, float>>> MetricsMap;

  std::mt19937 generator_{std::random_device{}()};
  std::map<std::string, BAGraph> embeddings_;
  std::vector<NodeProperties> nodes_properties_;

  [[nodiscard]]
  static MetricsMap calculate_all_paths_metric(BAGraph &graph);
  [[nodiscard]]
  static MetricsMap calculate_shortest_paths_metric(BAGraph
													&graph);
  [[nodiscard]]
  static float calculate_path_impact(BAGraph &graph, std::vector<size_t> &path, int max_degree);
};
} // namespace graph::random

#endif // GRAPH_RANDOM_AVOIDING_KURATOWSKI_GRAPHS_PAINTER_HPP
