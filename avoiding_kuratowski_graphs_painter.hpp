#ifndef GRAPH_RANDOM_AVOIDING_KURATOWSKI_GRAPHS_PAINTER_HPP
#define GRAPH_RANDOM_AVOIDING_KURATOWSKI_GRAPHS_PAINTER_HPP

#include <random>    // std::mt19937, std::random_device, std::uniform_real_distribution, std::uniform_int_distribution

#include "graph_painter.hpp"    // graph::random::GraphPainter

namespace graph::random {
class AvoidingKuratowskiGraphsPainter : public GraphPainter {
 public:
  AvoidingKuratowskiGraphsPainter() = default;
  explicit AvoidingKuratowskiGraphsPainter(std::vector<std::string> edges_colors);

  AvoidingKuratowskiGraphsPainter(const AvoidingKuratowskiGraphsPainter &) = default;
  AvoidingKuratowskiGraphsPainter(AvoidingKuratowskiGraphsPainter &&) = default;

  AvoidingKuratowskiGraphsPainter &operator=(const AvoidingKuratowskiGraphsPainter &) = default;
  AvoidingKuratowskiGraphsPainter &operator=(AvoidingKuratowskiGraphsPainter &&) = default;

  ~AvoidingKuratowskiGraphsPainter() override = default;

  void paintNode(BAGraph &graph, BANode &node) override {}
  void paintEdge(BAGraph &graph, BAEdge &edge) override;

  void reset() override;

 private:
  std::map<std::string, BAGraph> embeddings_;
};
} // namespace graph::random

#endif // GRAPH_RANDOM_AVOIDING_KURATOWSKI_GRAPHS_PAINTER_HPP
