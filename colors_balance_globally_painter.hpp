#ifndef GRAPH_RANDOM_COLORS_BALANCE_GLOBALLY_PAINTER_HPP
#define GRAPH_RANDOM_COLORS_BALANCE_GLOBALLY_PAINTER_HPP

#include <random>    // std::mt19937, std::random_device, std::uniform_real_distribution, std::uniform_int_distribution

#include "graph_painter.hpp"    // graph::random::GraphPainter

namespace graph::random {
class ColorsBalanceGloballyPainter : public GraphPainter {
 public:
  ColorsBalanceGloballyPainter() = default;
  explicit ColorsBalanceGloballyPainter(std::vector<std::string> edges_colors);

  ColorsBalanceGloballyPainter(const ColorsBalanceGloballyPainter &) = default;
  ColorsBalanceGloballyPainter(ColorsBalanceGloballyPainter &&) = default;

  ColorsBalanceGloballyPainter &operator=(const ColorsBalanceGloballyPainter &) = default;
  ColorsBalanceGloballyPainter &operator=(ColorsBalanceGloballyPainter &&) = default;

  ~ColorsBalanceGloballyPainter() override = default;

  [[nodiscard]]
  inline const std::map<std::string, size_t> &getColorsHistogram() const { return colors_histogram_; }

  void paintNode(BAGraph &graph, BANode &node) override {}
  void paintEdge(BAGraph &graph, BAEdge &edge) override;

  void reset() override;

 private:
  std::mt19937 generator_{std::random_device{}()};
  std::map<std::string, size_t> colors_histogram_{};
};
} // namespace graph::random

#endif // GRAPH_RANDOM_COLORS_BALANCE_GLOBALLY_PAINTER_HPP
