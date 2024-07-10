#ifndef GRAPH_RANDOM_COLORS_BALANCE_LOCALLY_PAINTER_HPP
#define GRAPH_RANDOM_COLORS_BALANCE_LOCALLY_PAINTER_HPP

#include <random>    			// std::mt19937, std::random_device, std::uniform_real_distribution, std::uniform_int_distribution

#include "graph_painter.hpp"    // graph::random::GraphPainter

namespace graph::random {
class ColorsBalanceLocallyPainter : public GraphPainter {
 public:
  ColorsBalanceLocallyPainter() = default;
  explicit ColorsBalanceLocallyPainter(std::vector<ColorType> edges_colors);

  ColorsBalanceLocallyPainter(const ColorsBalanceLocallyPainter &) = default;
  ColorsBalanceLocallyPainter(ColorsBalanceLocallyPainter &&) = default;

  ColorsBalanceLocallyPainter &operator=(const ColorsBalanceLocallyPainter &) = default;
  ColorsBalanceLocallyPainter &operator=(ColorsBalanceLocallyPainter &&) = default;

  ~ColorsBalanceLocallyPainter() override = default;

  [[nodiscard]]
  inline const std::map<ColorType, size_t> &getColorsHistogram() const { return colors_histogram_; }

  void paintNode(BAGraph &/*graph*/, BANode &/*node*/) override {}
  void paintEdge(BAGraph &/*graph*/, BAEdge &/*edge*/) override;

  void reset() override;
 private:
  std::mt19937 generator_{std::random_device{}()};
  std::map<ColorType , size_t> colors_histogram_{};
};
} // namespace graph::random

#endif // GRAPH_RANDOM_COLORS_BALANCE_LOCALLY_PAINTER_HPP
