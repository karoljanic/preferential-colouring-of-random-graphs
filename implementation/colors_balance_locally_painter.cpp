#include "colors_balance_locally_painter.hpp"

namespace graph::random {
ColorsBalanceLocallyPainter::ColorsBalanceLocallyPainter(std::vector<ColorType> edges_colors)
    : GraphPainter({}, std::move(edges_colors)) {
  for (const auto& color : edges_colors_) {
    colors_histogram_[color] = 0;
  }
}

void ColorsBalanceLocallyPainter::paintEdge(BAGraph& graph, BAEdge& edge) {
  std::map<ColorType, size_t> local_colors_histogram{};
  for (const auto& edges_color : edges_colors_) {
    local_colors_histogram[edges_color] = 0;
  }

  size_t total{0};
  for (const auto& adjacent_edge : graph.getAdjacentEdges(edge.source)) {
    ++local_colors_histogram[adjacent_edge.color];
    ++total;
  }
  for (const auto& adjacent_edge : graph.getAdjacentEdges(edge.target)) {
    ++local_colors_histogram[adjacent_edge.color];
    ++total;
  }

  if (total == 0) {
    std::uniform_int_distribution<size_t> distribution{0, edges_colors_.size() - 1};
    const ColorType color = edges_colors_[distribution(generator_)];
    edge.color = color;
    ++colors_histogram_[color];
  }
  else {
    std::uniform_real_distribution<double> distribution{0.0F, 1.0F};
    const double r = distribution(generator_);
    double cumulative_probability = 0.0F;
    for (const auto& color : edges_colors_) {
      const double probability =
          static_cast<double>(total - local_colors_histogram[color]) / static_cast<double>((edges_colors_.size() - 1) * total);
      cumulative_probability += probability;

      if (r < cumulative_probability) {
        edge.color = color;
        ++colors_histogram_[color];
        return;
      }
    }
  }
}

void ColorsBalanceLocallyPainter::reset() {
  for (auto& [color, count] : colors_histogram_) {
    count = 0;
  }
}
}  // namespace graph::random
