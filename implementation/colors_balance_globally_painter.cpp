#include "colors_balance_globally_painter.hpp"

namespace graph::random {
ColorsBalanceGloballyPainter::ColorsBalanceGloballyPainter(std::vector<ColorType> edges_colors)
    : GraphPainter({}, std::move(edges_colors)) {
  for (const auto& color : edges_colors_) {
    colors_histogram_[color] = 0;
  }
}

void ColorsBalanceGloballyPainter::paintEdge(BAGraph& graph, BAEdge& edge) {
  if (graph.getEdgesNumber() == 0) {
    std::uniform_int_distribution<size_t> distribution{0, edges_colors_.size() - 1};
    const ColorType color = edges_colors_[distribution(generator_)];
    edge.color = color;
    ++colors_histogram_[color];
  }
  else {
    std::uniform_real_distribution<double> distribution{0.0F, 1.0F};
    const double r = distribution(generator_);
    double cumulative_probability = 0.0F;
    for (const auto& [color, count] : colors_histogram_) {
      const double probability = static_cast<double>(graph.getEdgesNumber() - count) /
                          static_cast<double>((edges_colors_.size() - 1) * graph.getEdgesNumber());
      cumulative_probability += probability;

      if (r < cumulative_probability) {
        edge.color = color;
        ++colors_histogram_[color];
        return;
      }
    }
  }
}

void ColorsBalanceGloballyPainter::reset() {
  for (auto& [color, count] : colors_histogram_) {
    count = 0;
  }
}
}  // namespace graph::random
