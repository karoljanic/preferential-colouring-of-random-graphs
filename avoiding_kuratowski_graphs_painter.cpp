#include "avoiding_kuratowski_graphs_painter.hpp"

namespace graph::random {
AvoidingKuratowskiGraphsPainter::AvoidingKuratowskiGraphsPainter(std::vector<std::string> edges_colors)
	: GraphPainter({}, std::move(edges_colors)) {
}

void AvoidingKuratowskiGraphsPainter::paintEdge(BAGraph &graph, BAEdge &edge) {}

void AvoidingKuratowskiGraphsPainter::reset() {}
} // namespace graph::random
