#ifndef GRAPH_RANDOM_BA_GRAPH_HPP
#define GRAPH_RANDOM_BA_GRAPH_HPP

#include <cstring>

#include "sparse_graph.hpp"

namespace graph::random {
typedef struct ColorType { char hex[8]; } ColorType;

static bool operator<(const ColorType& lhs, const ColorType& rhs) {
  return std::strncmp(lhs.hex, rhs.hex, sizeof(lhs.hex)) < 0;
}

struct BANode {
  size_t id{0};
  ColorType color{"#000000"};
};

struct BAEdge {
  size_t source{0};
  size_t target{0};
  ColorType color{"#000000"};
};

typedef SparseGraph<BANode, BAEdge> BAGraph;
} // namespace graph::random

#endif // GRAPH_RANDOM_BA_GRAPH_HPP
