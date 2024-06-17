#ifndef GRAPH_RANDOM_BA_GRAPH_HPP
#define GRAPH_RANDOM_BA_GRAPH_HPP

#include "sparse_graph.hpp"

namespace graph::random {
struct BANode {
  size_t id{0};
  std::string color{"#000000"};
};

struct BAEdge {
  size_t source{0};
  size_t target{0};
  std::string color{"#000000"};
};

typedef SparseGraph<BANode, BAEdge> BAGraph;
} // namespace graph::random

#endif // GRAPH_RANDOM_BA_GRAPH_HPP
