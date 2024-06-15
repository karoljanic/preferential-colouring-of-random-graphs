#ifndef GRAPH_RANDOM_GRAPH_FACTORY_HPP
#define GRAPH_RANDOM_GRAPH_FACTORY_HPP

#include "sparse_graph.hpp"

namespace graph::random {
struct BANode {
  size_t id;
};

struct BAEdge {
  size_t source;
  size_t target;
};

typedef SparseGraph<BANode, BAEdge> BAGraph;

class RandomGraphFactory {
 public:
  RandomGraphFactory() = default;

  RandomGraphFactory(const RandomGraphFactory &) = default;
  RandomGraphFactory(RandomGraphFactory &&) = default;

  RandomGraphFactory &operator=(const RandomGraphFactory &) = default;
  RandomGraphFactory &operator=(RandomGraphFactory &&) = default;

  ~RandomGraphFactory() = default;

  BAGraph createBarabasiAlbertWithPreferentialAttachment(size_t nodes_number, size_t edges_number);
  BAGraph createBarabasiAlbertWithLinkSelection(size_t nodes_number, size_t edges_number);
  BAGraph createBarabasiAlbertWithCopyingModel(size_t nodes_number, size_t edges_number);
};
} // namespace graph::random

#endif // GRAPH_RANDOM_GRAPH_FACTORY_HPP
