// Based on: https://barabasi.com/f/622.pdf

#ifndef GRAPH_RANDOM_GRAPH_FACTORY_HPP
#define GRAPH_RANDOM_GRAPH_FACTORY_HPP

#include <cstddef>      // std::size_t
#include <string>       // std::string

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

class RandomGraphFactory {
 public:
  RandomGraphFactory() = default;

  RandomGraphFactory(const RandomGraphFactory &) = default;
  RandomGraphFactory(RandomGraphFactory &&) = default;

  RandomGraphFactory &operator=(const RandomGraphFactory &) = default;
  RandomGraphFactory &operator=(RandomGraphFactory &&) = default;

  ~RandomGraphFactory() = default;

  [[nodiscard]]
  static BAGraph createBarabasiAlbertWithPreferentialAttachment(size_t initial_nodes_number,
																size_t final_nodes_number,
																size_t edges_per_new_node_number,
																float exponent_parameter);
  [[nodiscard]]
  static BAGraph createBarabasiAlbertWithLinkSelection(size_t initial_nodes_number,
													   size_t final_nodes_number,
													   size_t edges_per_new_node_number);
  [[nodiscard]]
  static BAGraph createBarabasiAlbertWithCopyingModel(size_t initial_nodes_number,
													  size_t final_nodes_number,
													  size_t edges_per_new_node_number,
													  float copy_probability);
};
} // namespace graph::random

#endif // GRAPH_RANDOM_GRAPH_FACTORY_HPP
