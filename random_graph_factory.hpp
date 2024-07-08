// Based on: https://barabasi.com/f/622.pdf

#ifndef GRAPH_RANDOM_GRAPH_FACTORY_HPP
#define GRAPH_RANDOM_GRAPH_FACTORY_HPP

#include <cstddef>            // std::size_t
#include <memory>                // std::unique_ptr
#include <string>            // std::string

#include "ba_graph.hpp"        // graph::random::BAGraph
#include "graph_painter.hpp"    // graph::random::GraphPainter

namespace graph::random {
class RandomGraphFactory {
 public:
  RandomGraphFactory() = default;

  RandomGraphFactory(const RandomGraphFactory &) = default;
  RandomGraphFactory(RandomGraphFactory &&) = default;

  RandomGraphFactory &operator=(const RandomGraphFactory &) = default;
  RandomGraphFactory &operator=(RandomGraphFactory &&) = default;

  ~RandomGraphFactory() = default;

  [[nodiscard]]
  static BAGraph createBarabasiAlbertWithPreferentialAttachmentRepeatedNodes(size_t initial_nodes_number,
																			 size_t final_nodes_number,
																			 size_t edges_per_new_node_number,
																			 GraphPainter *painter = nullptr);

  [[nodiscard]]
  static BAGraph createBarabasiAlbertWithPreferentialAttachmentBatageljBrandes(size_t initial_nodes_number,
																			   size_t final_nodes_number,
																			   size_t edges_per_new_node_number,
																			   GraphPainter *painter = nullptr);

  [[nodiscard]]
  static BAGraph createBarabasiAlbertWithLinkSelection(size_t initial_nodes_number,
													   size_t final_nodes_number,
													   size_t edges_per_new_node_number,
													   GraphPainter *painter = nullptr);
  [[nodiscard]]
  static BAGraph createBarabasiAlbertWithCopyingModel(size_t initial_nodes_number,
													  size_t final_nodes_number,
													  size_t edges_per_new_node_number,
													  float copy_probability,
													  GraphPainter *painter = nullptr);
};
} // namespace graph::random

#endif // GRAPH_RANDOM_GRAPH_FACTORY_HPP
