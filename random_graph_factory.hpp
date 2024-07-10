// Based on: https://barabasi.com/f/622.pdf

#ifndef GRAPH_RANDOM_GRAPH_FACTORY_HPP
#define GRAPH_RANDOM_GRAPH_FACTORY_HPP

#include <cstddef>            	// std::size_t
#include <memory>             	// std::unique_ptr
#include <random>       		// std::mt19937, std::random_device, std::uniform_real_distribution, std::uniform_int_distribution

#include "ba_graph.hpp"       	// graph::random::BAGraph
#include "graph_painter.hpp"  	// graph::random::GraphPainter

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
  BAGraph createBarabasiAlbertWithPreferentialAttachmentRepeatedNodes(size_t initial_nodes_number,
																	  size_t final_nodes_number,
																	  size_t edges_per_new_node_number,
																	  GraphPainter *painter = nullptr);

  [[nodiscard]]
  BAGraph createBarabasiAlbertWithPreferentialAttachmentBatageljBrandes(size_t initial_nodes_number,
																		size_t final_nodes_number,
																		size_t edges_per_new_node_number,
																		GraphPainter *painter = nullptr);

  [[nodiscard]]
  BAGraph createBarabasiAlbertWithLinkSelection(size_t initial_nodes_number,
												size_t final_nodes_number,
												size_t edges_per_new_node_number,
												GraphPainter *painter = nullptr);
  [[nodiscard]]
  BAGraph createBarabasiAlbertWithCopyingModel(size_t initial_nodes_number,
											   size_t final_nodes_number,
											   size_t edges_per_new_node_number,
											   float copy_probability,
											   GraphPainter *painter = nullptr);

 private:
  std::mt19937 generator_{std::random_device{}()};
};
} // namespace graph::random

#endif // GRAPH_RANDOM_GRAPH_FACTORY_HPP
