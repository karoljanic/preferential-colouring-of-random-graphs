#include "random_graph_factory.hpp"

#include <random>       // std::mt19937, std::random_device, std::uniform_real_distribution, std::uniform_int_distribution
#include <stdexcept>    // std::invalid_argument

namespace graph::random {
BAGraph RandomGraphFactory::createBarabasiAlbertWithPreferentialAttachment(size_t initial_nodes_number,
																		   size_t final_nodes_number,
																		   size_t edges_per_new_node_number,
																		   float exponent_parameter) {
  if (edges_per_new_node_number < 1 || edges_per_new_node_number >= final_nodes_number) {
	throw std::invalid_argument("Edges per new node number must be greater than 0 and less than final nodes number");
  }

  if (initial_nodes_number < 2 || initial_nodes_number > final_nodes_number) {
	throw std::invalid_argument("Initial nodes number must be greater than 0 and less than final nodes number");
  }

  if (initial_nodes_number < edges_per_new_node_number) {
	throw std::invalid_argument("Initial nodes number must be greater than or equal to edges per new node number");
  }

  if (exponent_parameter <= 0.0F) {
	throw std::invalid_argument("Exponent parameter must be greater than 0");
  }

  BAGraph graph;
  BANode central_node{};
  graph.addNode(central_node);
  for (size_t i = 1; i < initial_nodes_number; ++i) {
	BANode node{};
	graph.addNode(node);
	BAEdge edge{node.id, central_node.id};
	graph.addEdge(edge);
  }

  std::mt19937 generator{std::random_device{}()};
  std::uniform_real_distribution<float> selection_distribution{0.0F, 1.0F};
  while (graph.getNodesNumber() < final_nodes_number) {
	BANode node{};
	graph.addNode(node);

	size_t edges_added{0};
	while (edges_added < edges_per_new_node_number) {
	  std::uniform_int_distribution<size_t> node_distribution{0, graph.getNodesNumber() - 1};
	  const BANode rand_node = graph.getNode(node_distribution(generator));
	  if (rand_node.id == node.id) {
		continue;
	  }

	  size_t total_degree = 2 * graph.getEdgesNumber();
	  float probability = std::pow(static_cast<float>(graph.getDegree(rand_node.id)) / static_cast<float>(total_degree),
								   exponent_parameter);
	  if (selection_distribution(generator) < probability) {
		BAEdge edge{node.id, rand_node.id};
		const bool edgeAdded = graph.addEdge(edge);

		if (edgeAdded) {
		  ++edges_added;
		}
	  }
	}
  }

  return graph;
}

BAGraph RandomGraphFactory::createBarabasiAlbertWithLinkSelection(size_t initial_nodes_number,
																  size_t final_nodes_number,
																  size_t edges_per_new_node_number) {
  if (edges_per_new_node_number < 1 || edges_per_new_node_number >= final_nodes_number) {
	throw std::invalid_argument("Edges per new node number must be greater than 0 and less than final nodes number");
  }

  if (initial_nodes_number < 2 || initial_nodes_number > final_nodes_number) {
	throw std::invalid_argument("Initial nodes number must be greater than 0 and less than final nodes number");
  }

  if (initial_nodes_number < edges_per_new_node_number) {
	throw std::invalid_argument("Initial nodes number must be greater than or equal to edges per new node number");
  }

  BAGraph graph;
  BANode central_node{};
  graph.addNode(central_node);
  for (size_t i = 1; i < initial_nodes_number; ++i) {
	BANode node{};
	graph.addNode(node);
	BAEdge edge{node.id, central_node.id};
	graph.addEdge(edge);
  }

  std::mt19937 generator{std::random_device{}()};
  while (graph.getNodesNumber() < final_nodes_number) {
	BANode node{};
	graph.addNode(node);

	size_t edges_added{0};
	while (edges_added < edges_per_new_node_number) {
	  std::uniform_int_distribution<size_t> node_distribution{0, graph.getNodesNumber() - 1};
	  const BANode rand_node = graph.getNode(node_distribution(generator));
	  if (rand_node.id == node.id) {
		continue;
	  }

	  const std::vector<BANode> rand_node_neighbours = graph.getNeighbours(rand_node.id);
	  std::uniform_int_distribution<size_t>
		  neighbour_distribution{0, rand_node_neighbours.size() - 1};
	  const BANode& rand_neighbour = rand_node_neighbours[neighbour_distribution(generator)];

	  BAEdge edge{node.id, rand_neighbour.id};
	  const bool edge_added = graph.addEdge(edge);

	  if (edge_added) {
		++edges_added;
	  }
	}
  }

  return graph;
}

BAGraph RandomGraphFactory::createBarabasiAlbertWithCopyingModel(size_t initial_nodes_number,
																 size_t final_nodes_number,
																 size_t edges_per_new_node_number,
																 float copy_probability) {
  if (edges_per_new_node_number < 1 || edges_per_new_node_number >= final_nodes_number) {
	throw std::invalid_argument("Edges per new node number must be greater than 0 and less than final nodes number");
  }

  if (initial_nodes_number < 2 || initial_nodes_number > final_nodes_number) {
	throw std::invalid_argument("Initial nodes number must be greater than 0 and less than final nodes number");
  }

  if (initial_nodes_number < edges_per_new_node_number) {
	throw std::invalid_argument("Initial nodes number must be greater than or equal to edges per new node number");
  }

  if (copy_probability < 0.0F || copy_probability > 1.0F) {
	throw std::invalid_argument("Copy probability must be between 0 and 1");
  }

  BAGraph graph;
  BANode central_node{};
  graph.addNode(central_node);
  for (size_t i = 1; i < initial_nodes_number; ++i) {
	BANode node{};
	graph.addNode(node);
	BAEdge edge{node.id, central_node.id};
	graph.addEdge(edge);
  }

  std::mt19937 generator{std::random_device{}()};
  std::uniform_real_distribution<float> option_distribution{0.0F, 1.0F};
  while (graph.getNodesNumber() < final_nodes_number) {
	BANode node{};
	graph.addNode(node);

	size_t edges_added{0};
	while (edges_added < edges_per_new_node_number) {
	  std::uniform_int_distribution<size_t> node_distribution{0, graph.getNodesNumber() - 1};
	  const BANode rand_node = graph.getNode(node_distribution(generator));
	  if (rand_node.id == node.id) {
		continue;
	  }

	  if (option_distribution(generator) < copy_probability) {
		BAEdge edge{node.id, rand_node.id};
		const bool edgeAdded = graph.addEdge(edge);

		if (edgeAdded) {
		  ++edges_added;
		}
	  } else {
		const std::vector<BANode> rand_node_neighbours = graph.getNeighbours(rand_node.id);
		std::uniform_int_distribution<size_t>
			neighbour_distribution{0, rand_node_neighbours.size() - 1};

		const BANode& rand_neighbour = rand_node_neighbours[neighbour_distribution(generator)];
		BAEdge edge{node.id, rand_neighbour.id};
		const bool edgeAdded = graph.addEdge(edge);

		if (edgeAdded) {
		  ++edges_added;
		}
	  }
	}
  }

  return graph;
}
} // namespace graph::random
