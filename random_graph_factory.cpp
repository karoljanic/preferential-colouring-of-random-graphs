#include "random_graph_factory.hpp"

#include <random>       // std::mt19937, std::random_device, std::uniform_real_distribution, std::uniform_int_distribution
#include <stdexcept>    // std::invalid_argument

#include <iostream>

namespace graph::random {
BAGraph RandomGraphFactory::createBarabasiAlbertWithPreferentialAttachmentRepeatedNodes(size_t initial_nodes_number,
																						size_t final_nodes_number,
																						size_t edges_per_new_node_number,
																						GraphPainter *painter) {
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
  std::vector<size_t> repeated_nodes;

  BANode central_node{};
  painter->paintNode(graph, central_node);
  graph.addNode(central_node);
  for (size_t i = 1; i < initial_nodes_number; ++i) {
	BANode node{};
	graph.addNode(node);
	painter->paintNode(graph, node);

	BAEdge edge{node.id, central_node.id};
	graph.addEdge(edge);
	painter->paintEdge(graph, edge);

	repeated_nodes.push_back(node.id);
	repeated_nodes.push_back(central_node.id);
  }

  std::mt19937 generator{std::random_device{}()};
  while (graph.getNodesNumber() < final_nodes_number) {
	BANode node{};
	graph.addNode(node);
	painter->paintNode(graph, node);

	size_t edges_added{0};
	std::vector<size_t> selected_nodes;
	std::uniform_int_distribution<size_t> repeated_node_distribution{0, repeated_nodes.size() - 1};
	while (edges_added < edges_per_new_node_number) {
	  size_t selected_node_id = repeated_nodes[repeated_node_distribution(generator)];
	  if (std::find(selected_nodes.begin(), selected_nodes.end(), selected_node_id) == selected_nodes.end()) {
		selected_nodes.push_back(selected_node_id);
		BAEdge edge{node.id, selected_node_id};
		if (!graph.edgeExists(edge)) {
		  graph.addEdge(edge);
		  painter->paintEdge(graph, edge);

		  ++edges_added;

		  repeated_nodes.push_back(node.id);
		  repeated_nodes.push_back(selected_node_id);
		}
	  }
	}
  }

  return graph;
}

BAGraph RandomGraphFactory::createBarabasiAlbertWithPreferentialAttachmentBatageljBrandes(size_t initial_nodes_number,
																						  size_t final_nodes_number,
																						  size_t edges_per_new_node_number,
																						  GraphPainter *painter) {
  if (edges_per_new_node_number < 1 || edges_per_new_node_number >= final_nodes_number) {
	throw std::invalid_argument("Edges per new node number must be greater than 0 and less than final nodes number");
  }

  if (initial_nodes_number < 2 || initial_nodes_number > final_nodes_number) {
	throw std::invalid_argument("Initial nodes number must be greater than 0 and less than final nodes number");
  }

  if (initial_nodes_number < edges_per_new_node_number) {
	throw std::invalid_argument("Initial nodes number must be greater than or equal to edges per new node number");
  }

  std::vector<size_t> nodes;
  nodes.reserve(2 * initial_nodes_number + 2 * (final_nodes_number - initial_nodes_number) * edges_per_new_node_number);

  std::vector<size_t> degrees(final_nodes_number, 0);

  auto addEdge = [&](size_t u, size_t v) {
	nodes.push_back(u);
	nodes.push_back(v);
	degrees[u]++;
	degrees[v]++;
  };

  BAGraph graph;
  for(size_t v = 0; v < final_nodes_number; ++v) {
	BANode node{};
	graph.addNode(node);
  }

  // initialize n0 connected nodes
  for (size_t v = 0; v < initial_nodes_number - 1; ++v) {
	addEdge(v, v + 1);
  }
  addEdge(0, initial_nodes_number - 1);


  std::mt19937 generator{std::random_device{}()};
  for (size_t v = initial_nodes_number; v < final_nodes_number; ++v) {
	// If we were to update the range in the next loop, the additionally available nodes
	// would only lead to self-loops are multi-edges.
	std::uniform_int_distribution<size_t> dist{0, nodes.size() - 1};
	auto firstNeighbor = nodes.size() + 1;

	for (size_t i = 0; i < edges_per_new_node_number; ++i) {
	  // let's sample a new neighbor and repeat if we're already connected to it
	  while (true) {
		const auto randomIndex = dist(generator);
		const auto newNeighbor = nodes[randomIndex];

		// the last 2*(i-1) positions contain all edges incident to v in to format
		//  Even  Odd   Even  Odd   Even  Odd
		// | v | Neigh | v | Neigh | v | Neigh ...
		// Hence, we need to compare the new neighbor to the previous (i-1) odd positions
		bool alreadyIncident = false;
		for (auto j = firstNeighbor; j < nodes.size(); j += 2) {
		  if (nodes[j] == newNeighbor) {
			alreadyIncident = true;
			break;
		  }
		}

		if (!alreadyIncident) {
		  addEdge(v, newNeighbor);
		  break;
		}
	  }
	}
  }

  std::cout << "11\n";
  std::cout << graph.getNodesNumber() << std::endl;

  for (size_t i = 0; i < nodes.size(); i += 2) {
	BAEdge edge{.source = nodes[i], .target = nodes[i+1]};
	graph.addEdge(edge);
  }

  return graph;
}

BAGraph RandomGraphFactory::createBarabasiAlbertWithLinkSelection(size_t initial_nodes_number,
																  size_t final_nodes_number,
																  size_t edges_per_new_node_number,
																  GraphPainter *painter) {
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
  painter->paintNode(graph, central_node);
  graph.addNode(central_node);
  for (size_t i = 1; i < initial_nodes_number; ++i) {
	BANode node{};
	graph.addNode(node);
	painter->paintNode(graph, node);

	BAEdge edge{node.id, central_node.id};
	painter->paintEdge(graph, edge);
	graph.addEdge(edge);
  }

  std::mt19937 generator{std::random_device{}()};
  while (graph.getNodesNumber() < final_nodes_number) {
	BANode node{};
	graph.addNode(node);
	painter->paintNode(graph, node);

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
	  const BANode &rand_neighbour = rand_node_neighbours[neighbour_distribution(generator)];

	  BAEdge edge{node.id, rand_neighbour.id};
	  if (!graph.edgeExists(edge)) {
		painter->paintEdge(graph, edge);
		graph.addEdge(edge);
		++edges_added;
	  }
	}
  }

  return graph;
}

BAGraph RandomGraphFactory::createBarabasiAlbertWithCopyingModel(size_t initial_nodes_number,
																 size_t final_nodes_number,
																 size_t edges_per_new_node_number,
																 float copy_probability,
																 GraphPainter *painter) {
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
  painter->paintNode(graph, central_node);
  graph.addNode(central_node);
  for (size_t i = 1; i < initial_nodes_number; ++i) {
	BANode node{};
	graph.addNode(node);
	painter->paintNode(graph, node);

	BAEdge edge{node.id, central_node.id};
	painter->paintEdge(graph, edge);
	graph.addEdge(edge);
  }

  std::mt19937 generator{std::random_device{}()};
  std::uniform_real_distribution<float> option_distribution{0.0F, 1.0F};
  while (graph.getNodesNumber() < final_nodes_number) {
	BANode node{};
	graph.addNode(node);
	painter->paintNode(graph, node);

	size_t edges_added{0};
	while (edges_added < edges_per_new_node_number) {
	  std::uniform_int_distribution<size_t> node_distribution{0, graph.getNodesNumber() - 1};
	  const BANode rand_node = graph.getNode(node_distribution(generator));
	  if (rand_node.id == node.id) {
		continue;
	  }

	  if (option_distribution(generator) < copy_probability) {
		BAEdge edge{node.id, rand_node.id};
		if (!graph.edgeExists(edge)) {
		  painter->paintEdge(graph, edge);
		  graph.addEdge(edge);
		  ++edges_added;
		}
	  } else {
		const std::vector<BANode> rand_node_neighbours = graph.getNeighbours(rand_node.id);
		std::uniform_int_distribution<size_t>
			neighbour_distribution{0, rand_node_neighbours.size() - 1};

		const BANode &rand_neighbour = rand_node_neighbours[neighbour_distribution(generator)];
		BAEdge edge{node.id, rand_neighbour.id};
		if (!graph.edgeExists(edge)) {
		  painter->paintEdge(graph, edge);
		  graph.addEdge(edge);
		  ++edges_added;
		}
	  }
	}
  }

  return graph;
}
} // namespace graph::random
