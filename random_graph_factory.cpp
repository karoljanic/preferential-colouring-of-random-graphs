#include "random_graph_factory.hpp"

#include <stdexcept>  // std::invalid_argument

namespace graph::random {
BAGraph RandomGraphFactory::createBarabasiAlbertWithPreferentialAttachmentRepeatedNodes(size_t initial_nodes_number,
                                                                                        size_t final_nodes_number,
                                                                                        size_t edges_per_new_node_number,
                                                                                        GraphPainter* painter) {
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

  size_t central_node_id = graph.addNode();
  painter->paintNode(graph, graph.getNode(central_node_id));
  for (size_t i = 1; i < initial_nodes_number; ++i) {
    size_t node_id = graph.addNode();
    painter->paintNode(graph, graph.getNode(node_id));

    graph.addEdge(node_id, central_node_id);
    painter->paintEdge(graph, graph.getEdge(node_id, central_node_id));

    repeated_nodes.push_back(node_id);
    repeated_nodes.push_back(central_node_id);
  }

  while (graph.getNodesNumber() < final_nodes_number) {
    size_t node_id = graph.addNode();
    painter->paintNode(graph, graph.getNode(node_id));

    size_t edges_added{0};
    std::vector<size_t> selected_nodes;
    std::uniform_int_distribution<size_t> repeated_node_distribution{0, repeated_nodes.size() - 1};
    while (edges_added < edges_per_new_node_number) {
      size_t selected_node_id = repeated_nodes[repeated_node_distribution(generator_)];
      if (std::find(selected_nodes.begin(), selected_nodes.end(), selected_node_id) == selected_nodes.end()) {
        selected_nodes.push_back(selected_node_id);
        if (!graph.edgeExists(node_id, selected_node_id)) {
          graph.addEdge(node_id, selected_node_id);
          painter->paintEdge(graph, graph.getEdge(node_id, selected_node_id));

          ++edges_added;

          repeated_nodes.push_back(node_id);
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
                                                                                          GraphPainter* painter) {
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

  // initialize n0 connected nodes
  for (size_t v = 0; v < initial_nodes_number; ++v) {
    nodes.push_back(v);
    nodes.push_back((v + 1) % initial_nodes_number);
    degrees[v]++;
    degrees[(v + 1) % initial_nodes_number]++;
  }

  for (size_t v = initial_nodes_number; v < final_nodes_number; ++v) {
    // If we were to update the range in the next loop, the additionally available nodes
    // would only lead to self-loops are multi-edges.
    std::uniform_int_distribution<size_t> dist{0, nodes.size() - 1};
    auto firstNeighbor = nodes.size() + 1;

    for (size_t i = 0; i < edges_per_new_node_number; ++i) {
      // let's sample a new neighbor and repeat if we're already connected to it
      while (true) {
        const auto randomIndex = dist(generator_);
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
          nodes.push_back(v);
          nodes.push_back(newNeighbor);
          degrees[v]++;
          degrees[newNeighbor]++;
          break;
        }
      }
    }
  }

  BAGraph graph;
  size_t node_id = graph.addNode();
  if (painter != nullptr) {
    painter->paintNode(graph, graph.getNode(node_id));
  }
  for (size_t i = 0; i < nodes.size(); i += 2) {
    while (node_id < nodes[i] || node_id < nodes[i + 1]) {
      node_id = graph.addNode();
      if (painter != nullptr) {
        painter->paintNode(graph, graph.getNode(graph.getNodesNumber() - 1));
      }
    }

    graph.addEdge(nodes[i], nodes[i + 1]);
    if (painter != nullptr) {
      painter->paintEdge(graph, graph.getEdge(nodes[i], nodes[i + 1]));
    }
  }

  return graph;
}

BAGraph RandomGraphFactory::createBarabasiAlbertWithLinkSelection(size_t initial_nodes_number, size_t final_nodes_number,
                                                                  size_t edges_per_new_node_number, GraphPainter* painter) {
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
  size_t central_node_id = graph.addNode();
  painter->paintNode(graph, graph.getNode(central_node_id));
  for (size_t i = 1; i < initial_nodes_number; ++i) {
    size_t node_id = graph.addNode();
    painter->paintNode(graph, graph.getNode(node_id));

    graph.addEdge(node_id, central_node_id);
    painter->paintEdge(graph, graph.getEdge(node_id, central_node_id));
  }

  while (graph.getNodesNumber() < final_nodes_number) {
    size_t node_id = graph.addNode();
    painter->paintNode(graph, graph.getNode(node_id));

    size_t edges_added{0};
    while (edges_added < edges_per_new_node_number) {
      std::uniform_int_distribution<size_t> node_distribution{0, graph.getNodesNumber() - 1};
      const BANode& rand_node = graph.getNode(node_distribution(generator_));
      if (rand_node.id == node_id) {
        continue;
      }

      const std::vector<BANode>& rand_node_neighbours = graph.getNeighbours(rand_node.id);
      std::uniform_int_distribution<size_t> neighbour_distribution{0, rand_node_neighbours.size() - 1};
      const BANode& rand_neighbour = rand_node_neighbours[neighbour_distribution(generator_)];

      if (!graph.edgeExists(node_id, rand_neighbour.id)) {
        graph.addEdge(node_id, rand_neighbour.id);
        painter->paintEdge(graph, graph.getEdge(node_id, rand_neighbour.id));
        ++edges_added;
      }
    }
  }

  return graph;
}

BAGraph RandomGraphFactory::createBarabasiAlbertWithCopyingModel(size_t initial_nodes_number, size_t final_nodes_number,
                                                                 size_t edges_per_new_node_number, float copy_probability,
                                                                 GraphPainter* painter) {
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
  size_t central_node_id = graph.addNode();
  painter->paintNode(graph, graph.getNode(central_node_id));
  for (size_t i = 1; i < initial_nodes_number; ++i) {
    size_t node_id = graph.addNode();
    painter->paintNode(graph, graph.getNode(node_id));

    graph.addEdge(node_id, central_node_id);
    painter->paintEdge(graph, graph.getEdge(node_id, central_node_id));
  }

  std::uniform_real_distribution<float> option_distribution{0.0F, 1.0F};
  while (graph.getNodesNumber() < final_nodes_number) {
    size_t node_id = graph.addNode();
    painter->paintNode(graph, graph.getNode(node_id));

    size_t edges_added{0};
    while (edges_added < edges_per_new_node_number) {
      std::uniform_int_distribution<size_t> node_distribution{0, graph.getNodesNumber() - 1};
      const BANode& rand_node = graph.getNode(node_distribution(generator_));
      if (rand_node.id == node_id) {
        continue;
      }

      if (option_distribution(generator_) < copy_probability) {
        if (!graph.edgeExists(node_id, rand_node.id)) {
          graph.addEdge(node_id, rand_node.id);
          painter->paintEdge(graph, graph.getEdge(node_id, rand_node.id));
          ++edges_added;
        }
      }
      else {
        const std::vector<BANode>& rand_node_neighbours = graph.getNeighbours(rand_node.id);
        std::uniform_int_distribution<size_t> neighbour_distribution{0, rand_node_neighbours.size() - 1};

        const BANode& rand_neighbour = rand_node_neighbours[neighbour_distribution(generator_)];
        if (!graph.edgeExists(node_id, rand_neighbour.id)) {
          graph.addEdge(node_id, rand_neighbour.id);
          painter->paintEdge(graph, graph.getEdge(node_id, rand_neighbour.id));
          ++edges_added;
        }
      }
    }
  }

  return graph;
}
}  // namespace graph::random
