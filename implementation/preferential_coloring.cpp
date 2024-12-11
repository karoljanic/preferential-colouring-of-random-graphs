#include "preferential_coloring.hpp"

namespace graph::random {
std::vector<BAGraph> PreferentialColoring::color(const BAGraph& graph, size_t colors_count) {
  std::vector<BAGraph> subgraphs(colors_count);
  std::vector<DisjointSet> components(colors_count);
  for (size_t color = 0; color < colors_count; ++color) {
    components[color] = DisjointSet(graph.getNodesNumber());
    for (size_t node = 0; node < graph.getNodesNumber(); ++node) {
      components[color].makeSet(node);
    }
  }

  std::vector<size_t> sum_of_edge_degrees_products(colors_count, 0);

  std::vector<size_t> nodes_order(graph.getNodesNumber());
  std::iota(nodes_order.begin(), nodes_order.end(), 0);

  for (const auto node : nodes_order) {
    for (size_t color = 0; color < colors_count; ++color) {
      subgraphs[color].addNode();
    }

    for (const auto neighbour : graph.getNeighbours(node)) {
      if (neighbour.id < node) {
        std::vector<size_t> sum_of_edge_degrees_products_copy = sum_of_edge_degrees_products;

        for (size_t color = 0; color < colors_count; ++color) {
          for (const auto next_node : graph.getNeighbours(neighbour.id)) {
            sum_of_edge_degrees_products_copy[color] -=
                subgraphs[color].getDegree(next_node.id) * subgraphs[color].getDegree(neighbour.id);
            sum_of_edge_degrees_products_copy[color] +=
                subgraphs[color].getDegree(next_node.id) * (subgraphs[color].getDegree(neighbour.id) + 1);
          }

          for (const auto next_node : graph.getNeighbours(node)) {
            sum_of_edge_degrees_products_copy[color] -=
                subgraphs[color].getDegree(next_node.id) * subgraphs[color].getDegree(node);
            sum_of_edge_degrees_products_copy[color] +=
                subgraphs[color].getDegree(next_node.id) * (subgraphs[color].getDegree(node) + 1);
          }

          sum_of_edge_degrees_products_copy[color] +=
              (subgraphs[color].getDegree(node) + 1) * (subgraphs[color].getDegree(neighbour.id) + 1);
        }

        bool found = false;
        size_t new_color = 0;
        for (size_t color = 0; color < colors_count; ++color) {
          if (components[color].findSet(node) != components[color].findSet(neighbour.id)) {
            components[color].unionSets(node, neighbour.id);
            new_color = color;
            found = true;
            break;
          }
        }

        if (!found) {
          double min_metric = std::numeric_limits<double>::max();
          for (size_t color = 0; color < colors_count; ++color) {
            const double metric = static_cast<double>(sum_of_edge_degrees_products_copy[color]) *
                            (subgraphs[color].getEdgesNumber() + 1) * std::log(subgraphs[color].getNodesNumber() + 1);
            if (metric < min_metric) {
              min_metric = metric;
              new_color = color;
            }
          }
        }

        subgraphs[new_color].addEdge(node, neighbour.id);
        sum_of_edge_degrees_products[new_color] = sum_of_edge_degrees_products_copy[new_color];
      }
    }
  }

  return subgraphs;
}
}  // namespace graph::random