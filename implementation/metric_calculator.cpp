#include "metric_calculator.hpp"

namespace graph::random {
void MetricCalculator::calculate(MetricType metric_type, const BAGraph& graph, MetricsMap& metrics_map) {
  switch (metric_type) {
    case MetricType::ALL_PATHS:
      allPathsMetric(graph, metrics_map);
      break;
    case MetricType::GIVEN_SIZE_PATHS:
      givenSizePathsMetric(graph, metrics_map);
      break;
    case MetricType::MAX:
      std::cout << "Use MAX_FASTER instead of MAX metric\n";
      maxMetric(graph, metrics_map);
      break;
    case MetricType::GIVEN_SIZE_MAX:
      givenSizeMaxMetric(graph, metrics_map);
      break;
    case MetricType::MAX_FASTER:
      maxMetricFaster(graph, metrics_map);
      break;
    case MetricType::EXPECTED:
      expectedMetric(graph, metrics_map);
      break;
    case MetricType::GIVEN_SIZE_EXPECTED:
      givenSizeExpectedMetric(graph, metrics_map);
      break;
  }
}
void MetricCalculator::allPathsMetric(const BAGraph& graph, MetricsMap& metrics_map) {
  for (size_t node1 = 0; node1 < graph.getNodesNumber(); ++node1) {
    for (size_t node2 = 0; node2 < graph.getNodesNumber(); ++node2) {
      metrics_map[node1][node2] = Metric{};
    }
  }

  for (size_t source = 0; source < graph.getNodesNumber(); ++source) {
    std::vector<bool> visited(graph.getNodesNumber(), false);
    std::vector<size_t> path;

    getAllPathsUtil(graph, source, visited, path, [&](const std::vector<size_t>& path) {
      const float k33_factor = calculatePathImpact<K33_MAX_DEG>(graph, path);
      metrics_map[path.front()][path.back()].k33 += k33_factor;
      metrics_map[path.back()][path.front()].k33 += k33_factor;

      const float k5_factor = calculatePathImpact<K5_MAX_DEG>(graph, path);
      metrics_map[path.front()][path.back()].k5 += k5_factor;
      metrics_map[path.back()][path.front()].k5 += k5_factor;
    });
  }
}

void MetricCalculator::givenSizePathsMetric(const BAGraph& graph, MetricsMap& metrics_map) {
  for (size_t node1 = 0; node1 < graph.getNodesNumber(); ++node1) {
    for (size_t node2 = 0; node2 < graph.getNodesNumber(); ++node2) {
      metrics_map[node1][node2] = Metric{};
    }
  }

  const size_t max_length = 7;

  for (size_t source = 0; source < graph.getNodesNumber(); ++source) {
    std::vector<bool> visited(graph.getNodesNumber(), false);
    std::vector<size_t> path;

    getAllPathsWithMaxLengthUtil(graph, source, visited, path, max_length, [&](const std::vector<size_t>& path) {
      const float k33_factor = calculatePathImpact<K33_MAX_DEG>(graph, path);
      metrics_map[path.front()][path.back()].k33 += k33_factor;
      metrics_map[path.back()][path.front()].k33 += k33_factor;

      const float k5_factor = calculatePathImpact<K5_MAX_DEG>(graph, path);
      metrics_map[path.front()][path.back()].k5 += k5_factor;
      metrics_map[path.back()][path.front()].k5 += k5_factor;
    });
  }
}

void MetricCalculator::maxMetric(const BAGraph& graph, MetricsMap& metrics_map) {
  for (size_t node1 = 0; node1 < graph.getNodesNumber(); ++node1) {
    for (size_t node2 = 0; node2 < graph.getNodesNumber(); ++node2) {
      metrics_map[node1][node2] = Metric{};
    }
  }

  for (size_t source = 0; source < graph.getNodesNumber(); ++source) {
    std::vector<bool> visited(graph.getNodesNumber(), false);
    std::vector<size_t> path;

    getAllPathsUtil(graph, source, visited, path, [&](const std::vector<size_t>& path) {
      const float k33_factor = calculatePathImpact<K33_MAX_DEG>(graph, path);
      metrics_map[path.front()][path.back()].k33 = std::max(metrics_map[path.front()][path.back()].k33, k33_factor);
      metrics_map[path.back()][path.front()].k33 = std::max(metrics_map[path.back()][path.front()].k33, k33_factor);

      const float k5_factor = calculatePathImpact<K5_MAX_DEG>(graph, path);
      metrics_map[path.front()][path.back()].k5 = std::max(metrics_map[path.front()][path.back()].k5, k5_factor);
      metrics_map[path.back()][path.front()].k5 = std::max(metrics_map[path.back()][path.front()].k5, k5_factor);
    });
  }
}

void MetricCalculator::givenSizeMaxMetric(const BAGraph& graph, MetricsMap& metrics_map) {
  for (size_t node1 = 0; node1 < graph.getNodesNumber(); ++node1) {
    for (size_t node2 = 0; node2 < graph.getNodesNumber(); ++node2) {
      metrics_map[node1][node2] = Metric{};
    }
  }

  const size_t max_length = 7;

  for (size_t source = 0; source < graph.getNodesNumber(); ++source) {
    std::vector<bool> visited(graph.getNodesNumber(), false);
    std::vector<size_t> path;

    getAllPathsWithMaxLengthUtil(graph, source, visited, path, max_length, [&](const std::vector<size_t>& path) {
      const float k33_factor = calculatePathImpact<K33_MAX_DEG>(graph, path);
      metrics_map[path.front()][path.back()].k33 = std::max(metrics_map[path.front()][path.back()].k33, k33_factor);
      metrics_map[path.back()][path.front()].k33 = std::max(metrics_map[path.back()][path.front()].k33, k33_factor);

      const float k5_factor = calculatePathImpact<K5_MAX_DEG>(graph, path);
      metrics_map[path.front()][path.back()].k5 = std::max(metrics_map[path.front()][path.back()].k5, k5_factor);
      metrics_map[path.back()][path.front()].k5 = std::max(metrics_map[path.back()][path.front()].k5, k5_factor);
    });
  }
}

void MetricCalculator::maxMetricFaster(const BAGraph& graph, MetricsMap& metrics_map) {
  for (size_t node1 = 0; node1 < graph.getNodesNumber(); ++node1) {
    for (size_t node2 = 0; node2 < graph.getNodesNumber(); ++node2) {
      metrics_map[node1][node2] = Metric{};
    }
  }

  for (size_t node = 0; node < graph.getNodesNumber(); ++node) {
    metrics_map[node][node].k33 = psi<K33_MAX_DEG>(graph.getDegree(node));
    metrics_map[node][node].k5 = psi<K5_MAX_DEG>(graph.getDegree(node));
  }

  for (size_t repeat = 0; repeat < graph.getNodesNumber(); ++repeat) {
    for (size_t source = 0; source < graph.getNodesNumber(); ++source) {
      for (size_t destination = 0; destination < graph.getNodesNumber(); ++destination) {
        for (const auto& neighbor : graph.getNeighbours(destination)) {
          if (neighbor.id == source) {
            const float new_k33 = psi<K33_MAX_DEG>(graph.getDegree(source)) * psi<K33_MAX_DEG>(graph.getDegree(destination));
            metrics_map[source][destination].k33 = std::max(metrics_map[source][destination].k33, new_k33);
            metrics_map[destination][source].k33 = std::max(metrics_map[destination][source].k33, new_k33);

            const float new_k5 = psi<K5_MAX_DEG>(graph.getDegree(source)) * psi<K5_MAX_DEG>(graph.getDegree(destination));
            metrics_map[source][destination].k5 = std::max(metrics_map[source][destination].k5, new_k5);
            metrics_map[destination][source].k5 = std::max(metrics_map[destination][source].k5, new_k5);
          }
          else {
            const float new_k33 = metrics_map[source][neighbor.id].k33 / psi<K33_MAX_DEG>(graph.getDegree(neighbor.id)) *
                                  phi(graph.getDegree(destination)) * psi<K33_MAX_DEG>(graph.getDegree(destination));
            metrics_map[source][destination].k33 = std::max(metrics_map[source][destination].k33, new_k33);
            metrics_map[destination][source].k33 = std::max(metrics_map[destination][source].k33, new_k33);

            const float new_k5 = metrics_map[source][neighbor.id].k5 / psi<K5_MAX_DEG>(graph.getDegree(neighbor.id)) *
                                 phi(graph.getDegree(destination)) * psi<K5_MAX_DEG>(graph.getDegree(destination));
            metrics_map[source][destination].k5 = std::max(metrics_map[source][destination].k5, new_k5);
            metrics_map[destination][source].k5 = std::max(metrics_map[destination][source].k5, new_k5);
          }
        }
      }
    }
  }

  // generate all vertex to vertex shortest paths using queue
  std::map<size_t, std::map<size_t, size_t>> distances;
  for (size_t node1 = 0; node1 < graph.getNodesNumber(); ++node1) {
    for (size_t node2 = 0; node2 < graph.getNodesNumber(); ++node2) {
      distances[node1][node2] = std::numeric_limits<size_t>::max();
    }
  }

  for (size_t source = 0; source < graph.getNodesNumber(); ++source) {
    std::queue<size_t> queue;
    queue.push(source);
    distances[source][source] = 0;

    while (!queue.empty()) {
      const size_t current_node = queue.front();
      queue.pop();

      for (const auto& neighbor : graph.getNeighbours(current_node)) {
        if (distances[source][neighbor.id] > distances[source][current_node] + 1) {
          distances[source][neighbor.id] = distances[source][current_node] + 1;
          queue.push(neighbor.id);
        }
      }
    }
  }

  for (size_t source = 0; source < graph.getNodesNumber(); ++source) {
    for (size_t destination = 0; destination < graph.getNodesNumber(); ++destination) {
      metrics_map[source][destination].k33 *= static_cast<float>(distances[source][destination]);
      metrics_map[source][destination].k5 *= static_cast<float>(distances[source][destination]);
    }
  }
}

void MetricCalculator::expectedMetric(const BAGraph& graph, MetricsMap& metrics_map) {
  for (size_t node1 = 0; node1 < graph.getNodesNumber(); ++node1) {
    for (size_t node2 = 0; node2 < graph.getNodesNumber(); ++node2) {
      metrics_map[node1][node2] = Metric{};
    }
  }

  size_t degrees_sum = 0;
  for (size_t node = 0; node < graph.getNodesNumber(); ++node) {
    degrees_sum += graph.getDegree(node);
  }
  const float average_degree = static_cast<float>(degrees_sum) / static_cast<float>(graph.getNodesNumber());

  for (size_t node = 0; node < graph.getNodesNumber(); ++node) {
    if (graph.getDegree(node) == 0) {
      continue;
    }
    metrics_map[node][node].k33 =
        0.5F * average_degree * std::pow(graph.getDegree(node), -3.0F) * std::log(1.0F + psi<K33_MAX_DEG>(graph.getDegree(node)));
    metrics_map[node][node].k5 =
        0.5F * average_degree * std::pow(graph.getDegree(node), -3.0F) * std::log(1.0F + psi<K5_MAX_DEG>(graph.getDegree(node)));
  }

  for (const auto& edge : graph.getEdges()) {
    std::set<size_t> edge_neighbors;
    for (const auto& neighbor : graph.getNeighbours(edge.source)) {
      edge_neighbors.insert(neighbor.id);
    }
    for (const auto& neighbor : graph.getNeighbours(edge.target)) {
      edge_neighbors.insert(neighbor.id);
    }
    edge_neighbors.erase(edge.source);
    edge_neighbors.erase(edge.target);

    size_t sum_of_neighbors_degrees = 0;
    for (const auto& neighbor : edge_neighbors) {
      std::set<size_t> banned_nodes{edge.source, edge.target};
      sum_of_neighbors_degrees += getNearestNoTwoDegreeNode(graph, neighbor, banned_nodes);
    }

    float edge_impact = 1.0F / std::pow(0.5F * average_degree, 0.1F * graph.getNodesNumber()) -
                        sum_of_neighbors_degrees / (0.4F * graph.getNodesNumber());
    if (graph.getDegree(edge.source) != 0) {
      edge_impact *= 0.5F * std::log(1.0F + phi(graph.getDegree(edge.source)));
    }

    if (graph.getDegree(edge.target) != 0) {
      edge_impact *= 0.5F * std::log(1.0F + phi(graph.getDegree(edge.target)));
    }

    metrics_map[edge.source][edge.target].k33 += edge_impact;
    metrics_map[edge.target][edge.source].k33 += edge_impact;

    metrics_map[edge.source][edge.target].k5 += edge_impact;
    metrics_map[edge.target][edge.source].k5 += edge_impact;
  }
}

void MetricCalculator::givenSizeExpectedMetric(const BAGraph& /*graph*/, MetricsMap& /*metrics_map*/) {}

void MetricCalculator::getAllPathsUtil(const BAGraph& graph, size_t node, std::vector<bool>& visited, std::vector<size_t>& path,
                                       const std::function<void(const std::vector<size_t>&)>& callback) {
  visited[node] = true;
  path.push_back(node);

  callback(path);

  for (const auto& neighbor : graph.getNeighbours(node)) {
    if (!visited[neighbor.id]) {
      getAllPathsUtil(graph, neighbor.id, visited, path, callback);
    }
  }

  path.pop_back();
  visited[node] = false;
}

void MetricCalculator::getAllPathsWithMaxLengthUtil(const BAGraph& graph, size_t node, std::vector<bool>& visited,
                                                    std::vector<size_t>& path, size_t max_length,
                                                    const std::function<void(const std::vector<size_t>&)>& callback) {
  visited[node] = true;
  path.push_back(node);

  if (path.size() <= max_length) {
    callback(path);
  }

  if (path.size() < max_length) {
    for (const auto& neighbor : graph.getNeighbours(node)) {
      if (!visited[neighbor.id]) {
        getAllPathsWithMaxLengthUtil(graph, neighbor.id, visited, path, max_length, callback);
      }
    }
  }

  path.pop_back();
  visited[node] = false;
}

size_t MetricCalculator::getNearestNoTwoDegreeNode(const BAGraph& graph, size_t node, std::set<size_t>& banned_nodes) {
  return graph.getDegree(node);
  //  if (graph.getDegree(node) != 2) {
  //    return graph.getDegree(node);
  //  }
  //
  //  banned_nodes.insert(node);
  //  const auto& neighbors = graph.getNeighbours(node);
  //  if (banned_nodes.find(neighbors.front().id) == banned_nodes.end()) {
  //    return getNearestNoTwoDegreeNode(graph, neighbors.front().id, banned_nodes);
  //  }
  //  return getNearestNoTwoDegreeNode(graph, neighbors.back().id, banned_nodes);
}

template <size_t MAX_DEG>
float MetricCalculator::calculatePathImpact(const BAGraph& graph, const std::vector<size_t>& path) {
  if (path.size() < 2) {
    return 0;
  }

  // standard
  //  float result{1.0F};
  //  result *= psi<MAX_DEG>(graph.getDegree(path.front()));
  //  result *= psi<MAX_DEG>(graph.getDegree(path.back()));
  //
  //  for (size_t i = 1; i < (path.size() - 1); ++i) {
  //    result *= phi(graph.getDegree(path.at(i)));
  //  }
  //  return result;

  //   logarithmic
  float result{1.0F};
  result *= psi<MAX_DEG>(graph.getDegree(path.front()));
  result *= psi<MAX_DEG>(graph.getDegree(path.back()));

  for (size_t i = 1; i < (path.size() - 1); ++i) {
    result *= phi(graph.getDegree(path.at(i)));
  }
  return std::log(1.0F + result);

  // logarithmic with sum
  //  float result{0.0F};
  //  result += std::log(psi<MAX_DEG>(graph.getDegree(path.front())));
  //  result += std::log(psi<MAX_DEG>(graph.getDegree(path.back())));
  //
  //  for (size_t i = 1; i < (path.size() - 1); ++i) {
  //    result += std::log(phi(graph.getDegree(path.at(i))));
  //  }
  //
  //  return result;
}

float MetricCalculator::phi(size_t node_degree) {
  if (node_degree < 2) {
    return 1;
  }

  return 1.0F / static_cast<float>(node_degree - 1);
  //  return 2.0F / static_cast<float>(node_degree * (node_degree - 1));
}

template <size_t MAX_DEG>
float MetricCalculator::psi(size_t node_degree) {
  return static_cast<float>(std::min(node_degree, MAX_DEG)) / static_cast<float>(MAX_DEG);
  //  return 1.0F;
}

}  // namespace graph::random