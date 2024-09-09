#include "avoiding_kuratowski_graphs_painter.hpp"

#include <chrono>
#include <limits>  // std::numeric_limits
#include <queue>   // std::queue
#include <stack>   // std::stack
#include <vector>  // std::vector

namespace graph::random {
AvoidingKuratowskiGraphsPainter::AvoidingKuratowskiGraphsPainter(std::vector<ColorType> edges_colors)
    : GraphPainter({}, std::move(edges_colors)) {

  for (const auto& edges_color : edges_colors_) {
    embeddings_[edges_color] = BAGraph{};
    metrics_map_[edges_color] = MetricsMap{};
  }
}

void AvoidingKuratowskiGraphsPainter::paintNode(BAGraph& /*graph*/, BANode& node) {
  for (const auto& edges_color : edges_colors_) {
    embeddings_[edges_color].addNode();

    for (auto& iter : metrics_map_[edges_color]) {
      iter.second[node.id] = Metric{};
    }

    metrics_map_[edges_color][node.id] = std::map<size_t, Metric>{};
    for (size_t node_v = 0; node_v < node.id; ++node_v) {
      metrics_map_[edges_color][node.id][node_v] = Metric{};
    }
  }
}

void AvoidingKuratowskiGraphsPainter::paintEdge(BAGraph& graph, BAEdge& edge) {
  auto start_coloring = std::chrono::high_resolution_clock::now();

  std::map<ColorType, float> metrics;
  std::map<ColorType, MetricsMap> metrics_map_copy;
  for (const auto& edges_color : edges_colors_) {
    metrics_map_copy[edges_color] = MetricsMap{metrics_map_[edges_color]};

    BAGraph embedding_copy{embeddings_[edges_color]};
    embedding_copy.addEdge(edge.source, edge.target);

    //	updateShortestPathsMetric(embedding_copy, metrics_map_copy);
    updateAllPathsMetric(embedding_copy, metrics_map_copy[edges_color]);
    //	updateAllPathsMetricSmart(embedding_copy, metrics_map_copy);
    const Metric metric = productSubgraphMetric(metrics_map_copy[edges_color]);
//    const Metric metric = sumMetric(metrics_map_copy[edges_color]);
    metrics[edges_color] = metric.k33 + metric.k5;
  }

  float max_metric_value{std::numeric_limits<float>::max()};
  for (const auto& edges_color : edges_colors_) {
    if (max_metric_value > metrics[edges_color]) {
      max_metric_value = metrics[edges_color];
    }
  }

  std::vector<ColorType> max_colors;
  for (const auto& edges_color : edges_colors_) {
    if (std::fabs(max_metric_value - metrics[edges_color]) < EPSILON) {
      max_colors.push_back(edges_color);
    }
  }

  std::uniform_int_distribution<size_t> distribution(0, max_colors.size() - 1);
  const size_t index = distribution(generator_);
  embeddings_[max_colors[index]].addEdge(edge.source, edge.target);
  edge.color = max_colors[index];
  embeddings_[max_colors[index]].getEdge(edge.source, edge.target).color = max_colors[index];
  metrics_map_[max_colors[index]] = metrics_map_copy[max_colors[index]];

  auto end_coloring = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_coloring - start_coloring);
  coloring_times_vector_.emplace_back(graph.getEdgesNumber(), duration.count());

  std::cout << "Painting edge: " << edge.source << " " << edge.target << " to " << max_colors[index].hex << std::endl;
}

void AvoidingKuratowskiGraphsPainter::reset() {}

void AvoidingKuratowskiGraphsPainter::updateAllPathsMetric(BAGraph& graph,
                                                           AvoidingKuratowskiGraphsPainter::MetricsMap& metrics_map) {
  metrics_map.clear();
  for (size_t node_v = 0; node_v < graph.getNodesNumber(); ++node_v) {
    metrics_map[node_v] = std::map<size_t, Metric>();
    for (size_t node_u = 0; node_u < graph.getNodesNumber(); ++node_u) {
      metrics_map[node_v][node_u] = Metric{};
    }
  }

  for (size_t source = 0; source < graph.getNodesNumber(); ++source) {
    std::stack<std::vector<size_t>> paths;
    paths.push({source});

    while (!paths.empty()) {
      std::vector<size_t> path = paths.top();
      paths.pop();
      const size_t last = path.back();

      const float k33_factor = calculatePathImpact<K33_MAX_DEG>(graph, path) * 0.5F;
      metrics_map[path.front()][path.back()].k33 += k33_factor;
      metrics_map[path.back()][path.front()].k33 += k33_factor;

      const float k5_factor = calculatePathImpact<K5_MAX_DEG>(graph, path) * 0.5F;
      metrics_map[path.front()][path.back()].k5 += k5_factor;
      metrics_map[path.back()][path.front()].k5 += k5_factor;

      for (const auto& neighbor : graph.getNeighbours(last)) {
        if (std::find(path.begin(), path.end(), neighbor.id) == path.end()) {
          std::vector<size_t> newPath = path;
          newPath.push_back(neighbor.id);
          paths.push(newPath);
        }
      }
    }
  }
}

void AvoidingKuratowskiGraphsPainter::updateAllPathsMetricSmart(BAGraph& graph,
                                                                AvoidingKuratowskiGraphsPainter::MetricsMap& metrics_map) {
  MetricsMap new_metrics_map{metrics_map};
  const size_t mm_size{metrics_map.size()};

  const BAEdge& last_added_edge = graph.getLastAddedEdge();  // source is added earlier than target
  const size_t source{last_added_edge.target};
  const size_t target{last_added_edge.source};

  std::cout << "New edge: " << source << " -> " << target << std::endl;

  const std::vector<BANode> source_neighbours = graph.getNeighbours(source);
  const std::vector<BANode> target_neighbours = graph.getNeighbours(target);

  std::cout << "Source(" << source << ")"
            << " neighbours: " << source_neighbours.size() << std::endl;
  std::cout << "Target(" << target << ")"
            << " neighbours: " << target_neighbours.size() << std::endl;
  std::cout << std::endl;

  if (source_neighbours.size() == 1 && target_neighbours.size() == 1) {
    const float k33 = psi<K33_MAX_DEG>(1) * psi<K33_MAX_DEG>(1);
    const float k5 = psi<K5_MAX_DEG>(1) * psi<K5_MAX_DEG>(1);

    new_metrics_map[source][target].k33 = k33;
    new_metrics_map[target][source].k33 = k33;
    new_metrics_map[source][target].k5 = k5;
    new_metrics_map[target][source].k5 = k5;
  }
  else if (source_neighbours.size() == 1) {
    // update earlier calculated metrics of paths ended with target

    // expand metric for paths ended in target to source node
    for (size_t node = 0; node < mm_size; ++node) {
      if (node == target || node == source) {
        continue;
      }

      const float val = metrics_map[node][target].k33 / psi<K33_MAX_DEG>(graph.getDegree(target) - 1) *
                  psi<K33_MAX_DEG>(graph.getDegree(source)) * phi(graph.getDegree(target));
      new_metrics_map[node][source].k33 += val;
      new_metrics_map[source][node].k33 += val;
    }

    // expand metric for path source -> target
    const float val = psi<K33_MAX_DEG>(1) * psi<K33_MAX_DEG>(1);
    new_metrics_map[source][target].k33 += val;
    new_metrics_map[target][source].k33 += val;
  }
  else {}

  metrics_map = std::move(new_metrics_map);
}

void AvoidingKuratowskiGraphsPainter::updateShortestPathsMetric(BAGraph& graph,
                                                                AvoidingKuratowskiGraphsPainter::MetricsMap& metrics_map) {
  metrics_map.clear();
  for (size_t node_v = 0; node_v < graph.getNodesNumber(); ++node_v) {
    metrics_map[node_v] = std::map<size_t, Metric>();
    for (size_t node_u = 0; node_u < graph.getNodesNumber(); ++node_u) {
      metrics_map[node_v][node_u] = Metric{};
    }
  }

  for (size_t source = 0; source < graph.getNodesNumber(); ++source) {
    std::vector<size_t> distances(graph.getNodesNumber(), std::numeric_limits<size_t>::max());
    std::vector<std::vector<size_t>> predecessors(graph.getNodesNumber(), std::vector<size_t>());
    std::queue<size_t> queue;

    distances[source] = 0;
    queue.push(source);

    while (!queue.empty()) {
      const size_t node = queue.front();
      queue.pop();

      for (const auto& neighbor : graph.getNeighbours(node)) {
        if (distances[neighbor.id] == std::numeric_limits<size_t>::max()) {
          distances[neighbor.id] = distances[node] + 1;
          predecessors[neighbor.id].push_back(node);
          queue.push(neighbor.id);
        }
        else if (distances[neighbor.id] == distances[node] + 1) {
          predecessors[neighbor.id].push_back(node);
        }
      }
    }

    for (size_t target = 0; target < graph.getNodesNumber(); ++target) {
      std::stack<std::vector<size_t>> paths;
      paths.push({target});

      while (!paths.empty()) {
        std::vector<size_t> path = paths.top();
        paths.pop();
        const size_t last = path.back();

        if (last == source) {
          const float k33_factor = calculatePathImpact<K33_MAX_DEG>(graph, path) * 0.5F;
          metrics_map[path.front()][path.back()].k33 += k33_factor;
          metrics_map[path.back()][path.front()].k33 += k33_factor;

          const float k5_factor = calculatePathImpact<K5_MAX_DEG>(graph, path) * 0.5F;
          metrics_map[path.front()][path.back()].k5 += k5_factor;
          metrics_map[path.back()][path.front()].k5 += k5_factor;
        }
        else {
          for (const size_t pred : predecessors[last]) {
            std::vector<size_t> new_path = path;
            new_path.push_back(pred);
            paths.push(new_path);
          }
        }
      }
    }
  }
}

AvoidingKuratowskiGraphsPainter::Metric AvoidingKuratowskiGraphsPainter::sumMetric(
    const AvoidingKuratowskiGraphsPainter::MetricsMap& metrics_map) {
  float k33_sum{0.0F};
  float k5_sum{0.0F};

  for (const auto& outerPair : metrics_map) {
    for (const auto& innerPair : outerPair.second) {
      k33_sum += innerPair.second.k33;
      k5_sum += innerPair.second.k5;
    }
  }

  return {.k33 = k33_sum, .k5 = k5_sum};
}

AvoidingKuratowskiGraphsPainter::Metric AvoidingKuratowskiGraphsPainter::productSubgraphMetric(
    const AvoidingKuratowskiGraphsPainter::MetricsMap& metrics_map) {
  constexpr float THRESHOLD{0.25F};

  float k33_sum{0.0F};
  float k5_sum{0.0F};

  const size_t nodes_number{metrics_map.size()};

  // check for K3,3
  for (size_t node1 = 0; node1 < nodes_number; ++node1) {
    for (size_t node2 = node1 + 1; node2 < nodes_number; ++node2) {
      for (size_t node3 = node2 + 1; node3 < nodes_number; ++node3) {
        for (size_t node4 = 0; node4 < nodes_number; ++node4) {
          if (node1 == node4 || node2 == node4 || node3 == node4) {
            continue;
          }
          if (metrics_map.at(node1).at(node4).k33 < THRESHOLD || metrics_map.at(node2).at(node4).k33 < THRESHOLD ||
              metrics_map.at(node3).at(node4).k33 < THRESHOLD) {
            continue;
          }

          for (size_t node5 = node4 + 1; node5 < nodes_number; ++node5) {
            if (node1 == node5 || node2 == node5 || node3 == node5) {
              continue;
            }
            if (metrics_map.at(node1).at(node5).k33 < THRESHOLD || metrics_map.at(node2).at(node5).k33 < THRESHOLD ||
                metrics_map.at(node3).at(node5).k33 < THRESHOLD) {
              continue;
            }

            for (size_t node6 = node5 + 1; node6 < nodes_number; ++node6) {
              if (node1 == node6 || node2 == node6 || node3 == node6) {
                continue;
              }
              if (metrics_map.at(node1).at(node6).k33 < THRESHOLD || metrics_map.at(node2).at(node6).k33 < THRESHOLD ||
                  metrics_map.at(node3).at(node6).k33 < THRESHOLD) {
                continue;
              }

              k33_sum += 0.5F * metrics_map.at(node1).at(node4).k33 * metrics_map.at(node1).at((node5)).k33 *
                         metrics_map.at(node1).at(node6).k33 * metrics_map.at(node2).at(node4).k33 *
                         metrics_map.at(node2).at((node5)).k33 * metrics_map.at(node2).at(node6).k33 *
                         metrics_map.at(node3).at(node4).k33 * metrics_map.at(node3).at((node5)).k33 *
                         metrics_map.at(node3).at(node6).k33;
            }
          }
        }
      }
    }
  }

  // check for K5
  for (size_t node1 = 0; node1 < nodes_number; ++node1) {
    for (size_t node2 = node1 + 1; node2 < nodes_number; ++node2) {
      if (metrics_map.at(node1).at(node2).k5 < THRESHOLD) {
        continue;
      }

      for (size_t node3 = node2 + 1; node3 < nodes_number; ++node3) {
        if (metrics_map.at(node1).at(node3).k5 < THRESHOLD || metrics_map.at(node2).at(node3).k5 < THRESHOLD) {
          continue;
        }

        for (size_t node4 = node3 + 1; node4 < nodes_number; ++node4) {
          if (metrics_map.at(node1).at(node4).k5 < THRESHOLD || metrics_map.at(node2).at(node4).k5 < THRESHOLD ||
              metrics_map.at(node3).at(node4).k5 < THRESHOLD) {
            continue;
          }

          for (size_t node5 = node4 + 1; node5 < nodes_number; ++node5) {
            if (metrics_map.at(node1).at(node5).k5 < THRESHOLD || metrics_map.at(node2).at(node5).k5 < THRESHOLD ||
                metrics_map.at(node3).at(node5).k5 < THRESHOLD || metrics_map.at(node4).at(node5).k5 < THRESHOLD) {
              continue;
            }

            k5_sum += metrics_map.at(node1).at(node2).k5 * metrics_map.at(node1).at(node3).k5 *
                      metrics_map.at(node1).at(node4).k5 * metrics_map.at(node1).at(node5).k5 *
                      metrics_map.at(node2).at(node3).k5 * metrics_map.at(node2).at(node4).k5 *
                      metrics_map.at(node2).at(node5).k5 * metrics_map.at(node3).at(node4).k5 *
                      metrics_map.at(node3).at(node5).k5 * metrics_map.at(node4).at(node5).k5;
          }
        }
      }
    }
  }

  return {.k33 = k33_sum, .k5 = k5_sum};
}

template <size_t MAX_DEG>
float AvoidingKuratowskiGraphsPainter::calculatePathImpact(BAGraph& graph, std::vector<size_t>& path) {
  if (path.size() < 2) {
    return 0;
  }

  float result{1.0F};
  result *= psi<MAX_DEG>(graph.getDegree(path.front()));
  result *= psi<MAX_DEG>(graph.getDegree(path.back()));

  for (size_t i = 1; i < (path.size() - 1); ++i) {
    result *= phi(graph.getDegree(path.at(i)));
  }

  return result;
}

float AvoidingKuratowskiGraphsPainter::phi(size_t node_degree) {
  if (node_degree < 2) {
    return 1;
  }

  return 1.0F / static_cast<float>(node_degree - 1);
  //  return 2.0F / static_cast<float>(node_degree * (node_degree - 1));
  //  return 1;
}

template <size_t MAX_DEG>
float AvoidingKuratowskiGraphsPainter::psi(size_t node_degree) {
  //  return static_cast<float>(node_degree) / static_cast<float>(MAX_DEG);
  return static_cast<float>(std::min(node_degree, MAX_DEG)) / static_cast<float>(MAX_DEG);
}

size_t AvoidingKuratowskiGraphsPainter::binomial(size_t n_value, size_t k_value) {
  if (k_value > n_value) {
    return 0;
  }

  if (2 * k_value > n_value) {
    return binomial(n_value, n_value - k_value);
  }

  size_t result{n_value};
  for (size_t i = 1; i < k_value; ++i) {
    result *= (n_value - i);
    result /= (i + 1);
  }

  return result;
}

size_t AvoidingKuratowskiGraphsPainter::factorial(size_t n_value) {
  size_t result{1};
  for (size_t i = 1; i <= n_value; ++i) {
    result *= i;
  }

  return result;
}

void AvoidingKuratowskiGraphsPainter::iterateOverSubsets(const std::vector<size_t>& elements, size_t target_size,
                                                         const std::function<void(const std::vector<size_t>&)>& callback) {

  std::vector<bool> mask(elements.size());
  std::fill(mask.end() - static_cast<std::iterator_traits<std::vector<bool>::iterator>::difference_type>(target_size), mask.end(),
            true);

  do {
    std::vector<size_t> currentSubset;
    for (size_t i = 0; i < elements.size(); ++i) {
      if (mask[i]) {
        currentSubset.push_back(elements[i]);
      }
    }
    callback(currentSubset);
  } while (std::next_permutation(mask.begin(), mask.end()));
}
}  // namespace graph::random
