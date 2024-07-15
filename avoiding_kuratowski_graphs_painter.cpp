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

void AvoidingKuratowskiGraphsPainter::paintNode(BAGraph& graph, BANode& node) {
  for (const auto& edges_color : edges_colors_) {
    embeddings_[edges_color].addNode();

    for (auto& it : metrics_map_[edges_color]) {
      it.second[node.id] = Metric{};
    }

    metrics_map_[edges_color][node.id] = std::map<size_t, Metric>{};
    for (size_t v = 0; v < node.id; ++v) {
      metrics_map_[edges_color][node.id][v] = Metric{};
    }
  }
}

void AvoidingKuratowskiGraphsPainter::paintEdge(BAGraph& graph, BAEdge& edge) {
  auto start_coloring = std::chrono::high_resolution_clock::now();

  std::map<ColorType, float> metrics;
  for (const auto& edges_color : edges_colors_) {
    BAGraph embedding_copy{embeddings_[edges_color]};
    MetricsMap metrics_map_copy{metrics_map_[edges_color]};

    embedding_copy.addEdge(edge.source, edge.target);
    //	updateShortestPathsMetric(embedding_copy, metrics_map_copy);
    updateAllPathsMetric(embedding_copy, metrics_map_copy);
    //	updateAllPathsMetricSmart(embedding_copy, metrics_map_copy);
    Metric metric = sumMetric(metrics_map_copy);
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
    if (std::fabs(max_metric_value - metrics[edges_color]) < 0.0001F) {
      max_colors.push_back(edges_color);
    }
  }

  std::uniform_int_distribution<size_t> distribution(0, max_colors.size() - 1);
  size_t index = distribution(generator_);
  embeddings_[max_colors[index]].addEdge(edge.source, edge.target);
  edge.color = max_colors[index];
  embeddings_[max_colors[index]].getEdge(edge.source, edge.target).color = max_colors[index];

  auto end_coloring = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_coloring - start_coloring);
  coloring_times_vector_.emplace_back(graph.getEdgesNumber(), duration.count());

  std::cout << "Painting edge: " << edge.source << " " << edge.target << " to " << max_colors[index].hex << std::endl;
}

void AvoidingKuratowskiGraphsPainter::reset() {}

void AvoidingKuratowskiGraphsPainter::updateAllPathsMetric(BAGraph& graph,
                                                           AvoidingKuratowskiGraphsPainter::MetricsMap& metrics_map) {
  metrics_map.clear();
  for (size_t v = 0; v < graph.getNodesNumber(); ++v) {
    metrics_map[v] = std::map<size_t, Metric>();
    for (size_t u = 0; u < graph.getNodesNumber(); ++u) {
      metrics_map[v][u] = Metric{};
    }
  }

  for (size_t source = 0; source < graph.getNodesNumber(); ++source) {
    std::stack<std::vector<size_t>> paths;
    paths.push({source});

    while (!paths.empty()) {
      std::vector<size_t> path = paths.top();
      paths.pop();
      size_t last = path.back();

      float k33_factor = calculatePathImpact<K33_MAX_DEG>(graph, path) * 0.5F;
      metrics_map[path.front()][path.back()].k33 += k33_factor;
      metrics_map[path.back()][path.front()].k33 += k33_factor;

      float k5_factor = calculatePathImpact<K5_MAX_DEG>(graph, path) * 0.5F;
      metrics_map[path.front()][path.back()].k5 += k5_factor;
      metrics_map[path.back()][path.front()].k5 += k5_factor;
      ;

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
  size_t mm_size{metrics_map.size()};

  BAEdge& last_added_edge = graph.getLastAddedEdge();  // source is added earlier than target
  size_t source{last_added_edge.target};
  size_t target{last_added_edge.source};

  std::cout << "New edge: " << source << " -> " << target << std::endl;

  std::vector<BANode> source_neighbours = graph.getNeighbours(source);
  std::vector<BANode> target_neighbours = graph.getNeighbours(target);

  std::cout << "Source(" << source << ")"
            << " neighbours: " << source_neighbours.size() << std::endl;
  std::cout << "Target(" << target << ")"
            << " neighbours: " << target_neighbours.size() << std::endl;
  std::cout << std::endl;

  if (source_neighbours.size() == 1 && target_neighbours.size() == 1) {
    float k33 = psi<K33_MAX_DEG>(1) * psi<K33_MAX_DEG>(1);
    float k5 = psi<K5_MAX_DEG>(1) * psi<K5_MAX_DEG>(1);

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

      float val = metrics_map[node][target].k33 / psi<K33_MAX_DEG>(graph.getDegree(target) - 1) *
                  psi<K33_MAX_DEG>(graph.getDegree(source)) * phi(graph.getDegree(target));
      new_metrics_map[node][source].k33 += val;
      new_metrics_map[source][node].k33 += val;
    }

    // expand metric for path source -> target
    float val = psi<K33_MAX_DEG>(1) * psi<K33_MAX_DEG>(1);
    new_metrics_map[source][target].k33 += val;
    new_metrics_map[target][source].k33 += val;
  }
  else {}

  metrics_map = std::move(new_metrics_map);
}

void AvoidingKuratowskiGraphsPainter::updateShortestPathsMetric(BAGraph& graph,
                                                                AvoidingKuratowskiGraphsPainter::MetricsMap& metrics_map) {
  metrics_map.clear();
  for (size_t v = 0; v < graph.getNodesNumber(); ++v) {
    metrics_map[v] = std::map<size_t, Metric>();
    for (size_t u = 0; u < graph.getNodesNumber(); ++u) {
      metrics_map[v][u] = Metric{};
    }
  }

  for (size_t source = 0; source < graph.getNodesNumber(); ++source) {
    std::vector<size_t> distances(graph.getNodesNumber(), std::numeric_limits<size_t>::max());
    std::vector<std::vector<size_t>> predecessors(graph.getNodesNumber(), std::vector<size_t>());
    std::queue<size_t> queue;

    distances[source] = 0;
    queue.push(source);

    while (!queue.empty()) {
      size_t node = queue.front();
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
        size_t last = path.back();

        if (last == source) {
          float k33_factor = calculatePathImpact<K33_MAX_DEG>(graph, path) * 0.5F;
          metrics_map[path.front()][path.back()].k33 += k33_factor;
          metrics_map[path.back()][path.front()].k33 += k33_factor;

          float k5_factor = calculatePathImpact<K5_MAX_DEG>(graph, path) * 0.5F;
          metrics_map[path.front()][path.back()].k5 += k5_factor;
          metrics_map[path.back()][path.front()].k5 += k5_factor;
        }
        else {
          for (size_t pred : predecessors[last]) {
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
    AvoidingKuratowskiGraphsPainter::MetricsMap& metrics_map) {
  float k33{0.0F};
  float k5{0.0F};

  for (const auto& outerPair : metrics_map) {
    for (const auto& innerPair : outerPair.second) {
      k33 += innerPair.second.k33;
      k5 += innerPair.second.k5;
    }
  }

  return {.k33 = k33, .k5 = k5};
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

float AvoidingKuratowskiGraphsPainter::phi(size_t k) {
  //  if (k < 2) {
  //	return 1;
  //  }
  //
  //  return 2.0F / static_cast<float>(k * (k - 1));

  //    return 1;

  if (k < 2) {
    return 1;
  }

  return 1.0F / static_cast<float>(k - 1);  // TODO: change to 1 / (k choose 2)
}

template <size_t MAX_DEG>
float AvoidingKuratowskiGraphsPainter::psi(size_t k) {
  return static_cast<float>(k) / static_cast<float>(MAX_DEG);
  //  return static_cast<float>(std::min(k, MAX_DEG)) / static_cast<float>(MAX_DEG);
}

}  // namespace graph::random
