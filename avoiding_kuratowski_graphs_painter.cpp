#include "avoiding_kuratowski_graphs_painter.hpp"

#include <chrono>
#include <limits>   // std::numeric_limits
#include <stack>    // std::stack
#include <queue>    // std::queue
#include <vector>   // std::vector

namespace graph::random {
AvoidingKuratowskiGraphsPainter::AvoidingKuratowskiGraphsPainter(std::vector<ColorType> edges_colors)
	: GraphPainter({}, std::move(edges_colors)) {

  for (const auto &edges_color : edges_colors_) {
	embeddings_[edges_color] = BAGraph{};
	metrics_map_[edges_color] = MetricsMap{};
  }
}

void AvoidingKuratowskiGraphsPainter::paintNode(BAGraph &graph, BANode &node) {
  for (const auto &edges_color : edges_colors_) {
	embeddings_[edges_color].addNode();

	for (auto &it : metrics_map_[edges_color]) {
	  it.second[node.id] = Metric{};
	}

	metrics_map_[edges_color][node.id] = std::map<size_t, Metric>{};
	for (size_t v = 0; v < node.id; ++v) {
	  metrics_map_[edges_color][node.id][v] = Metric{};
	}
  }
}

void AvoidingKuratowskiGraphsPainter::paintEdge(BAGraph &graph, BAEdge &edge) {
  auto start_coloring = std::chrono::high_resolution_clock::now();

  std::cout << "Painting edge: " << edge.source << " " << edge.target << std::endl;

  std::map<ColorType, float> metrics;
  for (const auto &edges_color : edges_colors_) {
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
  for (const auto &edges_color : edges_colors_) {
	if (max_metric_value > metrics[edges_color]) {
	  max_metric_value = metrics[edges_color];
	}
  }

  std::vector<ColorType> max_colors;
  for (const auto &edges_color : edges_colors_) {
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
}

void AvoidingKuratowskiGraphsPainter::reset() {}

void AvoidingKuratowskiGraphsPainter::updateAllPathsMetric(
	BAGraph &graph, AvoidingKuratowskiGraphsPainter::MetricsMap &metrics_map) {
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

	  float k33_factor = calculatePathImpact(graph, path, 3) * 0.5F;
	  metrics_map[path.front()][path.back()].k33 += k33_factor;
	  metrics_map[path.back()][path.front()].k33 += k33_factor;

	  float k5_factor = calculatePathImpact(graph, path, 4) * 0.5F;
	  metrics_map[path.front()][path.back()].k5 += k5_factor;
	  metrics_map[path.back()][path.front()].k5 += k5_factor;;

	  for (const auto &neighbor : graph.getNeighbours(last)) {
		if (std::find(path.begin(), path.end(), neighbor.id) == path.end()) {
		  std::vector<size_t> newPath = path;
		  newPath.push_back(neighbor.id);
		  paths.push(newPath);
		}
	  }
	}
  }
}

void AvoidingKuratowskiGraphsPainter::updateAllPathsMetricSmart(
	BAGraph &graph, AvoidingKuratowskiGraphsPainter::MetricsMap &metrics_map) {

}

void AvoidingKuratowskiGraphsPainter::updateShortestPathsMetric(
	BAGraph &graph, AvoidingKuratowskiGraphsPainter::MetricsMap &metrics_map) {
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

	  for (const auto &neighbor : graph.getNeighbours(node)) {
		if (distances[neighbor.id] == std::numeric_limits<size_t>::max()) {
		  distances[neighbor.id] = distances[node] + 1;
		  predecessors[neighbor.id].push_back(node);
		  queue.push(neighbor.id);
		} else if (distances[neighbor.id] == distances[node] + 1) {
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
		  float k33_factor = calculatePathImpact(graph, path, 3) * 0.5F;
		  metrics_map[path.front()][path.back()].k33 += k33_factor;
		  metrics_map[path.back()][path.front()].k33 += k33_factor;

		  float k5_factor = calculatePathImpact(graph, path, 4) * 0.5F;
		  metrics_map[path.front()][path.back()].k5 += k5_factor;
		  metrics_map[path.back()][path.front()].k5 += k5_factor;
		} else {
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

void AvoidingKuratowskiGraphsPainter::normalizeMetric(std::vector<Metric> &metrics) {
  float k33_sum{0.0F};
  float k5_sum{0.0F};

  for (const auto &metric : metrics) {
	k33_sum += metric.k33;
	k5_sum += metric.k5;
  }

  for (auto &metric : metrics) {
	metric.k33 /= k33_sum;
	metric.k5 /= k5_sum;
  }
}

AvoidingKuratowskiGraphsPainter::Metric
AvoidingKuratowskiGraphsPainter::sumMetric(AvoidingKuratowskiGraphsPainter::MetricsMap &metrics_map) {
  float k33{0.0F};
  float k5{0.0F};

  for (const auto &outerPair : metrics_map) {
	for (const auto &innerPair : outerPair.second) {
	  k33 += innerPair.second.k33;
	  k5 += innerPair.second.k5;
	}
  }

  return {.k33 = k33, .k5 = k5};
}

float AvoidingKuratowskiGraphsPainter::calculatePathImpact(BAGraph &graph,
															 std::vector<size_t> &path,
															 int max_degree) {
  if (path.size() < 2) {
	return 0;
  }

  size_t in_deg{graph.getDegree(path.front())};
  size_t out_deg{graph.getDegree(path.back())};

  float result{1.0F};
  result *= static_cast<float>(std::min<size_t>(max_degree, in_deg)) / static_cast<float>(max_degree);
  result *= static_cast<float>(std::min<size_t>(max_degree, out_deg)) / static_cast<float>(max_degree);

  for (size_t i = 1; i < (path.size() - 1); ++i) {
	float factor{1.0F / static_cast<float>(graph.getDegree(path.at(i)) - 1)};
	result *= factor;
  }

  return result;
}

} // namespace graph::random
