#include "avoiding_kuratowski_graphs_painter.hpp"

#include <limits>    // std::numeric_limits
#include <stack>    // std::stack
#include <queue>    // std::queue
#include <vector>   // std::vector

namespace graph::random {
AvoidingKuratowskiGraphsPainter::AvoidingKuratowskiGraphsPainter(std::vector<std::string> edges_colors)
	: GraphPainter({}, std::move(edges_colors)) {
}

void AvoidingKuratowskiGraphsPainter::paintNode(BAGraph &graph, BANode &node) {
  nodes_properties_.emplace_back(NodeProperties{0.0f, 0.0f});
}

void AvoidingKuratowskiGraphsPainter::paintEdge(BAGraph &graph, BAEdge &edge) {
//  update_metric(graph);
//
  // normalizing the metric
//  float sum_metric{0.0f};
//  for (const auto &node_properties : nodes_properties_) {
//	sum_metric += node_properties.k_5_probability;
//  }
//  for(auto &node_properties : nodes_properties_) {
//	node_properties.k_5_probability /= sum_metric;
//  }
//
//  graph.saveToFile("graph1-" + std::to_string(graph.getNodesNumber()) + ".txt");
//
//  std::ofstream file("graph1-metric-" + std::to_string(graph.getNodesNumber()) + ".txt");
//
//  for(size_t i = 0; i < graph.getNodesNumber(); ++i) {
//	file << "Node: " << i << " K_5: " << nodes_properties_[i].k_5_probability << std::endl;
//  }
//
//  std::cout << "Done\n";
}

void AvoidingKuratowskiGraphsPainter::reset() {}

AvoidingKuratowskiGraphsPainter::MetricsMap AvoidingKuratowskiGraphsPainter::calculate_all_paths_metric(BAGraph &graph) {
  MetricsMap metrics_map;
  for (size_t v = 0; v < graph.getNodesNumber(); ++v) {
	metrics_map[v] = std::map<size_t, std::pair<float, float>>();
	for (size_t u = 0; u < graph.getNodesNumber(); ++u) {
	  metrics_map[v][u] = std::make_pair(0.0F, 0.0F);
	}
  }

  auto update_k33_metric = [&](size_t v, size_t u, float m) {
	metrics_map[v][u].first += m;
	metrics_map[u][v].first += m;

  };

  auto update_k5_metric = [&](size_t v, size_t u, float m) {
	metrics_map[v][u].second += m;
	metrics_map[u][v].second += m;
  };

  for (size_t source = 0; source < graph.getNodesNumber(); ++source) {
	std::stack<std::vector<size_t>> paths;
	paths.push({source});

	while (!paths.empty()) {
	  std::vector<size_t> path = paths.top();
	  paths.pop();
	  size_t last = path.back();

	  // Print the current path
//	  for (size_t vertex : path) {
//		std::cout << vertex << " ";
//	  }
//	  std::cout << std::endl;

	  update_k33_metric(path.front(), path.back(), calculate_path_impact(graph, path, 3) * 0.5F);
	  update_k5_metric(path.front(), path.back(), calculate_path_impact(graph, path, 4) * 0.5F);

	  for (const auto &neighbor : graph.getNeighbours(last)) {
		if (std::find(path.begin(), path.end(), neighbor.id) == path.end()) {
		  std::vector<size_t> newPath = path;
		  newPath.push_back(neighbor.id);
		  paths.push(newPath);
		}
	  }
	}
  }

  return metrics_map;
}

AvoidingKuratowskiGraphsPainter::MetricsMap AvoidingKuratowskiGraphsPainter::calculate_shortest_paths_metric(BAGraph &graph) {
  MetricsMap metrics_map;
  for (size_t v = 0; v < graph.getNodesNumber(); ++v) {
	metrics_map[v] = std::map<size_t, std::pair<float, float>>();
	for (size_t u = 0; u < graph.getNodesNumber(); ++u) {
	  metrics_map[v][u] = std::make_pair(0.0F, 0.0F);
	}
  }

  auto update_k33_metric = [&](size_t v, size_t u, float m) {
	metrics_map[v][u].first += m;
	metrics_map[u][v].first += m;

  };

  auto update_k5_metric = [&](size_t v, size_t u, float m) {
	metrics_map[v][u].second += m;
	metrics_map[u][v].second += m;
  };

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
		  update_k33_metric(path.front(), path.back(), calculate_path_impact(graph, path, 3) * 0.5F);
		  update_k5_metric(path.front(), path.back(), calculate_path_impact(graph, path, 4) * 0.5F);
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

  return metrics_map;
}

float AvoidingKuratowskiGraphsPainter::calculate_path_impact(BAGraph &graph,
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
