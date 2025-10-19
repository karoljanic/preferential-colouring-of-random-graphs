#include <fstream>   // std::ofstream
#include <iostream>  // std::cerr, std::cout

#include "../include/max_planar_subgraph.hpp"
#include "../include/preferential_coloring.hpp"
#include "../include/random_graph_factory.hpp"

struct TestResult {
  double crossingEdges{0};
  double time{0.0};
};

struct TestResults {
  std::vector<size_t> crossingEdges;
  std::vector<double> times;
};

struct AveragedTestResults {
  double crossingEdges{0.0};
  double time{0.0};
};

template <typename T>
double mean(const std::vector<T>& values) {
  return std::accumulate(values.begin(), values.end(), 0.0) / values.size();
}

template <typename T>
double stddev(const std::vector<T>& values, double mean) {
  double sum = 0;
  for (auto val : values) {
    sum += (val - mean) * (val - mean);
  }

  if (values.size() <= 1) {
    return 0.0;
  }

  return std::sqrt(sum / (values.size() - 1));
}

void appendTestResult(TestResults& results, const TestResult& result) {
  results.crossingEdges.push_back(result.crossingEdges);
  results.times.push_back(result.time);
}

void summarizeTestResults(const TestResults& results, AveragedTestResults& avg, AveragedTestResults& stddevRes) {
  avg.crossingEdges = mean(results.crossingEdges);
  avg.time = mean(results.times);

  stddevRes.crossingEdges = stddev(results.crossingEdges, avg.crossingEdges);
  stddevRes.time = stddev(results.times, avg.time);
}

graph::random::BAGraph generateGraph(int method, size_t initial_vertices, size_t edges_per_vertex, size_t final_vertices) {
  graph::random::RandomGraphFactory random_graph_factory;
  switch (method) {
    case 1:
      return random_graph_factory.createBarabasiAlbertWithPreferentialAttachmentRepeatedNodes(initial_vertices, final_vertices,
                                                                                              edges_per_vertex);
    case 2:
      return random_graph_factory.createBarabasiAlbertWithPreferentialAttachmentBatageljBrandes(initial_vertices, final_vertices,
                                                                                                edges_per_vertex);
    case 3:
      return random_graph_factory.createBarabasiAlbertWithLinkSelection(initial_vertices, final_vertices, edges_per_vertex);
    case 4:
      return random_graph_factory.createBarabasiAlbertWithCopyingModel(initial_vertices, final_vertices, edges_per_vertex, 0.5);
    case 5:
      return random_graph_factory.createBarabasiAlbertWithLCDModel(final_vertices, edges_per_vertex);
    default:
      throw std::invalid_argument("Invalid method. Choose between 1 and 5.");
  }
}

TestResult runTestAllAlgorithms(const graph::random::BAGraph& graph, const size_t colors_number) {
  TestResult result;

  auto start = std::chrono::high_resolution_clock::now();
  const std::vector<graph::random::BAGraph> subgraphs = graph::random::PreferentialColoring::color(graph, colors_number);
  auto end = std::chrono::high_resolution_clock::now();

  size_t crossing_edges_sum = 0;
  for (const auto& subgraph : subgraphs) {
    crossing_edges_sum += graph::random::MaxPlanarSubgraph::crossingEdges(subgraph);
  }
  result.crossingEdges = static_cast<double>(crossing_edges_sum) / static_cast<double>(subgraphs.size());

  result.time = std::chrono::duration<double, std::milli>(end - start).count();

  return result;
}

int main(int argc, char* argv[]) {
  if (argc != 7) {
    std::cerr << "Usage: " << argv[0]
              << " <method> <initial_vertices_number> <edges_per_vertex> <final_vertices_number> <colors_number> <repetitions>\n";
    return 0;
  }

  const int method = std::stoi(argv[1]);
  const size_t initial_vertices_number = static_cast<size_t>(std::stoul(argv[2]));
  const size_t edges_per_vertex = static_cast<size_t>(std::stoul(argv[3]));
  const size_t final_vertices_number = static_cast<size_t>(std::stoul(argv[4]));
  const size_t colors_number = static_cast<size_t>(std::stoul(argv[5]));
  const size_t repetitions = static_cast<size_t>(std::stoul(argv[6]));

  TestResults totalResult;

  for (size_t rep = 0; rep < repetitions; ++rep) {
    graph::random::BAGraph graph = generateGraph(method, initial_vertices_number, edges_per_vertex, final_vertices_number);

    const TestResult result = runTestAllAlgorithms(graph, colors_number);
    appendTestResult(totalResult, result);
  }

  AveragedTestResults averageResult, stddevResult;
  summarizeTestResults(totalResult, averageResult, stddevResult);

  std::cout << averageResult.crossingEdges << " " << stddevResult.crossingEdges << " " << averageResult.time << " "
            << stddevResult.time << std::endl;

  return 0;
}