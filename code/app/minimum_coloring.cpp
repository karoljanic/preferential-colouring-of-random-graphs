#include <iostream>  // std::cerr, std::cout

#include "../include/max_planar_subgraph.hpp"
#include "../include/random_graph_factory.hpp"

struct TestResult {
  size_t boyerMyrvoldColors{0}, mstColors{0}, cactusColors{0};
  double boyerMyrvoldTime{0.0}, mstTime{0.0}, cactusTime{0.0};
};

struct TestResults {
  std::vector<size_t> boyerMyrvoldColors{}, mstColors{}, cactusColors{};
  std::vector<double> boyerMyrvoldTime{}, mstTime{}, cactusTime{};
};

struct AveragedTestResults {
  double boyerMyrvoldColors{0.0}, mstColors{0.0}, cactusColors{0.0};
  double boyerMyrvoldTime{0.0}, mstTime{0.0}, cactusTime{0.0};
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
  results.boyerMyrvoldColors.push_back(result.boyerMyrvoldColors);
  results.mstColors.push_back(result.mstColors);
  results.cactusColors.push_back(result.cactusColors);

  results.boyerMyrvoldTime.push_back(result.boyerMyrvoldTime);
  results.mstTime.push_back(result.mstTime);
  results.cactusTime.push_back(result.cactusTime);
}

void summarizeTestResults(const TestResults& results, AveragedTestResults& avg, AveragedTestResults& stddevRes) {
  avg.boyerMyrvoldColors = mean(results.boyerMyrvoldColors);
  avg.mstColors = mean(results.mstColors);
  avg.cactusColors = mean(results.cactusColors);

  avg.boyerMyrvoldTime = mean(results.boyerMyrvoldTime);
  avg.mstTime = mean(results.mstTime);
  avg.cactusTime = mean(results.cactusTime);

  stddevRes.boyerMyrvoldColors = stddev(results.boyerMyrvoldColors, avg.boyerMyrvoldColors);
  stddevRes.mstColors = stddev(results.mstColors, avg.mstColors);
  stddevRes.cactusColors = stddev(results.cactusColors, avg.cactusColors);

  stddevRes.boyerMyrvoldTime = stddev(results.boyerMyrvoldTime, avg.boyerMyrvoldTime);
  stddevRes.mstTime = stddev(results.mstTime, avg.mstTime);
  stddevRes.cactusTime = stddev(results.cactusTime, avg.cactusTime);
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

TestResult runTestAllAlgorithms(const graph::random::BAGraph& graph) {
  TestResult result;

  {
    graph::random::BAGraph graph_copy = graph;
    const auto start = std::chrono::high_resolution_clock::now();
    while (graph_copy.getEdgesNumber() > 0) {
      graph::random::BAGraph max_planar_subgraph;
      graph::random::PlanarityTest::boyerMyrvoldPlanarSubgraph(graph_copy, max_planar_subgraph);
      graph::random::MaxPlanarSubgraph::maximizeSubgraph(graph_copy, max_planar_subgraph);
      result.boyerMyrvoldColors += 1.0;

      for (const auto& edge : max_planar_subgraph.getEdges()) {
        graph_copy.removeEdge(edge.source, edge.target);
      }
    }
    const auto end = std::chrono::high_resolution_clock::now();

    result.boyerMyrvoldTime = std::chrono::duration<double>(end - start).count();
  }

  {
    graph::random::BAGraph graph_copy = graph;
    const auto start = std::chrono::high_resolution_clock::now();
    while (graph_copy.getEdgesNumber() > 0) {
      graph::random::BAGraph max_planar_subgraph;
      graph::random::MaxPlanarSubgraph::mstBased(graph_copy, max_planar_subgraph);
      graph::random::MaxPlanarSubgraph::maximizeSubgraph(graph_copy, max_planar_subgraph);
      result.mstColors += 1.0;

      for (const auto& edge : max_planar_subgraph.getEdges()) {
        graph_copy.removeEdge(edge.source, edge.target);
      }
    }
    const auto end = std::chrono::high_resolution_clock::now();

    result.mstTime = std::chrono::duration<double>(end - start).count();
  }

  {
    graph::random::BAGraph graph_copy = graph;
    const auto start = std::chrono::high_resolution_clock::now();
    while (graph_copy.getEdgesNumber() > 0) {
      graph::random::BAGraph max_planar_subgraph;
      graph::random::MaxPlanarSubgraph::cactusBased(graph_copy, max_planar_subgraph);
      graph::random::MaxPlanarSubgraph::maximizeSubgraph(graph_copy, max_planar_subgraph);
      result.cactusColors += 1.0;

      for (const auto& edge : max_planar_subgraph.getEdges()) {
        graph_copy.removeEdge(edge.source, edge.target);
      }
    }
    const auto end = std::chrono::high_resolution_clock::now();

    result.cactusTime = std::chrono::duration<double>(end - start).count();
  }

  return result;
}

TestResult runTestAllWeightedAlgorithms(const graph::random::BAGraph& graph) {
  TestResult result;

  {
    graph::random::BAGraph graph_copy = graph;
    const auto start = std::chrono::high_resolution_clock::now();
    while (graph_copy.getEdgesNumber() > 0) {
      graph::random::BAGraph max_planar_subgraph;
      graph::random::PlanarityTest::boyerMyrvoldPlanarSubgraph(graph_copy, max_planar_subgraph);
      graph::random::MaxPlanarSubgraph::weightedMaximizeSubgraph(graph_copy, max_planar_subgraph);
      result.boyerMyrvoldColors += 1.0;

      for (const auto& edge : max_planar_subgraph.getEdges()) {
        graph_copy.removeEdge(edge.source, edge.target);
      }
    }
    const auto end = std::chrono::high_resolution_clock::now();

    result.boyerMyrvoldTime = std::chrono::duration<double>(end - start).count();
  }

  {
    graph::random::BAGraph graph_copy = graph;
    const auto start = std::chrono::high_resolution_clock::now();
    while (graph_copy.getEdgesNumber() > 0) {
      graph::random::BAGraph max_planar_subgraph;
      graph::random::MaxPlanarSubgraph::weightedMstBased(graph_copy, max_planar_subgraph);
      graph::random::MaxPlanarSubgraph::weightedMaximizeSubgraph(graph_copy, max_planar_subgraph);
      result.mstColors += 1.0;

      for (const auto& edge : max_planar_subgraph.getEdges()) {
        graph_copy.removeEdge(edge.source, edge.target);
      }
    }
    const auto end = std::chrono::high_resolution_clock::now();

    result.mstTime = std::chrono::duration<double>(end - start).count();
  }

  {
    graph::random::BAGraph graph_copy = graph;
    const auto start = std::chrono::high_resolution_clock::now();
    while (graph_copy.getEdgesNumber() > 0) {
      graph::random::BAGraph max_planar_subgraph;
      graph::random::MaxPlanarSubgraph::weightedCactusBased(graph_copy, max_planar_subgraph);
      graph::random::MaxPlanarSubgraph::weightedMaximizeSubgraph(graph_copy, max_planar_subgraph);
      result.cactusColors += 1.0;

      for (const auto& edge : max_planar_subgraph.getEdges()) {
        graph_copy.removeEdge(edge.source, edge.target);
      }
    }
    const auto end = std::chrono::high_resolution_clock::now();

    result.cactusTime = std::chrono::duration<double>(end - start).count();
  }

  return result;
}

int main(int argc, char* argv[]) {
  if (argc != 6) {
    std::cerr << "Usage: " << argv[0]
              << " <method> <initial_vertices_number> <edges_per_vertex> <final_vertices_number> <repetitions>\n";
    return 0;
  }

  const int method = std::stoi(argv[1]);
  const size_t initial_vertices_number = static_cast<size_t>(std::stoul(argv[2]));
  const size_t edges_per_vertex = static_cast<size_t>(std::stoul(argv[3]));
  const size_t final_vertices_number = static_cast<size_t>(std::stoul(argv[4]));
  const size_t repetitions = static_cast<size_t>(std::stoul(argv[5]));

  TestResults basicTotal, weightedTotal;

  for (size_t rep = 0; rep < repetitions; ++rep) {
    const graph::random::BAGraph graph = generateGraph(method, initial_vertices_number, edges_per_vertex, final_vertices_number);

    const TestResult basicResult = runTestAllAlgorithms(graph);
    appendTestResult(basicTotal, basicResult);

    const TestResult weightedResult = runTestAllWeightedAlgorithms(graph);
    appendTestResult(weightedTotal, weightedResult);
  }

  AveragedTestResults basicAverage, basicStddev;
  AveragedTestResults weightedAverage, weightedStddev;

  summarizeTestResults(basicTotal, basicAverage, basicStddev);
  summarizeTestResults(weightedTotal, weightedAverage, weightedStddev);

  std::cout << basicAverage.boyerMyrvoldColors << " " << basicStddev.boyerMyrvoldColors << " "
            << weightedAverage.boyerMyrvoldColors << " " << weightedStddev.boyerMyrvoldColors << " " << basicAverage.mstColors
            << " " << basicStddev.mstColors << " " << weightedAverage.mstColors << " " << weightedStddev.mstColors << " "
            << basicAverage.cactusColors << " " << basicStddev.cactusColors << " " << weightedAverage.cactusColors << " "
            << weightedStddev.cactusColors << std::endl;

  std::cout << basicAverage.boyerMyrvoldTime << " " << basicStddev.boyerMyrvoldTime << " " << weightedAverage.boyerMyrvoldTime
            << " " << weightedStddev.boyerMyrvoldTime << " " << basicAverage.mstTime << " " << basicStddev.mstTime << " "
            << weightedAverage.mstTime << " " << weightedStddev.mstTime << " " << basicAverage.cactusTime << " "
            << basicStddev.cactusTime << " " << weightedAverage.cactusTime << " " << weightedStddev.cactusTime << std::endl;

  return 0;
}