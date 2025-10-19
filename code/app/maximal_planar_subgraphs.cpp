#include <iostream>  // std::cerr, std::cout

#include "../include/max_planar_subgraph.hpp"
#include "../include/random_graph_factory.hpp"

struct TestResult {
  size_t mstEdges{0}, cactusEdges{0}, boyerMyrvoldEdges{0};
  double mstRatio{0.0}, cactusRatio{0.0}, boyerMyrvoldRatio{0.0};
  double mstTime{0.0}, cactusTime{0.0}, boyerMyrvoldTime{0.0};
};

struct TestResults {
  std::vector<size_t> mstEdges{}, cactusEdges{}, boyerMyrvoldEdges{};
  std::vector<double> mstRatio{}, cactusRatio{}, boyerMyrvoldRatio{};
  std::vector<double> mstTime{}, cactusTime{}, boyerMyrvoldTime{};
};

struct AveragedTestResults {
  double mstEdges{0}, cactusEdges{0}, boyerMyrvoldEdges{0};
  double mstRatio{0}, cactusRatio{0}, boyerMyrvoldRatio{0};
  double mstTime{0}, cactusTime{0}, boyerMyrvoldTime{0};
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
  results.boyerMyrvoldEdges.push_back(result.boyerMyrvoldEdges);
  results.mstEdges.push_back(result.mstEdges);
  results.cactusEdges.push_back(result.cactusEdges);

  results.boyerMyrvoldRatio.push_back(result.boyerMyrvoldRatio);
  results.mstRatio.push_back(result.mstRatio);
  results.cactusRatio.push_back(result.cactusRatio);

  results.boyerMyrvoldTime.push_back(result.boyerMyrvoldTime);
  results.mstTime.push_back(result.mstTime);
  results.cactusTime.push_back(result.cactusTime);
}

void summarizeTestResults(const TestResults& results, AveragedTestResults& avg, AveragedTestResults& stddevRes) {
  avg.mstEdges = mean(results.mstEdges);
  avg.cactusEdges = mean(results.cactusEdges);
  avg.boyerMyrvoldEdges = mean(results.boyerMyrvoldEdges);

  avg.mstRatio = mean(results.mstRatio);
  avg.cactusRatio = mean(results.cactusRatio);
  avg.boyerMyrvoldRatio = mean(results.boyerMyrvoldRatio);

  avg.mstTime = mean(results.mstTime);
  avg.cactusTime = mean(results.cactusTime);
  avg.boyerMyrvoldTime = mean(results.boyerMyrvoldTime);

  stddevRes.mstEdges = stddev(results.mstEdges, avg.mstEdges);
  stddevRes.cactusEdges = stddev(results.cactusEdges, avg.cactusEdges);
  stddevRes.boyerMyrvoldEdges = stddev(results.boyerMyrvoldEdges, avg.boyerMyrvoldEdges);

  stddevRes.mstRatio = stddev(results.mstRatio, avg.mstRatio);
  stddevRes.cactusRatio = stddev(results.cactusRatio, avg.cactusRatio);
  stddevRes.boyerMyrvoldRatio = stddev(results.boyerMyrvoldRatio, avg.boyerMyrvoldRatio);

  stddevRes.mstTime = stddev(results.mstTime, avg.mstTime);
  stddevRes.cactusTime = stddev(results.cactusTime, avg.cactusTime);
  stddevRes.boyerMyrvoldTime = stddev(results.boyerMyrvoldTime, avg.boyerMyrvoldTime);
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

TestResult runTestAllAlgorithms(const graph::random::BAGraph& graph, bool maximize) {
  TestResult result;

  {
    graph::random::BAGraph subgraph;
    const auto start = std::chrono::high_resolution_clock::now();
    graph::random::MaxPlanarSubgraph::mstBased(graph, subgraph);
    if (maximize) {
      graph::random::MaxPlanarSubgraph::maximizeSubgraph(graph, subgraph);
    }
    const auto end = std::chrono::high_resolution_clock::now();
    result.mstEdges = subgraph.getEdgesNumber();
    result.mstRatio = static_cast<double>(result.mstEdges) / static_cast<double>(3 * graph.getNodesNumber() - 6);
    result.mstTime = std::chrono::duration<double, std::milli>(end - start).count();
  }

  {
    graph::random::BAGraph subgraph;
    const auto start = std::chrono::high_resolution_clock::now();
    graph::random::MaxPlanarSubgraph::cactusBased(graph, subgraph);
    if (maximize) {
      graph::random::MaxPlanarSubgraph::maximizeSubgraph(graph, subgraph);
    }
    const auto end = std::chrono::high_resolution_clock::now();
    result.cactusEdges = subgraph.getEdgesNumber();
    result.cactusRatio = static_cast<double>(result.cactusEdges) / static_cast<double>(3 * graph.getNodesNumber() - 6);
    result.cactusTime = std::chrono::duration<double, std::milli>(end - start).count();
  }

  {
    graph::random::BAGraph subgraph;
    const auto start = std::chrono::high_resolution_clock::now();
    graph::random::PlanarityTest::boyerMyrvoldPlanarSubgraph(graph, subgraph);
    if (maximize) {
      graph::random::MaxPlanarSubgraph::maximizeSubgraph(graph, subgraph);
    }
    const auto end = std::chrono::high_resolution_clock::now();
    result.boyerMyrvoldEdges = subgraph.getEdgesNumber();
    result.boyerMyrvoldRatio =
        static_cast<double>(result.boyerMyrvoldEdges) / static_cast<double>(3 * graph.getNodesNumber() - 6);
    result.boyerMyrvoldTime = std::chrono::duration<double, std::milli>(end - start).count();
  }

  return result;
}

TestResult runTestAllWeightedAlgorithms(const graph::random::BAGraph& graph, bool maximize) {
  TestResult result;

  {
    graph::random::BAGraph subgraph;
    const auto start = std::chrono::high_resolution_clock::now();
    graph::random::MaxPlanarSubgraph::weightedMstBased(graph, subgraph);
    if (maximize) {
      graph::random::MaxPlanarSubgraph::weightedMaximizeSubgraph(graph, subgraph);
    }
    const auto end = std::chrono::high_resolution_clock::now();
    result.mstEdges = subgraph.getEdgesNumber();
    result.mstRatio = static_cast<double>(result.mstEdges) / static_cast<double>(3 * graph.getNodesNumber() - 6);
    result.mstTime = std::chrono::duration<double, std::milli>(end - start).count();
  }

  {
    graph::random::BAGraph subgraph;
    const auto start = std::chrono::high_resolution_clock::now();
    graph::random::MaxPlanarSubgraph::weightedCactusBased(graph, subgraph);
    if (maximize) {
      graph::random::MaxPlanarSubgraph::weightedMaximizeSubgraph(graph, subgraph);
    }
    const auto end = std::chrono::high_resolution_clock::now();
    result.cactusEdges = subgraph.getEdgesNumber();
    result.cactusRatio = static_cast<double>(result.cactusEdges) / static_cast<double>(3 * graph.getNodesNumber() - 6);
    result.cactusTime = std::chrono::duration<double, std::milli>(end - start).count();
  }

  {
    graph::random::BAGraph subgraph;
    const auto start = std::chrono::high_resolution_clock::now();
    graph::random::PlanarityTest::boyerMyrvoldPlanarSubgraph(graph, subgraph);
    if (maximize) {
      graph::random::MaxPlanarSubgraph::weightedMaximizeSubgraph(graph, subgraph);
    }
    const auto end = std::chrono::high_resolution_clock::now();
    result.boyerMyrvoldEdges = subgraph.getEdgesNumber();
    result.boyerMyrvoldRatio =
        static_cast<double>(result.boyerMyrvoldEdges) / static_cast<double>(3 * graph.getNodesNumber() - 6);
    result.boyerMyrvoldTime = std::chrono::duration<double, std::milli>(end - start).count();
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

  TestResults basicTotal, weightedTotal, maximizedBasicTotal, maximizedWeightedTotal;

  for (size_t rep = 0; rep < repetitions; ++rep) {
    const graph::random::BAGraph graph = generateGraph(method, initial_vertices_number, edges_per_vertex, final_vertices_number);

    const TestResult basicResult = runTestAllAlgorithms(graph, false);
    appendTestResult(basicTotal, basicResult);

    const TestResult weightedResult = runTestAllWeightedAlgorithms(graph, false);
    appendTestResult(weightedTotal, weightedResult);

    const TestResult maximizedBasicResult = runTestAllAlgorithms(graph, true);
    appendTestResult(maximizedBasicTotal, maximizedBasicResult);

    const TestResult maximizedWeightedResult = runTestAllWeightedAlgorithms(graph, true);
    appendTestResult(maximizedWeightedTotal, maximizedWeightedResult);
  }

  AveragedTestResults basicAvg, basicStd, weightedAvg, weightedStd;
  AveragedTestResults maxBasicAvg, maxBasicStd, maxWeightedAvg, maxWeightedStd;

  summarizeTestResults(basicTotal, basicAvg, basicStd);
  summarizeTestResults(weightedTotal, weightedAvg, weightedStd);
  summarizeTestResults(maximizedBasicTotal, maxBasicAvg, maxBasicStd);
  summarizeTestResults(maximizedWeightedTotal, maxWeightedAvg, maxWeightedStd);

  std::cout << basicAvg.boyerMyrvoldEdges << " " << basicStd.boyerMyrvoldEdges << " " << weightedAvg.boyerMyrvoldEdges << " "
            << weightedStd.boyerMyrvoldEdges << " " << maxBasicAvg.boyerMyrvoldEdges << " " << maxBasicStd.boyerMyrvoldEdges
            << " " << maxWeightedAvg.boyerMyrvoldEdges << " " << maxWeightedStd.boyerMyrvoldEdges << " " << basicAvg.mstEdges
            << " " << basicStd.mstEdges << " " << weightedAvg.mstEdges << " " << weightedStd.mstEdges << " "
            << maxBasicAvg.mstEdges << " " << maxBasicStd.mstEdges << " " << maxWeightedAvg.mstEdges << " "
            << maxWeightedStd.mstEdges << " " << basicAvg.cactusEdges << " " << basicStd.cactusEdges << " "
            << weightedAvg.cactusEdges << " " << weightedStd.cactusEdges << " " << maxBasicAvg.cactusEdges << " "
            << maxBasicStd.cactusEdges << " " << maxWeightedAvg.cactusEdges << " " << maxWeightedStd.cactusEdges << std::endl;

  std::cout << basicAvg.boyerMyrvoldRatio << " " << basicStd.boyerMyrvoldRatio << " " << weightedAvg.boyerMyrvoldRatio << " "
            << weightedStd.boyerMyrvoldRatio << " " << maxBasicAvg.boyerMyrvoldRatio << " " << maxBasicStd.boyerMyrvoldRatio
            << " " << maxWeightedAvg.boyerMyrvoldRatio << " " << maxWeightedStd.boyerMyrvoldRatio << " " << basicAvg.mstRatio
            << " " << basicStd.mstRatio << " " << weightedAvg.mstRatio << " " << weightedStd.mstRatio << " "
            << maxBasicAvg.mstRatio << " " << maxBasicStd.mstRatio << " " << maxWeightedAvg.mstRatio << " "
            << maxWeightedStd.mstRatio << " " << basicAvg.cactusRatio << " " << basicStd.cactusRatio << " "
            << weightedAvg.cactusRatio << " " << weightedStd.cactusRatio << " " << maxBasicAvg.cactusRatio << " "
            << maxBasicStd.cactusRatio << " " << maxWeightedAvg.cactusRatio << " " << maxWeightedStd.cactusRatio << std::endl;

  std::cout << basicAvg.boyerMyrvoldTime << " " << basicStd.boyerMyrvoldTime << " " << weightedAvg.boyerMyrvoldTime << " "
            << weightedStd.boyerMyrvoldTime << " " << maxBasicAvg.boyerMyrvoldTime << " " << maxBasicStd.boyerMyrvoldTime << " "
            << maxWeightedAvg.boyerMyrvoldTime << " " << maxWeightedStd.boyerMyrvoldTime << " " << basicAvg.mstTime << " "
            << basicStd.mstTime << " " << weightedAvg.mstTime << " " << weightedStd.mstTime << " " << maxBasicAvg.mstTime << " "
            << maxBasicStd.mstTime << " " << maxWeightedAvg.mstTime << " " << maxWeightedStd.mstTime << " " << basicAvg.cactusTime
            << " " << basicStd.cactusTime << " " << weightedAvg.cactusTime << " " << weightedStd.cactusTime << " "
            << maxBasicAvg.cactusTime << " " << maxBasicStd.cactusTime << " " << maxWeightedAvg.cactusTime << " "
            << maxWeightedStd.cactusTime << std::endl;

  return 0;
}