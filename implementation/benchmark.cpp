#include <chrono>
#include <iostream>

#include "avoiding_kuratowski_graphs_painter.hpp"  // graph::random::AvoidingKuratowskiGraphsPainter
#include "ba_graph.hpp"                            // graph::random::BAGraph
#include "max_planar_subgraph.hpp"                 // graph::random::MaxPlanarSubgraph
#include "metric_calculator.hpp"                   // graph::random::MetricCalculator
#include "preferential_coloring.hpp"               // graph::random::PreferentialColoring
#include "random_graph_factory.hpp"                // graph::random::RandomGraphFactory

using namespace graph::random;

constexpr ColorType Blue{"#1f77b4"};
constexpr ColorType Orange{"#ff7f0e"};
constexpr ColorType Green{"#2ca02c"};
constexpr ColorType Red{"#d62728"};
constexpr ColorType Purple{"#9467bd"};
constexpr ColorType Brown{"#8c564b"};
constexpr ColorType Pink{"#e377c2"};
constexpr ColorType Gray{"#7f7f7f"};
constexpr ColorType Olive{"#bcbd22"};
constexpr ColorType Cyan{"#17becf"};

const std::vector<ColorType> allColors{Blue, Orange, Green, Red, Purple, Brown, Pink, Gray, Olive, Cyan};

struct TestResult {
  double boyerMyrvoldEdges{0};
  double mstEdges{0};
  double cactusEdges{0};

  double boyerMyrvoldRatio{0.0};
  double mstRatio{0.0};
  double cactusRatio{0.0};

  double boyerMyrvoldTime{0.0};
  double mstTime{0.0};
  double cactusTime{0.0};
};

TestResult operator+(const TestResult& lhs, const TestResult& rhs) {
  return TestResult{lhs.boyerMyrvoldEdges + rhs.boyerMyrvoldEdges, lhs.mstEdges + rhs.mstEdges, lhs.cactusEdges + rhs.cactusEdges,
                    lhs.boyerMyrvoldRatio + rhs.boyerMyrvoldRatio, lhs.mstRatio + rhs.mstRatio, lhs.cactusRatio + rhs.cactusRatio,
                    lhs.boyerMyrvoldTime + rhs.boyerMyrvoldTime,   lhs.mstTime + rhs.mstTime,   lhs.cactusTime + rhs.cactusTime};
}

TestResult operator/(const TestResult& lhs, double rhs) {
  return TestResult{lhs.boyerMyrvoldEdges / rhs, lhs.mstEdges / rhs, lhs.cactusEdges / rhs,
                    lhs.boyerMyrvoldRatio / rhs, lhs.mstRatio / rhs, lhs.cactusRatio / rhs,
                    lhs.boyerMyrvoldTime / rhs,  lhs.mstTime / rhs,  lhs.cactusTime / rhs};
}

TestResult testBasicMethods(const BAGraph& graph, bool maximize) {
  TestResult result;

  auto start = std::chrono::high_resolution_clock::now();
  BAGraph max_planar_subgraph_mst;
  MaxPlanarSubgraph::mstBased(graph, max_planar_subgraph_mst);
  if (maximize) {
    MaxPlanarSubgraph::maximizeSubgraph(graph, max_planar_subgraph_mst);
  }
  auto end = std::chrono::high_resolution_clock::now();

  result.mstEdges = max_planar_subgraph_mst.getEdgesNumber();
  result.mstRatio = static_cast<double>(result.mstEdges) / static_cast<double>(3 * graph.getNodesNumber() - 6);
  result.mstTime = std::chrono::duration<double, std::milli>(end - start).count();

  start = std::chrono::high_resolution_clock::now();
  BAGraph max_planar_subgraph_cactus;
  MaxPlanarSubgraph::cactusBased(graph, max_planar_subgraph_cactus);
  if (maximize) {
    MaxPlanarSubgraph::maximizeSubgraph(graph, max_planar_subgraph_cactus);
  }
  end = std::chrono::high_resolution_clock::now();

  result.cactusEdges = max_planar_subgraph_cactus.getEdgesNumber();
  result.cactusRatio = static_cast<double>(result.cactusEdges) / static_cast<double>(3 * graph.getNodesNumber() - 6);
  result.cactusTime = std::chrono::duration<double, std::milli>(end - start).count();

  start = std::chrono::high_resolution_clock::now();
  BAGraph max_planar_subgraph_boyerMyrvold;
  PlanarityTest::boyerMyrvoldPlanarSubgraph(graph, max_planar_subgraph_boyerMyrvold);
  if (maximize) {
    MaxPlanarSubgraph::maximizeSubgraph(graph, max_planar_subgraph_boyerMyrvold);
  }
  end = std::chrono::high_resolution_clock::now();

  result.boyerMyrvoldEdges = max_planar_subgraph_boyerMyrvold.getEdgesNumber();
  result.boyerMyrvoldRatio = static_cast<double>(result.boyerMyrvoldEdges) / static_cast<double>(3 * graph.getNodesNumber() - 6);
  result.boyerMyrvoldTime = std::chrono::duration<double, std::milli>(end - start).count();

  return result;
}

TestResult testWeightedMethods(const BAGraph& graph, bool maximize) {
  TestResult result;

  auto start = std::chrono::high_resolution_clock::now();
  BAGraph max_planar_subgraph_mst;
  MaxPlanarSubgraph::weightedMstBased(graph, max_planar_subgraph_mst);
  if (maximize) {
    MaxPlanarSubgraph::weightedMaximizeSubgraph(graph, max_planar_subgraph_mst);
  }
  auto end = std::chrono::high_resolution_clock::now();

  result.mstEdges = max_planar_subgraph_mst.getEdgesNumber();
  result.mstRatio = static_cast<double>(result.mstEdges) / static_cast<double>(3 * graph.getNodesNumber() - 6);
  result.mstTime = std::chrono::duration<double, std::milli>(end - start).count();

  start = std::chrono::high_resolution_clock::now();
  BAGraph max_planar_subgraph_cactus;
  MaxPlanarSubgraph::weightedCactusBased(graph, max_planar_subgraph_cactus);
  if (maximize) {
    MaxPlanarSubgraph::weightedMaximizeSubgraph(graph, max_planar_subgraph_cactus);
  }
  end = std::chrono::high_resolution_clock::now();

  result.cactusEdges = max_planar_subgraph_cactus.getEdgesNumber();
  result.cactusRatio = static_cast<double>(result.cactusEdges) / static_cast<double>(3 * graph.getNodesNumber() - 6);
  result.cactusTime = std::chrono::duration<double, std::milli>(end - start).count();

  start = std::chrono::high_resolution_clock::now();
  BAGraph max_planar_subgraph_boyerMyrvold;
  PlanarityTest::boyerMyrvoldPlanarSubgraph(graph, max_planar_subgraph_boyerMyrvold);
  if (maximize) {
    MaxPlanarSubgraph::weightedMaximizeSubgraph(graph, max_planar_subgraph_boyerMyrvold);
  }
  end = std::chrono::high_resolution_clock::now();

  result.boyerMyrvoldEdges = max_planar_subgraph_boyerMyrvold.getEdgesNumber();
  result.boyerMyrvoldRatio = static_cast<double>(result.boyerMyrvoldEdges) / static_cast<double>(3 * graph.getNodesNumber() - 6);
  result.boyerMyrvoldTime = std::chrono::duration<double, std::milli>(end - start).count();

  return result;
}

TestResult testColoring(const BAGraph& graph) {
  TestResult result;

  BAGraph graph1 = graph;
  while (graph1.getEdgesNumber() > 0) {
    BAGraph max_planar_subgraph_mst;
    MaxPlanarSubgraph::mstBased(graph1, max_planar_subgraph_mst);
    MaxPlanarSubgraph::maximizeSubgraph(graph1, max_planar_subgraph_mst);
    result.mstRatio += 1.0;

    for (const auto& edge : max_planar_subgraph_mst.getEdges()) {
      graph1.removeEdge(edge.source, edge.target);
    }
  }

  BAGraph graph2 = graph;
  while (graph2.getEdgesNumber() > 0) {
    BAGraph max_planar_subgraph_cactus;
    MaxPlanarSubgraph::cactusBased(graph2, max_planar_subgraph_cactus);
    MaxPlanarSubgraph::maximizeSubgraph(graph2, max_planar_subgraph_cactus);
    result.cactusRatio += 1.0;

    for (const auto& edge : max_planar_subgraph_cactus.getEdges()) {
      graph2.removeEdge(edge.source, edge.target);
    }
  }

  BAGraph graph3 = graph;
  while (graph3.getEdgesNumber() > 0) {
    BAGraph max_planar_subgraph_boyerMyrvold;
    PlanarityTest::boyerMyrvoldPlanarSubgraph(graph3, max_planar_subgraph_boyerMyrvold);
    MaxPlanarSubgraph::maximizeSubgraph(graph3, max_planar_subgraph_boyerMyrvold);
    result.boyerMyrvoldRatio += 1.0;

    for (const auto& edge : max_planar_subgraph_boyerMyrvold.getEdges()) {
      graph3.removeEdge(edge.source, edge.target);
    }
  }

  return result;
}

TestResult testWeightedColoring(const BAGraph& graph) {
  TestResult result;

  BAGraph graph1 = graph;
  while (graph1.getEdgesNumber() > 0) {
    BAGraph max_planar_subgraph_mst;
    MaxPlanarSubgraph::weightedMstBased(graph1, max_planar_subgraph_mst);
    MaxPlanarSubgraph::weightedMaximizeSubgraph(graph1, max_planar_subgraph_mst);
    result.mstRatio += 1.0;

    for (const auto& edge : max_planar_subgraph_mst.getEdges()) {
      graph1.removeEdge(edge.source, edge.target);
    }
  }

  BAGraph graph2 = graph;
  while (graph2.getEdgesNumber() > 0) {
    BAGraph max_planar_subgraph_cactus;
    MaxPlanarSubgraph::weightedCactusBased(graph2, max_planar_subgraph_cactus);
    MaxPlanarSubgraph::weightedMaximizeSubgraph(graph2, max_planar_subgraph_cactus);
    result.cactusRatio += 1.0;

    for (const auto& edge : max_planar_subgraph_cactus.getEdges()) {
      graph2.removeEdge(edge.source, edge.target);
    }
  }

  BAGraph graph3 = graph;
  while (graph3.getEdgesNumber() > 0) {
    BAGraph max_planar_subgraph_boyerMyrvold;
    PlanarityTest::boyerMyrvoldPlanarSubgraph(graph3, max_planar_subgraph_boyerMyrvold);
    MaxPlanarSubgraph::weightedMaximizeSubgraph(graph3, max_planar_subgraph_boyerMyrvold);
    result.boyerMyrvoldRatio += 1.0;

    for (const auto& edge : max_planar_subgraph_boyerMyrvold.getEdges()) {
      graph3.removeEdge(edge.source, edge.target);
    }
  }

  return result;
}

TestResult testLiveColoring(size_t initial_nodes, size_t final_nodes, size_t edges_per_vertex,
                            const std::vector<ColorType>& colors, MetricCalculator::MetricType metric_type) {
  TestResult result;

  AvoidingKuratowskiGraphsPainter painter{colors, metric_type};
  RandomGraphFactory random_graph_factory;

  auto start = std::chrono::high_resolution_clock::now();
  const BAGraph graph = random_graph_factory.createBarabasiAlbertWithPreferentialAttachmentBatageljBrandes(
      initial_nodes, final_nodes, edges_per_vertex, &painter);
  auto end = std::chrono::high_resolution_clock::now();

  result.mstTime = std::chrono::duration<double, std::milli>(end - start).count();
  size_t crossing_edges_sum = 0;
  for (const auto& embedding : painter.getEmbeddings()) {
    crossing_edges_sum += MaxPlanarSubgraph::crossingEdges(embedding.second);
  }
  result.mstEdges = static_cast<double>(crossing_edges_sum) / static_cast<double>(painter.getEmbeddings().size());

  return result;
}

TestResult testLiveColoring2(size_t initial_nodes, size_t final_nodes, size_t edges_per_vertex,
                             const std::vector<ColorType>& colors) {
  TestResult result;

  BAGraph graph = RandomGraphFactory().createBarabasiAlbertWithPreferentialAttachmentBatageljBrandes(initial_nodes, final_nodes,
                                                                                                     edges_per_vertex);

  auto start = std::chrono::high_resolution_clock::now();
  std::vector<BAGraph> subgraphs = PreferentialColoring::color(graph, colors.size());
  auto end = std::chrono::high_resolution_clock::now();

  result.mstTime = std::chrono::duration<double, std::milli>(end - start).count();
  size_t crossing_edges_sum = 0;
  for (const auto& subgraph : subgraphs) {
    crossing_edges_sum += MaxPlanarSubgraph::crossingEdges(subgraph);
  }
  result.mstEdges = static_cast<double>(crossing_edges_sum) / static_cast<double>(subgraphs.size());

  return result;
}

void baGenerationBenchmark() {
  constexpr size_t initial_nodes = 10;

  constexpr size_t min_m = 2;
  constexpr size_t max_m = 5;
  constexpr size_t step_m = 1;

  constexpr size_t min_n = 50;
  constexpr size_t max_n = 1000;
  constexpr size_t step_n = 50;

  constexpr size_t repetitions = 100000;

  for (size_t m = min_m; m <= max_m; m += step_m) {
    std::cout << "m = " << m << std::endl;

    std::ofstream ba_times_file("ba_time_m_" + std::to_string(m) + ".txt");

    for (size_t n = min_n; n <= max_n; n += step_n) {
      std::cout << "n = " << n << std::endl;

      auto start = std::chrono::high_resolution_clock::now();
      for (size_t i = 0; i < repetitions; ++i) {
        RandomGraphFactory random_graph_factory;
        const BAGraph graph =
            random_graph_factory.createBarabasiAlbertWithPreferentialAttachmentBatageljBrandes(initial_nodes, n, m);
      }
      auto end = std::chrono::high_resolution_clock::now();

      const double average = std::chrono::duration<double, std::milli>(end - start).count() / static_cast<double>(repetitions);

      ba_times_file << n << " " << average << std::endl;
    }
  }
}

void maxPlanarSubgraphBenchmark() {
  constexpr size_t initial_nodes = 10;

  constexpr size_t min_m = 2;
  constexpr size_t max_m = 5;
  constexpr size_t step_m = 1;

  constexpr size_t min_n = 50;
  constexpr size_t max_n = 500;
  constexpr size_t step_n = 50;

  constexpr size_t repetitions = 100;

  for (size_t m = min_m; m <= max_m; m += step_m) {
    std::cout << "m = " << m << std::endl;

    std::ofstream base_edges_file("base_max_planar_subgraph_edges_m_" + std::to_string(m) + ".txt");
    std::ofstream base_ratios_file("base_max_planar_subgraph_ratios_m_" + std::to_string(m) + ".txt");
    std::ofstream base_times_file("base_max_planar_subgraph_times_m_" + std::to_string(m) + ".txt");

    std::ofstream edges_file("max_planar_subgraph_edges_m_" + std::to_string(m) + ".txt");
    std::ofstream ratios_file("max_planar_subgraph_ratios_m_" + std::to_string(m) + ".txt");
    std::ofstream times_file("max_planar_subgraph_times_m_" + std::to_string(m) + ".txt");

    for (size_t n = min_n; n <= max_n; n += step_n) {
      std::cout << "n = " << n << std::endl;
      TestResult baseBasicTotal;
      TestResult baseWeightedTotal;
      TestResult basicTotal;
      TestResult weightedTotal;

      for (size_t i = 0; i < repetitions; ++i) {
        RandomGraphFactory random_graph_factory;
        const BAGraph graph =
            random_graph_factory.createBarabasiAlbertWithPreferentialAttachmentBatageljBrandes(initial_nodes, n, m);

        const TestResult baseBasicResult = testBasicMethods(graph, false);
        baseBasicTotal = baseBasicTotal + baseBasicResult;

        const TestResult baseWeightedResult = testWeightedMethods(graph, false);
        baseWeightedTotal = baseWeightedTotal + baseWeightedResult;

        const TestResult basicResult = testBasicMethods(graph, true);
        basicTotal = basicTotal + basicResult;

        const TestResult weightedResult = testWeightedMethods(graph, true);
        weightedTotal = weightedTotal + weightedResult;
      }

      const TestResult baseBasicAverage = baseBasicTotal / static_cast<double>(repetitions);
      const TestResult baseWeightedAverage = baseWeightedTotal / static_cast<double>(repetitions);

      const TestResult basicAverage = basicTotal / static_cast<double>(repetitions);
      const TestResult weightedAverage = weightedTotal / static_cast<double>(repetitions);

      base_edges_file << n << " " << baseBasicAverage.boyerMyrvoldEdges << " " << baseWeightedAverage.boyerMyrvoldEdges << " "
                      << baseBasicAverage.mstEdges << " " << baseWeightedAverage.mstEdges << " " << baseBasicAverage.cactusEdges
                      << " " << baseWeightedAverage.cactusEdges << std::endl;
      base_ratios_file << n << " " << baseBasicAverage.boyerMyrvoldRatio << " " << baseWeightedAverage.boyerMyrvoldRatio << " "
                       << baseBasicAverage.mstRatio << " " << baseWeightedAverage.mstRatio << " " << baseBasicAverage.cactusRatio
                       << " " << baseWeightedAverage.cactusRatio << std::endl;
      base_times_file << n << " " << baseBasicAverage.boyerMyrvoldTime << " " << baseWeightedAverage.boyerMyrvoldTime << " "
                      << baseBasicAverage.mstTime << " " << baseWeightedAverage.mstTime << " " << baseBasicAverage.cactusTime
                      << " " << baseWeightedAverage.cactusTime << std::endl;

      edges_file << n << " " << basicAverage.boyerMyrvoldEdges << " " << weightedAverage.boyerMyrvoldEdges << " "
                 << basicAverage.mstEdges << " " << weightedAverage.mstEdges << " " << basicAverage.cactusEdges << " "
                 << weightedAverage.cactusEdges << std::endl;
      ratios_file << n << " " << basicAverage.boyerMyrvoldRatio << " " << weightedAverage.boyerMyrvoldRatio << " "
                  << basicAverage.mstRatio << " " << weightedAverage.mstRatio << " " << basicAverage.cactusRatio << " "
                  << weightedAverage.cactusRatio << std::endl;
      times_file << n << " " << basicAverage.boyerMyrvoldTime << " " << weightedAverage.boyerMyrvoldTime << " "
                 << basicAverage.mstTime << " " << weightedAverage.mstTime << " " << basicAverage.cactusTime << " "
                 << weightedAverage.cactusTime << std::endl;
    }
  }
}

void coloringBenchmark() {
  constexpr size_t initial_nodes = 10;

  constexpr size_t min_m = 2;
  constexpr size_t max_m = 5;
  constexpr size_t step_m = 1;

  constexpr size_t min_n = 50;
  constexpr size_t max_n = 500;
  constexpr size_t step_n = 25;

  constexpr size_t repetitions = 100;

  for (size_t m = min_m; m <= max_m; m += step_m) {
    std::cout << "m = " << m << std::endl;

    std::ofstream ratios_file("coloring_ratios_m_" + std::to_string(m) + ".txt");

    for (size_t n = min_n; n <= max_n; n += step_n) {
      std::cout << "n = " << n << std::endl;
      TestResult basicTotal;
      TestResult weightedTotal;

      for (size_t i = 0; i < repetitions; ++i) {
        RandomGraphFactory random_graph_factory;
        const BAGraph graph =
            random_graph_factory.createBarabasiAlbertWithPreferentialAttachmentBatageljBrandes(initial_nodes, n, m);

        const TestResult basicResult = testColoring(graph);
        basicTotal = basicTotal + basicResult;

        const TestResult weightedResult = testWeightedColoring(graph);
        weightedTotal = weightedTotal + weightedResult;
      }

      const TestResult basicAverage = basicTotal / static_cast<double>(repetitions);
      const TestResult weightedAverage = weightedTotal / static_cast<double>(repetitions);

      ratios_file << n << " " << basicAverage.boyerMyrvoldRatio << " " << weightedAverage.boyerMyrvoldRatio << " "
                  << basicAverage.mstRatio << " " << weightedAverage.mstRatio << " " << basicAverage.cactusRatio << " "
                  << weightedAverage.cactusRatio << std::endl;
    }
  }
}

void liveColoringBenchmark() {
  constexpr size_t initial_nodes = 5;

  constexpr size_t min_m = 2;
  constexpr size_t max_m = 5;
  constexpr size_t step_m = 1;

  constexpr size_t min_n = 100;
  constexpr size_t max_n = 500;
  constexpr size_t step_n = 50;

  constexpr size_t repetitions = 100;

  for (size_t m = min_m; m <= max_m; m += step_m) {
    std::cout << "m = " << m << std::endl;

    std::ofstream times_file("live_coloring_times_m_" + std::to_string(m) + ".txt");
    std::ofstream crossing_edges_file("live_coloring_crossing_edges_m_" + std::to_string(m) + ".txt");

    for (size_t n = min_n; n <= max_n; n += step_n) {
      std::cout << "n = " << n << std::endl;
      TestResult total;

      for (size_t i = 0; i < repetitions; ++i) {
        std::vector<ColorType> colors(allColors.begin(), allColors.begin() + m);
        TestResult result = testLiveColoring(initial_nodes, n, m, colors, MetricCalculator::MetricType::EXPECTED);

        total = total + result;
      }

      const TestResult average = total / static_cast<double>(repetitions);
      times_file << n << " " << average.mstTime << std::endl;
      crossing_edges_file << n << " " << average.mstEdges << std::endl;
    }
  }
}

void liveColoring2Benchmark() {
  constexpr size_t initial_nodes = 5;

  constexpr size_t min_m = 3;
  constexpr size_t max_m = 5;
  constexpr size_t step_m = 1;

  constexpr size_t min_n = 50;
  constexpr size_t max_n = 500;
  constexpr size_t step_n = 5;

  constexpr size_t repetitions = 5;

  for (size_t m = min_m; m <= max_m; m += step_m) {
    std::cout << "m = " << m << std::endl;

    std::ofstream times_file("live_coloring2_times_m_" + std::to_string(m) + ".txt");
    std::ofstream crossing_edges_file("live_coloring2_crossing_edges_m_" + std::to_string(m) + ".txt");

    for (size_t n = min_n; n <= max_n; n += step_n) {
      std::cout << "n = " << n << std::endl;
      TestResult total;

      for (size_t i = 0; i < repetitions; ++i) {
        std::vector<ColorType> colors(allColors.begin(), allColors.begin() + 2);
        TestResult result = testLiveColoring2(initial_nodes, n, m, colors);

        total = total + result;
      }

      const TestResult average = total / static_cast<double>(repetitions);
      times_file << n << " " << average.mstTime << std::endl;
      crossing_edges_file << n << " " << average.mstEdges << std::endl;
    }
  }
}

void crossingEdgesEstimation() {
  constexpr size_t min_m = 3;
  constexpr size_t max_m = 5;
  constexpr size_t step_m = 1;

  constexpr size_t min_n = 100;
  constexpr size_t max_n = 10000;
  constexpr size_t step_n = 500;

  constexpr size_t repetitions = 1000;

  std::random_device rd;

  for (size_t m = min_m; m <= max_m; m += step_m) {
    std::cout << "m = " << m << std::endl;

    for (size_t n = min_n; n <= max_n; n += step_n) {
      std::vector<size_t> indexes(m * n);
      std::iota(indexes.begin(), indexes.end(), 1);

      double total = 0.0;
      for (size_t rep = 0; rep < repetitions; ++rep) {
        std::shuffle(indexes.begin(), indexes.end(), rd);

        std::vector<size_t> five_indexes(indexes.begin(), indexes.begin() + 5);
        std::sort(five_indexes.begin(), five_indexes.end());

        double sum = 0.0;
        sum += 16.0 / static_cast<double>(five_indexes[0]);
        sum += 49.0 / static_cast<double>(five_indexes[1]);
        sum += 64.0 / static_cast<double>(five_indexes[2]);
        sum += 49.0 / static_cast<double>(five_indexes[3]);
        sum += 16.0 / static_cast<double>(five_indexes[4]);

        total += sum;
      }

      std::cout << n << " " << total / static_cast<double>(repetitions) << std::endl;
    }
  }
}

int main() {
  baGenerationBenchmark();
  maxPlanarSubgraphBenchmark();
  coloringBenchmark();
  liveColoringBenchmark();
  liveColoring2Benchmark();

  crossingEdgesEstimation();

  return 0;
}