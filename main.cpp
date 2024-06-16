#include <iostream>

#include "sparse_graph.hpp"
#include "random_graph_factory.hpp"

void printHistogram(const std::map<size_t, size_t> &histogram) {
  for (const auto &[degree, count] : histogram) {
	std::cout << "Degree: " << degree << " Count: " << count << std::endl;
  }
  std::cout << std::endl;
}

int main() {
  std::cout << "Hello, World!" << std::endl;

	const graph::random::BAGraph
		graph1 = graph::random::RandomGraphFactory::createBarabasiAlbertWithPreferentialAttachment(10, 50, 2, 1.2F);
	const graph::random::BAGraph
		graph2 = graph::random::RandomGraphFactory::createBarabasiAlbertWithLinkSelection(10, 50, 2);
	const graph::random::BAGraph
		graph3 = graph::random::RandomGraphFactory::createBarabasiAlbertWithCopyingModel(10, 50, 2, 0.5F);

  std::cout << "Graph 1 density: " << graph1.getDensity() << std::endl;
  std::cout << "Graph 1 degrees histogram:" << std::endl;
  printHistogram(graph1.getDegreesHistogram());

  std::cout << "Graph 2 density: " << graph2.getDensity() << std::endl;
  std::cout << "Graph 2 degrees histogram:" << std::endl;
  printHistogram(graph2.getDegreesHistogram());

  std::cout << "Graph 3 density: " << graph3.getDensity() << std::endl;
  std::cout << "Graph 3 degrees histogram:" << std::endl;
  printHistogram(graph3.getDegreesHistogram());

  return 0;
}
