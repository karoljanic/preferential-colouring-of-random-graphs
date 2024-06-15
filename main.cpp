#include <iostream>

#include "sparse_graph.hpp"
#include "random_graph_factory.hpp"

int main() {
  std::cout << "Hello, World!" << std::endl;

  graph::random::RandomGraphFactory factory;
  const graph::random::BAGraph graph = factory.createBarabasiAlbertWithPreferentialAttachment(10, 10);

  return 0;
}
