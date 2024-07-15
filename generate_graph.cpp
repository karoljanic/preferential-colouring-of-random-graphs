#include "ba_graph.hpp"              // graph::random::BAGraph
#include "random_graph_factory.hpp"  // graph::random::RandomGraphFactory

int main(int argc, char* argv[]) {
  if (argc != 5) {
    std::cerr << "Usage: " << argv[0] << " <initial_nodes> <edges_per_vertex> <final_nodes> <output_file>" << std::endl;
    return 1;
  }

  const size_t initial_nodes = std::stoi(argv[1]);
  const size_t edges_per_vertex = std::stoi(argv[2]);
  const size_t final_nodes = std::stoi(argv[3]);

  graph::random::RandomGraphFactory random_graph_factory;
  graph::random::BAGraph graph = random_graph_factory.createBarabasiAlbertWithPreferentialAttachmentBatageljBrandes(
      initial_nodes, final_nodes, edges_per_vertex);
  graph.saveToFile(argv[4]);

  return 0;
}
