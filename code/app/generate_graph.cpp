#include <iostream>  // std::cerr, std::cout
#include <fstream>   // std::ofstream

#include "../include/random_graph_factory.hpp"

int main(int argc, char* argv[]) {
  if (argc != 6) {
    std::cerr << "Usage: " << argv[0]
              << " <method> <initial_vertices_number> <edges_per_vertex> <final_vertices_number> <output_file>\n";
    return 0;
  }

  const int method = std::stoi(argv[1]);
  const size_t initial_vertices = std::stoi(argv[2]);
  const size_t edges_per_vertex = std::stoi(argv[3]);
  const size_t final_vertices = std::stoi(argv[4]);
  const std::string filename = argv[5];

  graph::random::BAGraph graph;
  graph::random::RandomGraphFactory random_graph_factory;
  switch (method) {
    case 1:
      graph = random_graph_factory.createBarabasiAlbertWithPreferentialAttachmentRepeatedNodes(initial_vertices, final_vertices,
                                                                                               edges_per_vertex);
      break;
    case 2:
      graph = random_graph_factory.createBarabasiAlbertWithPreferentialAttachmentBatageljBrandes(initial_vertices, final_vertices,
                                                                                                 edges_per_vertex);
      break;
    case 3:
      graph = random_graph_factory.createBarabasiAlbertWithLinkSelection(initial_vertices, final_vertices, edges_per_vertex);
      break;
    case 4:
      graph = random_graph_factory.createBarabasiAlbertWithCopyingModel(initial_vertices, final_vertices, edges_per_vertex, 0.5);
      break;
    case 5:
      graph = random_graph_factory.createBarabasiAlbertWithLCDModel(final_vertices, edges_per_vertex);
      break;
    default:
      std::cerr << "Invalid method. Choose between 1 and 5.\n";
      return 0;
  }

  std::ofstream file{filename};
  for(const auto& edge : graph.getEdges()) {
    file << edge.source << " " << edge.target << "\n";
  }
  file.close();

  return 0;
}