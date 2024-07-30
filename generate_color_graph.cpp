#include "avoiding_kuratowski_graphs_painter.hpp"  // graph::random::AvoidingKuratowskiGraphsPainter
#include "ba_graph.hpp"                            // graph::random::BAGraph
#include "random_graph_factory.hpp"                // graph::random::RandomGraphFactory

constexpr graph::random::ColorType Blue{"#1f77b4"};
constexpr graph::random::ColorType Orange{"#ff7f0e"};
constexpr graph::random::ColorType Green{"#2ca02c"};
constexpr graph::random::ColorType Red{"#d62728"};
constexpr graph::random::ColorType Purple{"#9467bd"};
constexpr graph::random::ColorType Brown{"#8c564b"};
constexpr graph::random::ColorType Pink{"#e377c2"};
constexpr graph::random::ColorType Gray{"#7f7f7f"};
constexpr graph::random::ColorType Olive{"#bcbd22"};
constexpr graph::random::ColorType Cyan{"#17becf"};

std::vector<graph::random::ColorType> allColors{Blue, Orange, Green, Red, Purple, Brown, Pink, Gray, Olive, Cyan};

void saveWithColors(const graph::Graph<graph::random::BANode, graph::random::BAEdge>& graph, const std::string& filename,
                    bool node_coloring, bool edge_coloring) {

  auto node_color = node_coloring ?
                                  [](std::ofstream &file, const graph::random::BANode &node) {
                                    file << R"(, "color": ")" << node.color.hex << "\"";
                                  } : [](std::ofstream &/*file*/, const graph::random::BANode &/*node*/) {};

  auto edge_color = edge_coloring ?
                                  [](std::ofstream &file, const graph::random::BAEdge &edge) {
                                    file << R"(, "color": ")" << edge.color.hex << "\"";
                                  } : [](std::ofstream &/*file*/, const graph::random::BAEdge &/*edge*/) {};

  graph.saveToFile(filename, node_color, edge_color);
}

int main(int argc, char* argv[]) {
  if (argc != 6) {
    std::cerr << "Usage: " << argv[0] << " <initial_nodes> <edges_per_vertex> <final_nodes> <colors_number> <output_dir>"
              << std::endl;
    return 1;
  }

  const size_t initial_nodes = std::stoi(argv[1]);
  const size_t edges_per_vertex = std::stoi(argv[2]);
  const size_t final_nodes = std::stoi(argv[3]);
  const size_t colors_numbers = std::stoi(argv[4]);
  const std::vector<graph::random::ColorType> colors(allColors.begin(), allColors.begin() + colors_numbers);

  graph::random::AvoidingKuratowskiGraphsPainter painter{colors};

  graph::random::RandomGraphFactory random_graph_factory;
  graph::random::BAGraph graph = random_graph_factory.createBarabasiAlbertWithPreferentialAttachmentBatageljBrandes(
      initial_nodes, final_nodes, edges_per_vertex, &painter);

  saveWithColors(graph, argv[5] + std::string("/graph.json"), true, true);
  for (const auto& embedding : painter.getEmbeddings()) {
    saveWithColors(embedding.second, argv[5] + std::string("/embedding_") + std::string{embedding.first.hex} + ".txt", true, true);
  }

  return 0;
}
