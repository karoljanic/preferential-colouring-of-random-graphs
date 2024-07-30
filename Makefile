CXX = clang++
CXXFLAGS = -O3

all: generategraph generatecolorgraph

generategraph: random_graph_factory.cpp generate_graph.cpp
	$(CXX) -std=c++20 -o generategraph $(CXXFLAGS) random_graph_factory.cpp generate_graph.cpp

generatecolorgraph: random_graph_factory.cpp avoiding_kuratowski_graphs_painter.cpp generate_color_graph.cpp
	$(CXX) -std=c++20 -o generatecolorgraph $(CXXFLAGS) random_graph_factory.cpp avoiding_kuratowski_graphs_painter.cpp generate_color_graph.cpp

clean:
	rm -f generategraph generatecolorgraph