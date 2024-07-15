CXX = clang++
CXXFLAGS = -O3

generategraph: random_graph_factory.cpp generate_graph.cpp
	$(CXX) -std=c++20 -o generategraph $(CXXFLAGS) random_graph_factory.cpp generate_graph.cpp

clean:
	rm -f generategraph