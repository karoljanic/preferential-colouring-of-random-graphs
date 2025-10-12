cmake_minimum_required(VERSION 3.5)

project(PreferentialColouringOfRandomGraphs)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

include(FetchContent)

FetchContent_Declare(json
        URL https://github.com/nlohmann/json/releases/download/v3.11.3/json.tar.xz)
FetchContent_MakeAvailable(json)

FetchContent_Declare(OGDF
        GIT_REPOSITORY https://github.com/ogdf/ogdf.git)
FetchContent_MakeAvailable(OGDF)

FetchContent_Declare(Peregrine
        GIT_REPOSITORY https://github.com/pdclab/peregrine.git)
FetchContent_MakeAvailable(Peregrine)

set(SOURCE_FILES
        src/random_graph_factory.cpp
        src/max_planar_subgraph.cpp
        src/preferential_coloring.cpp)

add_executable(GenerateGraph ${SOURCE_FILES} app/generate_graph.cpp)
# add_executable(MaximalPlanarSubgraphs ${SOURCE_FILES} app/maximal_planar_subgraphs.cpp)
# add_executable(MinimumColoring ${SOURCE_FILES} app/minimum_coloring.cpp)
# add_executable(MinimalCrossingEdges ${SOURCE_FILES} app/minimal_crossing_edges.cpp)

set_target_properties(GenerateGraph PROPERTIES OUTPUT_NAME "generate_graph")

target_link_libraries(GenerateGraph PRIVATE nlohmann_json::nlohmann_json OGDF)
# target_link_libraries(MaximalPlanarSubgraphs PRIVATE nlohmann_json::nlohmann_json OGDF)
# target_link_libraries(MinimumColoring PRIVATE nlohmann_json::nlohmann_json OGDF)
# target_link_libraries(MinimalCrossingEdges PRIVATE nlohmann_json::nlohmann_json OGDF)
