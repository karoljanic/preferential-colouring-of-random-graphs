cmake_minimum_required(VERSION 3.15)

project(Graphs)

set(OUTPUT "../graphs")

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

set(SOURCE_FILES
        random_graph_factory.cpp
        colors_balance_globally_painter.cpp
        colors_balance_locally_painter.cpp
        avoiding_kuratowski_graphs_painter.cpp
        metric_calculator.cpp
        max_planar_subgraph.cpp
        preferential_coloring.cpp)

add_executable(Graphs ${SOURCE_FILES} main.cpp)
add_executable(Benchmark ${SOURCE_FILES} benchmark.cpp)

set_target_properties(Graphs PROPERTIES OUTPUT_NAME ${OUTPUT})

target_link_libraries(Graphs PRIVATE nlohmann_json::nlohmann_json OGDF)
target_link_libraries(Benchmark PRIVATE nlohmann_json::nlohmann_json OGDF)

