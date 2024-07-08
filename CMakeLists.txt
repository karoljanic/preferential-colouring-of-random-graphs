cmake_minimum_required(VERSION 3.15)

project(Graphs)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(SOURCE_FILES
        main.cpp
        random_graph_factory.cpp
        colors_balance_globally_painter.cpp
        colors_balance_locally_painter.cpp
        avoiding_kuratowski_graphs_painter.cpp)

add_executable(Graphs ${SOURCE_FILES})
