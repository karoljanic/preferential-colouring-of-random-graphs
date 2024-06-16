cmake_minimum_required(VERSION 3.15)

project(Graphs)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(SOURCE_FILES main.cpp random_graph_factory.cpp)

add_executable(Graphs ${SOURCE_FILES})
