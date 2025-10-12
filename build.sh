#!/bin/bash

# You have to install manually Unittest++ library for your system
# eg. for Ubuntu: sudo apt-get install libunittest++-dev
# eg. for Pacman: sudo pacman -S unittestpp

# C++ implementation part building
rm -rf build
mkdir build
cd build
cmake ../code
make

# Peregrine project building
cd ../build/_deps/peregrine-src || exit
source tbb2020/bin/tbbvars.sh intel64
make -j CC=g++
