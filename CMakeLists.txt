cmake_minimum_required(VERSION 3.11)
project(GeneticTSP VERSION 1.0
        DESCRIPTION "TSP Genetic Algorithm"
        LANGUAGES CXX)

# find packages
find_package(Threads REQUIRED)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif ()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# different build requires different flags
set(CMAKE_CXX_FLAGS_DEBUG "-Wall -Wextra -O0 -fsanitize=thread -fstack-check -fstack-protector-strong -fno-omit-frame-pointer -g -DDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -ftree-vectorize") # vectorization info -fopt-info-vec

# creating targets GeneticTSP and GeneticAlgorithms
add_executable(GeneticTSP main.cpp)
add_library(GeneticAlgorithms
        STATIC
        src/graph/graph.h   
        src/tsp/sequentialTSP.h
        src/tsp/fastflowTSP.h
        src/tsp/geneticAlgorithm.h
        src/tsp/parallelTSP.h)

set_property(TARGET GeneticTSP PROPERTY CXX_STANDARD 17)
set_target_properties(GeneticAlgorithms PROPERTIES
        CXX_STANDARD 17
        LINKER_LANGUAGE CXX)

# link libraries to main target
target_link_libraries(GeneticTSP
        GeneticAlgorithms
        Threads::Threads
        )