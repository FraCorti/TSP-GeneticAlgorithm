cmake_minimum_required(VERSION 3.16)
project(GeneticTSP VERSION 1.0
        DESCRIPTION "Genetic Algorithms project"
        LANGUAGES CXX)

# find package needed by the project
find_package(Threads REQUIRED)
find_package(OpenMP REQUIRED)

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
        src/tsp/tsp.h
        src/tsp/sequentialTSP.h
        src/tsp/parallelTSP.h
        src/tsp/fastflowTSP.h)

# set compiler
target_compile_features(GeneticTSP PUBLIC cxx_std_17)
set_target_properties(GeneticAlgorithms PROPERTIES
        CXX_STANDARD 17
        LINKER_LANGUAGE CXX)

# link libraries to main target
target_link_libraries(GeneticTSP
        GeneticAlgorithms
        OpenMP::OpenMP_CXX
        Threads::Threads
        )
