cmake_minimum_required(VERSION 3.16)
project(ThirdAssignmentParallel VERSION 1.0
        DESCRIPTION "Third parallel Assignment
        LANGUAGES CXX")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -ftree-vectorize")

find_package(Threads REQUIRED)

add_library(DAC divideAndConquer.cpp)
add_executable(ThirdAssignmentParallel main.cpp)
target_link_libraries(ThirdAssignmentParallel
        DAC
        Threads::Threads
        )