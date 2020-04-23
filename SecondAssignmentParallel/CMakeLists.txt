cmake_minimum_required(VERSION 3.11)
project(SecondAssignmentParallel VERSION 1.0
        DESCRIPTION "Second Assignment"
                                 LANGUAGES CXX)

# find package needed by the project
find_package(Threads REQUIRED)
find_package(OpenMP REQUIRED)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# different build requires different flags
set(CMAKE_CXX_FLAGS_DEBUG "-Wall -Wextra -O0 -fsanitize=thread -fstack-check -fstack-protector-strong -fno-omit-frame-pointer -g -DDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -ftree-vectorize") # vectorization info -fopt-info-vec

# creating the target SecondAssignmentParallel
add_executable(SecondAssignmentParallel main.cpp)
add_library(MasterWorker MasterWorker.h)

# set compiler
target_compile_features(SecondAssignmentParallel PUBLIC cxx_std_17)
set_target_properties(MasterWorker PROPERTIES
        CXX_STANDARD 17
        LINKER_LANGUAGE CXX)

# link libraries
target_link_libraries(SecondAssignmentParallel
        MasterWorker
        OpenMP::OpenMP_CXX
        Threads::Threads
        )