cmake_minimum_required(VERSION 3.11)
project(FifthAssignmentParallel VERSION 1.0
        DESCRIPTION "Fifth Assignment"
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
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -finline-functions") # vectorization info -fopt-info-vec

# creating targets SecondAssignmentParallel and MasterWorker
add_library(FFCompressor ffc_pipe.cpp compdecomp.cpp utility.hpp)
add_executable(FifthAssignmentParallel main.cpp)

# add directory Miniz to miniz library
target_include_directories(FFCompressor
        PUBLIC
        miniz)

# set compiler
target_compile_features(FifthAssignmentParallel PUBLIC cxx_std_17)

# set c++17 for miniz
set_target_properties(FFCompressor PROPERTIES
        CXX_STANDARD 17
        LINKER_LANGUAGE CXX)

# link libraries to main target
target_link_libraries(FifthAssignmentParallel
        FFCompressor
        OpenMP::OpenMP_CXX
        Threads::Threads
        )