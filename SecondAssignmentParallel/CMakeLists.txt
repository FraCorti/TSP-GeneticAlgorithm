cmake_minimum_required(VERSION 3.15)
project(SecondAssignmentParallel)

set(CMAKE_CXX_STANDARD 17)
find_package(OpenMP REQUIRED)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

set(CMAKE_CXX_FLAGS_DEBUG "-Wall -Wextra -O0 -fsanitize=thread -fstack-check -fstack-protector-strong -fno-omit-frame-pointer -g -DDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -fopt-info-vec")
set(CMAKE_CXX_FLAGS -pthread)
set(THREADS_PREFER_PTHREAD_FLAG ON)



add_executable(SecondAssignmentParallel main.cpp)
target_link_libraries(SecondAssignmentParallel OpenMP::OpenMP_CXX)