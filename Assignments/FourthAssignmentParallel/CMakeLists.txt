cmake_minimum_required(VERSION 3.11)
project(FourthAssignmentParallel VERSION 1.0
        DESCRIPTION "Fourth Assignment"
        LANGUAGES CXX)

include_directories("/home/marcod/opencv4/include/opencv4" "/home/marcod/opencv4/lib")
add_executable(FourthAssignmentParallel main.cpp ImageProcessing.h)

find_package(OpenCV REQUIRED)
target_include_directories(FourthAssignmentParallel
        PUBLIC
        ${OpenCV_INCLUDE_DIRS}
        )
#list(INSERT CMAKE_MODULE_PATH 0 ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
#find_package(OpenCV REQUIRED)

# different build requires different flags
set(CMAKE_CXX_FLAGS_DEBUG "-Wall -Wextra -O0 -fsanitize=thread -fstack-check -fstack-protector-strong -fno-omit-frame-pointer -g -DDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -ftree-vectorize") # vectorization info -fopt-info-vec

set_target_properties(FourthAssignmentParallel PROPERTIES
        CXX_STANDARD 17
        LINKER_LANGUAGE CXX)

target_link_libraries(FourthAssignmentParallel
        ${OpenCV_LIBS}
        )