cmake_minimum_required(VERSION 3.10)

# Project name
project(test)

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Compiler optimization flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -fopenmp")

# Include directories (add header files path)
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(/usr/include/eigen3)

# Get all source files from src/ folder
file(GLOB SRC_LIST "${PROJECT_SOURCE_DIR}/src/*.cpp")

set(EXECUTABLE_NAME "test")

# Create executable
add_executable(${EXECUTABLE_NAME} ${SRC_LIST})

# Link GSL and CBLAS libraries
target_link_libraries(${EXECUTABLE_NAME} gsl gslcblas)
