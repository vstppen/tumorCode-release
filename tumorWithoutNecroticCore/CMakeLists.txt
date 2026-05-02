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

# Get all source files from src/ folder
file(GLOB SRC_LIST "${PROJECT_SOURCE_DIR}/src/*.cpp")

set(EXECUTABLE_NAME "test")
set(MAIN_FILE "${PROJECT_SOURCE_DIR}/src/test.cpp")

list(REMOVE_ITEM SRC_LIST "${PROJECT_SOURCE_DIR}/src/test.cpp")

# Add the selected main file
list(APPEND SRC_LIST "${MAIN_FILE}")

# Create executable
add_executable(${EXECUTABLE_NAME} ${SRC_LIST})

# Link GSL and CBLAS libraries
target_link_libraries(${EXECUTABLE_NAME} gsl gslcblas)
