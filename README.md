# A Stabilized Numerical Framework for Necrotic Tumor Growth

This repository contains the code for the paper *“A Stabilized Numerical Framework for Necrotic Tumor Growth via Coupled Boundary Integral and Obstacle Solvers”*.

## C++ Projects Collection

This repository contains multiple independent C++ projects. Each project is organized in its own directory and uses CMake for building.

### Project Structure
project1/
project2/
project3/

Each project contains a `build` directory for out-of-source builds.

## Requirements

Make sure you have the following tools installed on your `Linux' system:

- CMake
- GCC / G++ (or another C++ compiler)
- Make

## Build Instructions

For each project, navigate to its `build` directory and run:

```bash
cmake ..
make -j

After building, the executable file will be generated in the build directory:

./test
