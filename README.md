# A Stabilized Numerical Framework for Necrotic Tumor Growth

This repository contains the code for the paper *“A Stabilized Numerical Framework for Necrotic Tumor Growth via Coupled Boundary Integral and Obstacle Solvers”*.

## C++ Projects Collection

This repository contains multiple independent C++ projects. Each project is organized in its own directory and uses CMake for building.

### Project Structure
tumorWithoutNecroticCore/

tumorWithNecroticCore/ 

tumorWithNecroticCoreEmergency/

Each project contains a `build` directory for out-of-source builds.

## Requirements

Make sure you have the following tools installed on your `Linux` system:

- CMake
  
- GCC / G++ (or another C++ compiler)

- Make

## Build Instructions

For each project, navigate to its `build` directory and run:

`cmake .. && make -j`

After building, the executable file will be generated in the `build' directory:

`./test`

### Note

For each project, data for running `./test` will be stored in `build/result/points` and `build/result/solutions`, and there are several `.ipynb` Python files in `build` or `build/result` directory for visualization. 
