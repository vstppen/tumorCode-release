#ifndef PRIMALDUALALGORITHM_H
#define PRIMALDUALALGORITHM_H

class CartesianGridAndControlPoints;

int PrimalDualAlgorithm(CartesianGridAndControlPoints* G, double** source, double** numerical_solution, double* numerical_v1, const std::string& filename_sol, const std::string& filename_v1, double**& int_bd_points);


#endif