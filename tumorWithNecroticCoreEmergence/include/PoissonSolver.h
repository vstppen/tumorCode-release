# ifndef POISSONSOLVER_H
# define POISSONSOLVER_H

class CartesianGridAndControlPoints0;

void solvePoisson(CartesianGridAndControlPoints0* G, double** source, double** numerical_solution, double* numerical_v, const std::string& filename_sol, const std::string& filename_v);

# endif