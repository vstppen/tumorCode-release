# ifndef MHELMHOLTZSOLVER_H
# define MHELMHOLTZSOLVER_H

class CartesianGridAndControlPoints0;

void solveModifiedHelmholtz(CartesianGridAndControlPoints0* G, double** numerical_solution, const std::string& filename);

# endif