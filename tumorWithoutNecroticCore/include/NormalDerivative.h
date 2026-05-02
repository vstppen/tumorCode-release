# ifndef NORMALDERIVATIVE_H
# define NORMALDERIVATIVE_H

class CartesianGridAndControlPoints;

void findInteriorPoints(CartesianGridAndControlPoints* G, double (*g)(double, double), int bdry_idx, double** u, double& x, double& y, double** coord, double* grid_u, int& count);
void computeBoundaryValues(double x, double y, double** coord, double* grid_u, int count, double& u, double& ux, double& uy);
void extractboundaryData(CartesianGridAndControlPoints* G, double (*g)(double, double), double** w, double* ux, double* uy);
void getNormalDerivatives(CartesianGridAndControlPoints* G, double (*g)(double, double), double** u, double* v);

# endif