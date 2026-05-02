# ifndef NORMALDERIVATIVE_H
# define NORMALDERIVATIVE_H

class CartesianGridAndControlPoints;
class CartesianGridAndControlPoints0;

void getNormalDerivatives0(CartesianGridAndControlPoints0* G, double (*g)(double, double), double** u, double* v);
void getNormalDerivatives(CartesianGridAndControlPoints* G, double (*g)(double, double), double** u, double* v);

# endif