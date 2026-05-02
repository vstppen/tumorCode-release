# ifndef GMRES_H
# define GMRES_H

class CartesianGridAndControlPoints;

bool GMRES(CartesianGridAndControlPoints* G0, CartesianGridAndControlPoints* G1, const double *b, double *u, int max_m, int max_itr_num, double tol, int &itr_num);

# endif