# ifndef INTERFACEBVPSOLVER_H
# define INTERFACEBVPSOLVER_H

class CartesianGridAndControlPoints;

void Yi(CartesianGridAndControlPoints* G0, CartesianGridAndControlPoints* G1, double** w, double* bdry0_w, double* bdry0_wn);
void Ye(CartesianGridAndControlPoints* G0, CartesianGridAndControlPoints* G1, double** w, double* bdry0_w, double* bdry0_wn, double* bdry1_w, double* bdry1_wn);
void Vi(CartesianGridAndControlPoints* G0, CartesianGridAndControlPoints* G1, double* psi_vec, double** w, double* bdry0_w, double* bdry0_wn);
void Ve(CartesianGridAndControlPoints* G0, CartesianGridAndControlPoints* G1, int flag, double* psi_vec, double** w, double* bdry0_w, double* bdry0_wn, double* bdry1_w, double* bdry1_wn);
void Wi(CartesianGridAndControlPoints* G0, CartesianGridAndControlPoints* G1, double* phi_vec, double** w, double* bdry0_w, double* bdry0_wn);
void We(CartesianGridAndControlPoints* G0, CartesianGridAndControlPoints* G1, int flag, double* phi_vec, double** w, double* bdry0_w, double* bdry0_wn, double* bdry1_w, double* bdry1_wn);

void solveInterfacePDE(CartesianGridAndControlPoints* G0, CartesianGridAndControlPoints* G1, double** numerical_solution, const std::string& filename);

# endif