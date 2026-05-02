# ifndef POISSONSOLVER_H
# define POISSONSOLVER_H

class CartesianGridAndControlPoints;

extern double G0;

double g_Poisson(double x0, double y0); 
double F_Poisson(double x_min, double y_min, double dx, double dy, int I, int J, double** source, double x0, double y0);

void VolumeIntegralPoisson(CartesianGridAndControlPoints* G, double** source, double** w, double* bdry_w);

double KernelBdryPoisson(double x1_t, double x2_t, double dx1_t, double dx2_t, double ddx1_t, double ddx2_t, double x1_s, double x2_s, double dx1_s, double dx2_s, double ddx1_s, double ddx2_s);
double KernelInteriorPoisson(double x, double y, double x1_s, double x2_s, double dx1_s, double dx2_s, double ddx1_s, double ddx2_s); 
double BoundaryIntegralPoisson(CartesianGridAndControlPoints* G, double* phi, int idx_t); 
double BoundaryIntegralPoisson(CartesianGridAndControlPoints* G, double* phi, int idx_t, double** K);
void solveBIEPoisson(CartesianGridAndControlPoints* G, double** source, double* phi, double** w);

double computeNearBoundaryDoubleLayerPotential(CartesianGridAndControlPoints* G, double* phi, double* M__, double* alpha_, double* beta_, double x, double y);

void solvePoisson(CartesianGridAndControlPoints* G, double** source, double** numerical_solution, double* numerical_v, const std::string& filename_sol, const std::string& filename_v);

# endif