# ifndef MHELMHOLTZSOLVER_H
# define MHELMHOLTZSOLVER_H

class CartesianGridAndControlPoints;

extern double lambda;
extern double sqrt_lambda;
void update_sqrt_lambda();
extern double c_B;

double g_MHelmholtz(double x0, double y0);
double F_MHelmholtz(double x0, double y0);

// void VolumeIntegralMHelmholtz(CartesianGridAndControlPoints* G, double** w, double* bdry_w);

double KernelBdryMHelmholtz(double x1_t, double x2_t, double dx1_t, double dx2_t, double ddx1_t, double ddx2_t, double x1_s, double x2_s, double dx1_s, double dx2_s, double ddx1_s, double ddx2_s);
double KernelInteriorMHelmholtz(double x, double y, double x1_s, double x2_s, double dx1_s, double dx2_s);
double BoundaryIntegralMHelmholtz(CartesianGridAndControlPoints* G, double* phi, int idx_t);
double BoundaryIntegralMHelmholtz(CartesianGridAndControlPoints* G, double* phi, int idx_t, double** K);
void solveBIEMHelmholtz(CartesianGridAndControlPoints* G, double* phi, double** w);
void solveHomogeneousBIEMHelmholtz(CartesianGridAndControlPoints* G, double* phi);

double computeNearBoundaryValuesMHelmholtz(double** coor, double* u, double x, double y, int count);

void solveModifiedHelmholtz(CartesianGridAndControlPoints* G, double** numerical_solution, const std::string& filename);

# endif