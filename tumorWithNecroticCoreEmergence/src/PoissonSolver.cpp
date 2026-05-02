#include <iostream>
#include <math.h>
#include <string.h>
#include <cstring>
#include <omp.h>    // OpenMP
// #include <vector>      // 用于动态数组 std::vector
// #include <utility>     // 用于 std::pair
#include <chrono>
#include <algorithm>  // For std::sort

#include "LinearAlgebra.h"
#include "CartesianGridAndControlPoints0.h"
#include "PoissonSolver.h"
#include "FastAlgorithm.h"
#include "Const.h"
#include "LU.h"
#include "NormalDerivative.h"
#include "CubicSpline.h"
#include "Config.h"

const double beta = 0.8;

double g_Poisson(double x0, double y0) 
{   
    return 0;
}

double F_Poisson(double x_min, double y_min, double dx, double dy, int I, int J, double** source, double x0, double y0)
{   
    int i = static_cast<int>((x0 - x_min) / dx + 1e-8); 
    int j = static_cast<int>((y0 - y_min) / dy + 1e-8); 

    if (i < 0 or i > I or j < 0 or j > J) {
        std::cout << "Warning about getting F_Poisson!" << std::endl;
        exit(EXIT_FAILURE);
    }

    double x1 = x_min + i * dx;
    double x2 = x1 + dx;
    double y1 = y_min + j * dy;
    double y2 = y1 + dy;

    double f;
    
    if (i == I && j == J) {
        f = source[I][J];
    } else if (i == I && j != J) {
        f = source[i][j] * (y2 - y0) + source[i][j + 1] * (y0 - y1);
        f /= dy;
    } else if (j == J && i != I) {
        f = source[i][j] * (x2 - x0) + source[i + 1][j] * (x0 - x1);
        f /= dx;
    } else {
        f = source[i][j] * (x2 - x0) * (y2 - y0) + source[i + 1][j] * (x0 - x1) * (y2 - y0) + source[i + 1][j + 1] * (x0 - x1) * (y0 - y1) + source[i][j + 1] * (x2 - x0) * (y0 - y1);
        f /= (dx * dy);
    }

    return f;
}

void VolumeIntegralPoisson(CartesianGridAndControlPoints0* G, double** source, double** w, double* bdry_w)
{   
    int size[2] = {G->I + 1, G->J + 1};
    double h2 = G->dx * G->dy; 

    for (int i = 1; i < G->I; i++) {
        for (int j = 1; j < G->J; j++) {
            if (G->interior[i][j]) {
                w[i][j] = - h2 *  F_Poisson(G->x_min, G->y_min, G->dx, G->dy, G->I, G->J, source, G->Grid_x[i], G->Grid_y[j]); 
            }
        }
    }

    G->makeCorrection(F_Poisson, source, w); 
    
    FastPoissonSolver(w, size); // fast Fourier transform based solver

    G->extractDirichletBoundaryData(w, F_Poisson, source, bdry_w, G->M); 
}

double KernelBdryPoisson(double x1_t, double x2_t, double dx1_t, double dx2_t, double ddx1_t, double ddx2_t, double x1_s, double x2_s, double dx1_s, double dx2_s, double ddx1_s, double ddx2_s) 
{   
    double k = 0.;
    if (abs(x1_t - x1_s) > EPSILON or abs(x2_t - x2_s) > EPSILON) {
        double temp1 = x1_s - x1_t;
        double temp2 = x2_s - x2_t;
        k = (dx2_s * temp1 - dx1_s * temp2) / (temp1 * temp1 + temp2 * temp2);
    } else {
        k = 0.5 * (ddx2_t * dx1_t - ddx1_t * dx2_t) / (dx1_t * dx1_t + dx2_t * dx2_t);
    }
    return k / (2 * M_PI);
}

double KernelInteriorPoisson(double x, double y, double x1_s, double x2_s, double dx1_s, double dx2_s, double ddx1_s, double ddx2_s) 
{
    double r2 = (x1_s - x) * (x1_s - x) + (x2_s - y) * (x2_s - y);
    double result;
    if (r2 > EPSILON12) {
        result = (dx2_s * (x1_s - x) - dx1_s * (x2_s - y)) / r2;
    } else {
        // std::cout << "warning about computing dG " << x << " " << y << " " << x1_s << " " << x2_s << std::endl;
        result = KernelBdryPoisson(x1_s, x2_s, dx1_s, dx2_s, ddx1_s, ddx2_s, x1_s, x2_s, dx1_s, dx2_s, ddx1_s, ddx2_s) * 2 * M_PI;
    }
    return result / (2 * M_PI);
}

double BoundaryIntegralPoisson(CartesianGridAndControlPoints0* G, double* phi, int idx_t) 
{
    double result = 0.0;
    for (int i = 0; i < G->M; i++) {
        result += KernelBdryPoisson(G->xy[idx_t][0], G->xy[idx_t][1], G->dxdy[idx_t][0], G->dxdy[idx_t][1], 
                            G->ddxddy[idx_t][0], G->ddxddy[idx_t][1], G->xy[i][0], G->xy[i][1], 
                            G->dxdy[i][0], G->dxdy[i][1], G->ddxddy[i][0], G->ddxddy[i][1]) * phi[i];
    }
    return G->bdry_delta * result;
}

double BoundaryIntegralPoisson(CartesianGridAndControlPoints0* G, double* phi, int idx_t, double** K) 
{
    double result = 0.0;
    for (int i = 0; i < G->M; i++) {
        result += K[idx_t][i] * phi[i];
    }
    return G->bdry_delta * result;
}


void solveBIEPoisson(CartesianGridAndControlPoints0* G, double** source, double* phi, double** w)
{
    double** K = createMatrix(G->M, G->M);
    double* bdry_w = createVector(G->M);
    for (int i = 0; i < G->M; i++) {
        for (int j = 0; j < G->M; j++) {
            K[i][j] = KernelBdryPoisson(G->xy[i][0], G->xy[i][1], G->dxdy[i][0], G->dxdy[i][1], 
                            G->ddxddy[i][0], G->ddxddy[i][1], G->xy[j][0], G->xy[j][1], 
                            G->dxdy[j][0], G->dxdy[j][1], G->ddxddy[j][0], G->ddxddy[j][1]);
        }
    }

    double* gD = createVector(G->M);

    VolumeIntegralPoisson(G, source, w, bdry_w);

    for (int i = 0; i < G->M; i ++) {
        gD[i] = g_Poisson(G->xy[i][0], G->xy[i][1]) - bdry_w[i];
        phi[i] = 0.5 * gD[i];
    }

    double rtol = 1.0E-12; // relative tolerance
    double res_norm0 = computeMaxNorm_L(gD, G->M);

    double atol = res_norm0 * rtol; // absolute tolerance
    int max_itr_num = 10000;
    
    double* r = createVector(G->M);
    double* integral = createVector(G->M);
    int count = 0;
    double res_norm = 0.0;

    for (int i = 0; i < G->M; i ++) {
        integral[i] = BoundaryIntegralPoisson(G, phi, i, K);
    }

    for (int i = 0; i < G->M; i ++) {
        r[i] = gD[i] - 0.5 * phi[i] - integral[i];
    }

    do {
        count++;
        
        for (int i = 0; i < G->M; i ++) {
            phi[i] = phi[i] + 2 * beta * r[i];
        }

        for (int i = 0; i < G->M; i ++) {
            integral[i] = BoundaryIntegralPoisson(G, phi, i, K);
        }

        for (int i = 0; i < G->M; i ++) {
            r[i] = gD[i] - 0.5 * phi[i] - integral[i];
        }
        
        res_norm = computeMaxNorm_L(r, G->M);

        // std::cout << "n = " << n << ". After iterating " 
        // << count << " times, the residual-norm is " << res_norm << "." << std::endl;
            
    } while ((res_norm > atol) && (count < max_itr_num)) ;

    if (count == max_itr_num) {
        std::cout << "The Richardson iteration failed to converge!" << std::endl;
        exit(EXIT_FAILURE);
    }

    freeVector(r);
    freeVector(integral);
    freeVector(gD);
    freeVector(bdry_w);
    freeMatrix(K, G->M);
}

inline double E(double eta, double s)
{
    return exp(2 * eta * s) * std::erfc(eta + s) + exp(- 2 * eta * s) * std::erfc(- eta + s);
}

double computeNearBoundaryDoubleLayerPotential(CartesianGridAndControlPoints0* G, double* phi, double* M__, double* alpha_, double* beta_, double x, double y) {
    double delta = G->bdry_delta;
    double alpha0, x0, y0;
    G->findClosestPoint(x, y, alpha0, x0, y0);
    double dx0, dy0, ddx0, ddy0;
    G->getPoint3(alpha0, x0, y0, dx0, dy0, ddx0, ddy0);

    int idx = static_cast<int>(alpha0 / G->bdry_delta);
    // if (idx >= G->M) idx = G->M - 1; 
    int idx_next = (idx + 1) % G->M;
    int idx_before = (idx - 1 + G->M) % G->M;

    double t = (alpha0 - idx * G->bdry_delta) / G->bdry_delta;

    double phi0 = periodicCubicSplineGetS(M__, alpha_, beta_, G->M, alpha0);
    double dphi0 = periodicCubicSplineGetDS(M__, alpha_, beta_, G->M, alpha0);
    double ddphi0 = periodicCubicSplineGetDDS(M__, alpha_, beta_, G->M, alpha0);

    double b = - sqrt((x - x0) * (x - x0) + (y - y0) * (y - y0));
    double eta = b / delta;
    double tau0 = sqrt(dx0 * dx0 + dy0 * dy0);

    double T1 = - delta * delta / (4 * M_PI) * eta * (sqrt(M_PI) * exp(- eta * eta) - M_PI * abs(eta) * std::erfc(abs(eta))) 
                * (ddphi0 / (tau0 * tau0) - (ddx0 * dx0 + ddy0 * dy0) * dphi0 / (tau0 * tau0 * tau0 * tau0)); 
    double T2 = 0;
    for (int n = 1; n < 5; n++) {
        T2 += sin(2 * n * M_PI * alpha0 / G->bdry_delta) * E(eta, n * M_PI * delta / (G->bdry_delta * tau0));
    }
    T2 *= -0.5 * G->bdry_delta * dphi0 * eta * delta / (G->bdry_delta * tau0);

    double result = 0.;
    double r2;
    for (int k = 0; k < G->M; k++) {
        r2 = (x - G->xy[k][0]) * (x - G->xy[k][0]) + (y - G->xy[k][1]) * (y - G->xy[k][1]);
        result += G->bdry_delta * (1 - exp(- r2 / (delta * delta))) * KernelInteriorPoisson(x, y, G->xy[k][0], G->xy[k][1], G->dxdy[k][0], G->dxdy[k][1], G->ddxddy[k][0], G->ddxddy[k][1]) * (phi[k] - phi0);
    }
    result += phi0 + T1 + T2;

    return result;
}

void solvePoisson(CartesianGridAndControlPoints0* G, double** source, double** numerical_solution, double* numerical_v, const std::string& filename_sol, const std::string& filename_v) 
{       
    double* phi = createVector(G->M);
    double** w = createMatrix(G->I + 1, G->J + 1);
    
    solveBIEPoisson(G, source, phi, w);
    
    int total_size = (G->I + 1) * (G->J + 1); 
    #pragma omp parallel for num_threads(omp_get_max_threads())
    for (int idx = 0; idx < total_size; idx++) {
        int i = idx / (G->J + 1); 
        int j = idx % (G->J + 1); 
        numerical_solution[i][j] = 0;

        if (G->interior[i][j] & !G->near_bdry[i][j]) {
            for (int k = 0; k < G->M; k++) {
                numerical_solution[i][j] += G->bdry_delta * KernelInteriorPoisson(G->Grid_x[i], G->Grid_y[j], G->xy[k][0], G->xy[k][1], G->dxdy[k][0], G->dxdy[k][1], G->ddxddy[k][0], G->ddxddy[k][1]) * phi[k];
            }
            numerical_solution[i][j] += w[i][j];
        }
    }

    double* M__ = createVector(G->M);
    double* alpha_  = createVector(G->M);
    double* beta_ = createVector(G->M);
    periodicCubicSplineInterpolation(phi, M__, alpha_, beta_, G->M);

    #pragma omp parallel for num_threads(omp_get_max_threads())
    for (int idx = 0; idx < total_size; idx++) {
        int i = idx / (G->J + 1); 
        int j = idx % (G->J + 1);
        if (G->interior[i][j] & G->near_bdry[i][j]) {
            numerical_solution[i][j] = computeNearBoundaryDoubleLayerPotential(G, phi, M__, alpha_, beta_, G->Grid_x[i], G->Grid_y[j]) + w[i][j];
        }
    }

    freeVector(M__);
    freeVector(beta_);
    freeVector(alpha_);

    saveMatrixToFile(const_cast<const double**>(numerical_solution), G->I + 1, G->J + 1, filename_sol);
    
    for(int i = 0; i < G->n_ctrl; i++) {
        numerical_v[i] = 0.;
    }
    
    getNormalDerivatives0(G, g_Poisson, numerical_solution, numerical_v);

    for(int i = 0; i < G->n_ctrl; i++) {
        numerical_v[i] = -numerical_v[i];
    }
    
    saveVectorToFile(const_cast<const double*>(numerical_v), G->n_ctrl, filename_v);
    
    freeVector(phi);
    freeMatrix(w, G->I + 1);
}