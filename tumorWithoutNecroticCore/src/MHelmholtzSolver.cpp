#include <iostream>
#include <math.h>
#include <string.h>
#include <cstring>
#include <gsl/gsl_sf_bessel.h>
#include <omp.h>    // OpenMP
#include <chrono>
#include <algorithm>  // For std::sort

#include "LinearAlgebra.h"
#include "CartesianGridAndControlPoints.h"
#include "MHelmholtzSolver.h"
#include "FastAlgorithm.h"
#include "Const.h"
#include "LU.h"

const double beta = 0.8;

double lambda = 1.;
double sqrt_lambda = sqrt(lambda);
void update_sqrt_lambda() {
    sqrt_lambda = sqrt(lambda);
}
double c_B = 1;

double g_MHelmholtz(double x0, double y0) 
{
  return c_B; 
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double F_MHelmholtz(double x0, double y0)
{
    return 0;
}

// void VolumeIntegralMHelmholtz(CartesianGridAndControlPoints* G, double** w, double* bdry_w)
// {   
//     int size[2] = {G->I + 1, G->J + 1};
//     double h2 = G->dx * G->dy; 

//     for (int i = 1; i < G->I; i++) {
//         for (int j = 1; j < G->J; j++) {
//             if (G->interior[i][j]) {
//                 w[i][j] = h2 * F_MHelmholtz(G->Grid_x[i], G->Grid_y[j]); 
//             }
//         }
//     }

//     G->makeCorrection(F_MHelmholtz, w); 
    
//     FastMHelmholtzSolver(w, size, lambda * h2); // fast Fourier transform based solver

//     G->extractDirichletBoundaryData(w, F_MHelmholtz, bdry_w, G->M); 

// }

double KernelBdryMHelmholtz(double x1_t, double x2_t, double dx1_t, double dx2_t, double ddx1_t, double ddx2_t, double x1_s, double x2_s, double dx1_s, double dx2_s, double ddx1_s, double ddx2_s) 
{   
    double k;
    double temp1 = x1_s - x1_t;
    double temp2 = x2_s - x2_t;
    double r2 = temp1 * temp1 + temp2 * temp2;
    double r = sqrt(r2);

    if (abs(temp1) > EPSILON or abs(temp2) > EPSILON) {
        k = gsl_sf_bessel_K1(sqrt_lambda * r) * sqrt_lambda * (dx2_s * temp1 - dx1_s * temp2) / r;
    } else {
        k = 0.5 * (ddx2_t * dx1_t - ddx1_t * dx2_t) / (dx1_t * dx1_t + dx2_t * dx2_t);
    }
    return k / (2 * M_PI);
}

double KernelInteriorMHelmholtz(double x, double y, double x1_s, double x2_s, double dx1_s, double dx2_s) 
{
    double temp1 = x1_s - x;
    double temp2 = x2_s - y;
    double r2 = temp1 * temp1 + temp2 * temp2;
    double r = sqrt(r2);

    double k = gsl_sf_bessel_K1(sqrt_lambda * r) * sqrt_lambda * (dx2_s * temp1 - dx1_s * temp2) / r;
    return k / (2 * M_PI);
}

double BoundaryIntegralMHelmholtz(CartesianGridAndControlPoints* G, double* phi, int idx_t) 
{
    double result = 0.0;
    for (int i = 0; i < G->M; i++) {
        result += KernelBdryMHelmholtz(G->xy[idx_t][0], G->xy[idx_t][1], G->dxdy[idx_t][0], G->dxdy[idx_t][1], 
                            G->ddxddy[idx_t][0], G->ddxddy[idx_t][1], G->xy[i][0], G->xy[i][1], 
                            G->dxdy[i][0], G->dxdy[i][1], G->ddxddy[i][0], G->ddxddy[i][1]) * phi[i];
    }
    return G->bdry_delta * result;
}

double BoundaryIntegralMHelmholtz(CartesianGridAndControlPoints* G, double* phi, int idx_t, double** K) 
{
    double result = 0.0;
    for (int i = 0; i < G->M; i++) {
        result += K[idx_t][i] * phi[i];
    }
    return G->bdry_delta * result;
}

void solveBIEMHelmholtz(CartesianGridAndControlPoints* G, double* phi, double** w)
{
    double** K = createMatrix(G->M, G->M);
    double* bdry_w = createVector(G->M);

    for (int i = 0; i < G->M; i++) {
        #pragma omp parallel for num_threads(omp_get_max_threads())
        for (int j = 0; j < G->M; j++) {
            K[i][j] = KernelBdryMHelmholtz(G->xy[i][0], G->xy[i][1], G->dxdy[i][0], G->dxdy[i][1], 
                            G->ddxddy[i][0], G->ddxddy[i][1], G->xy[j][0], G->xy[j][1], 
                            G->dxdy[j][0], G->dxdy[j][1], G->ddxddy[j][0], G->ddxddy[j][1]);
        }
    }

    double* gD = createVector(G->M);

    // VolumeIntegralMHelmholtz(G, w, bdry_w);

    // printMatrix(w, G->I + 1, G->J +1);
    // printVector(bdry_w, G->M);

    for (int i = 0; i < G->M; i ++) {
        gD[i] = g_MHelmholtz(G->xy[i][0], G->xy[i][1]) - bdry_w[i];
        phi[i] = 0.5 * gD[i];
    }

    double rtol = 1.0E-12; // relative tolerance
    double res_norm0 = computeMaxNorm(gD, G->M);

    double atol = res_norm0 * rtol; // absolute tolerance
    int max_itr_num = 10000;
    
    double* r = createVector(G->M);
    double* integral = createVector(G->M);
    int count = 0;
    double res_norm = 0.0;

    for (int i = 0; i < G->M; i ++) {
        integral[i] = BoundaryIntegralMHelmholtz(G, phi, i, K);
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
            integral[i] = BoundaryIntegralMHelmholtz(G, phi, i, K);
        }

        for (int i = 0; i < G->M; i ++) {
            r[i] = gD[i] - 0.5 * phi[i] - integral[i];
        }
        
        res_norm = computeMaxNorm(r, G->M);

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

void solveHomogeneousBIEMHelmholtz(CartesianGridAndControlPoints* G, double* phi)
{
    double** K = createMatrix(G->M, G->M);

    for (int i = 0; i < G->M; i++) {
        #pragma omp parallel for num_threads(omp_get_max_threads())
        for (int j = 0; j < G->M; j++) {
            K[i][j] = KernelBdryMHelmholtz(G->xy[i][0], G->xy[i][1], G->dxdy[i][0], G->dxdy[i][1], 
                            G->ddxddy[i][0], G->ddxddy[i][1], G->xy[j][0], G->xy[j][1], 
                            G->dxdy[j][0], G->dxdy[j][1], G->ddxddy[j][0], G->ddxddy[j][1]);
        }
    }

    double* gD = createVector(G->M);

    // VolumeIntegralMHelmholtz(G, w, bdry_w);

    // printMatrix(w, G->I + 1, G->J +1);
    // printVector(bdry_w, G->M);

    for (int i = 0; i < G->M; i ++) {
        gD[i] = g_MHelmholtz(G->xy[i][0], G->xy[i][1]);
        phi[i] = 0.5 * gD[i];
    }

    double rtol = 1.0E-12; // relative tolerance
    double res_norm0 = computeMaxNorm(gD, G->M);

    double atol = res_norm0 * rtol; // absolute tolerance
    int max_itr_num = 10000;
    
    double* r = createVector(G->M);
    double* integral = createVector(G->M);
    int count = 0;
    double res_norm = 0.0;

    for (int i = 0; i < G->M; i ++) {
        integral[i] = BoundaryIntegralMHelmholtz(G, phi, i, K);
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
            integral[i] = BoundaryIntegralMHelmholtz(G, phi, i, K);
        }

        for (int i = 0; i < G->M; i ++) {
            r[i] = gD[i] - 0.5 * phi[i] - integral[i];
        }
        
        res_norm = computeMaxNorm(r, G->M);

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
    freeMatrix(K, G->M);
}


double computeNearBoundaryValuesMHelmholtz(double** coor, double* u, double x, double y, int count) {
    int valid_count = count;
    if (count > 6) {
        valid_count = 6;
    }

    // Arrays to store distances and indices of the points
    double* distances = createVector(count);
    int* indices = createIntVector(count);

    // Calculate distances from (x, y) to each point in coor
    for (int i = 0; i < count; i++) {
        distances[i] = (x - coor[i][0]) * (x - coor[i][0]) + (y - coor[i][1]) * (y - coor[i][1]);
        indices[i] = i;  // Store original index of each point
    }

    // Sort points by distance to (x, y)
    std::sort(indices, indices + count, [&](int i1, int i2) {
        return distances[i1] < distances[i2];
    });

    double** A = createMatrix(valid_count, 6);
    double** AT = createMatrix(6, valid_count);
    double** ATA = createMatrix(6, 6);
    double* b = createVector(valid_count);
    double* ATb = createVector(6);
    double* sol = createVector(6);
    for (int i = 0; i < valid_count; i++) {
        int idx = indices[i];  // Get the index of the nearest point
        A[i][0] = coor[idx][0] * coor[idx][0];
        A[i][1] = coor[idx][0] * coor[idx][1];
        A[i][2] = coor[idx][1] * coor[idx][1];
        A[i][3] = coor[idx][0];
        A[i][4] = coor[idx][1];
        A[i][5] = 1;
        b[i] = u[idx];  // Corresponding value of u at the point
    }
 
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < valid_count; j++) {
            AT[i][j] = A[j][i];
        }
    }

    MatrixMulMatrix(AT, A, ATA, 6, valid_count, 6);
    MatrixMulVector(AT, b, ATb, 6, valid_count);
    SolveLinearSystemByCG(ATA, ATb, sol, 6, 0);

    double result = sol[0] * x * x + sol[1] * x * y + sol[2] * y * y + sol[3] * x + sol[4] * y + sol[5];
    
    freeMatrix(A, valid_count);
    freeMatrix(AT, 6);
    freeMatrix(ATA, 6);
    freeVector(b);
    freeVector(ATb);
    freeVector(sol);
    freeIntVector(indices);
    freeVector(distances);

    return result;
}

void solveModifiedHelmholtz(CartesianGridAndControlPoints* G, double** numerical_solution, const std::string& filename) 
{       
    double* phi = createVector(G->M);
    double** w = createMatrix(G->I + 1, G->J + 1);
    
    solveBIEMHelmholtz(G, phi, w); 

    int total_size = (G->I + 1) * (G->J + 1);
    #pragma omp parallel for num_threads(omp_get_max_threads())
    for (int idx = 0; idx < total_size; idx++) {
        int i = idx / (G->J + 1);  
        int j = idx % (G->J + 1);  
        numerical_solution[i][j] = 0;
        
        if (G->interior[i][j] & !G->near_bdry[i][j]) {
            // exact_solution[i][j] = g_MHelmholtz(G->Grid_x[i], G->Grid_y[j]);
            for (int k = 0; k < G->M; k++) {
                numerical_solution[i][j] += G->bdry_delta * KernelInteriorMHelmholtz(G->Grid_x[i], G->Grid_y[j], G->xy[k][0], G->xy[k][1], G->dxdy[k][0], G->dxdy[k][1]) * phi[k];
            }
            numerical_solution[i][j] += w[i][j];
        }
    }
    
    double** coor = createMatrix(20, 2);
    double* u = createVector(20);
    for (int i = 0; i < G->I + 1; i++) {
        for (int j = 0; j < G->J + 1; j++) {
            if (G->interior[i][j] & G->near_bdry[i][j]) {
                int count = 0;
                
                G->findGoodPoints(i, j, coor, u, numerical_solution, g_MHelmholtz, 0, count);
                G->findGoodPoints(i + 1, j, coor, u, numerical_solution, g_MHelmholtz, 1, count);
                G->findGoodPoints(i - 1, j, coor, u, numerical_solution, g_MHelmholtz, 1, count);
                G->findGoodPoints(i, j + 1, coor, u, numerical_solution, g_MHelmholtz, 1, count);
                G->findGoodPoints(i, j - 1, coor, u, numerical_solution, g_MHelmholtz, 1, count);
                G->findGoodPoints(i + 1, j + 1, coor, u, numerical_solution, g_MHelmholtz, 1, count);
                G->findGoodPoints(i - 1, j + 1, coor, u, numerical_solution, g_MHelmholtz, 1, count);
                G->findGoodPoints(i - 1, j - 1, coor, u, numerical_solution, g_MHelmholtz, 1, count);
                G->findGoodPoints(i + 1, j - 1, coor, u, numerical_solution, g_MHelmholtz, 1, count);

                numerical_solution[i][j] = computeNearBoundaryValuesMHelmholtz(coor, u, G->Grid_x[i], G->Grid_y[j], count);
                // exact_solution[i][j] = g_MHelmholtz(G->Grid_x[i], G->Grid_y[j]);

            }
        }
    }
    freeVector(u);
    freeMatrix(coor, 20);

    saveMatrixToFile(const_cast<const double**>(numerical_solution), G->I + 1, G->J + 1, filename);

    freeVector(phi);
    freeMatrix(w, G->I + 1);
}