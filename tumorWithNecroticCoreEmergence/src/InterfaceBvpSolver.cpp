#include <iostream>
#include <math.h>
#include <string.h>
#include <cstring>
#include <omp.h>    // OpenMP
#include <chrono>
#include <gsl/gsl_sf_bessel.h>

#include "LinearAlgebra.h"
#include "CartesianGridAndControlPoints.h"
#include "InterfaceBvpSolver.h"
#include "FastAlgorithm.h"
#include "Const.h"
#include "LU.h"
#include "CubicSpline.h"
#include "GMRES.h"
#include "ExtractBoundaryData.h"
#include "Config.h"

double kappa_i = lambda * n_c;
double kappa_e = lambda;


// ============================================================
double U_e(double x, double y)
{
    return c_B;
}

double Jump_g(double x, double y)
{
    return 0.;
}

double Jump_j(double x, double y, double nx, double ny)
{
    return 0.;
}

double F_i(double x, double y)
{
    return 0.;
}

double F_e(double x, double y)
{
    return 0.;
}
// ============================================================


// ============================================================
void computeVolumeIntegral(CartesianGridAndControlPoints* G, double (*F)(double, double), double** w, double* bdry_w, double* bdry_wn, double kappa)
{
    int size[2] = {G->I + 1, G->J + 1};
    double h2 = G->dx * G->dy; 

    for (int i = 0; i < G->I + 1; i++) {
        for (int j = 0; j < G->J + 1; j++) {
            if (G->interior[i][j]) {
                w[i][j] = h2 * F(G->Grid_x[i], G->Grid_y[j]);
            } else {
                w[i][j] = 0.;
            }
        }
    }

    G->makeCorrection(F, w, kappa, true);

    FastMHelmholtzSolver2(w, size, kappa * h2); // fast Fourier transform based solver to compute the volume integral on the grid

    G->extractDirichletBoundaryData(w, F, bdry_w, G->M, kappa, true);
    G->extractNeumannBoundaryData(w, F, bdry_wn, G->M, kappa, true);
}

void computeSingleLayerPotential(CartesianGridAndControlPoints* G, double* phi_vec, double* psi_vec, double** w, double* bdry_w, double* bdry_wn, double kappa)
{
    int size[2] = {G->I + 1, G->J + 1};
    double h2 = G->dx * G->dy; 

    for (int i = 0; i < G->I + 1; i++) {
        for (int j = 0; j < G->J + 1; j++) {
            w[i][j] = 0.0;
        }
    }
    
    G->makeCorrection(phi_vec, psi_vec, G->M, w, kappa);

    FastMHelmholtzSolver2(w, size, kappa * h2); // fast Fourier transform based solver to compute the Single layer potential on the grid
    
    G->extractDirichletBoundaryData(w, phi_vec, psi_vec, bdry_w, G->M, kappa);
    G->extractNeumannBoundaryData(w, phi_vec, psi_vec, bdry_wn, G->M, kappa); 
}

void computeDoubleLayerPotential(CartesianGridAndControlPoints* G, double* phi_vec, double* psi_vec, double** w, double* bdry_w, double* bdry_wn, double kappa)
{
    int size[2] = {G->I + 1, G->J + 1};
    double h2 = G->dx * G->dy; 

    for (int i = 0; i < G->I + 1; i++) {
        for (int j = 0; j < G->J + 1; j++) {
            w[i][j] = 0.0;
        }
    }

    G->makeCorrection(phi_vec, psi_vec, G->M, w, kappa); 
    
    FastMHelmholtzSolver2(w, size, kappa * h2); // fast Fourier transform based solver to compute the Double layer potential on the grid
    
    G->extractDirichletBoundaryData(w, phi_vec, psi_vec, bdry_w, G->M, kappa);
    G->extractNeumannBoundaryData(w, phi_vec, psi_vec, bdry_wn, G->M, kappa); 
}
// ============================================================



// ============================================================
void Yi(CartesianGridAndControlPoints* G0, CartesianGridAndControlPoints* G1, double** w, double* bdry0_w, double* bdry0_wn)
{
    computeVolumeIntegral(G0, F_i, w, bdry0_w, bdry0_wn, kappa_i);
}

void Ye(CartesianGridAndControlPoints* G0, CartesianGridAndControlPoints* G1, double** w, double* bdry0_w, double* bdry0_wn, double* bdry1_w, double* bdry1_wn)
{
    double** w_part1 = createMatrix(G1->I + 1, G1->I + 1);
    double* bdry0_w_part1 = createVector(G1->M);
    double* bdry0_wn_part1 = createVector(G1->M);
    double* bdry1_w_part1 = createVector(G1->M);
    double* bdry1_wn_part1 = createVector(G1->M);

    double** w_part2 = createMatrix(G1->I + 1, G1->I + 1);
    double* bdry0_w_part2 = createVector(G1->M);
    double* bdry0_wn_part2 = createVector(G1->M);
    double* bdry1_w_part2 = createVector(G1->M);
    double* bdry1_wn_part2 = createVector(G1->M);
    

    computeVolumeIntegral(G1, F_e, w_part1, bdry1_w_part1, bdry1_wn_part1, kappa_e);
    extractBoundaryDataForContinousFunction(G0, w_part1, bdry0_w_part1, bdry0_wn_part1);

    computeVolumeIntegral(G0, F_e, w_part2, bdry0_w_part2, bdry0_wn_part2, kappa_e);
    extractBoundaryDataForContinousFunction(G1, w_part2, bdry1_w_part2, bdry1_wn_part2);

    for (int i = 0; i < G1->I + 1; i++) {
        for (int j = 0; j < G1->J + 1; j++) {
            w[i][j] = w_part1[i][j] - w_part2[i][j];
        }
    }

    for (int i = 0; i < G1->M; i++) {
        bdry0_w[i] = bdry0_w_part1[i] - bdry0_w_part2[i];
        bdry0_wn[i] = bdry0_wn_part1[i] - bdry0_wn_part2[i];
        bdry1_w[i] = bdry1_w_part1[i] - bdry1_w_part2[i];
        bdry1_wn[i] = bdry1_wn_part1[i] - bdry1_wn_part2[i];
    }


    freeVector(bdry1_wn_part2);
    freeVector(bdry1_w_part2);
    freeVector(bdry0_wn_part2);
    freeVector(bdry0_w_part2);
    freeMatrix(w_part2, G1->I + 1);
    freeVector(bdry1_wn_part1);
    freeVector(bdry1_w_part1);
    freeVector(bdry0_wn_part1);
    freeVector(bdry0_w_part1);
    freeMatrix(w_part1, G1->I + 1);
}

void Vi(CartesianGridAndControlPoints* G0, CartesianGridAndControlPoints* G1, double* psi_vec, double** w, double* bdry0_w, double* bdry0_wn)
{
    double* zeros = createVector(G1->M);
    
    computeSingleLayerPotential(G0, zeros, psi_vec, w, bdry0_w, bdry0_wn, kappa_i);

    freeVector(zeros);
}

void Ve(CartesianGridAndControlPoints* G0, CartesianGridAndControlPoints* G1, int flag, double* psi_vec, double** w, double* bdry0_w, double* bdry0_wn, double* bdry1_w, double* bdry1_wn)
{
    double* zeros = createVector(G1->M);

    if (flag == 0) {
        computeSingleLayerPotential(G0, zeros, psi_vec, w, bdry0_w, bdry0_wn, kappa_e);
        extractBoundaryDataForContinousFunction(G1, w, bdry1_w, bdry1_wn);
    } else if (flag == 1) {
        computeSingleLayerPotential(G1, zeros, psi_vec, w, bdry1_w, bdry1_wn, kappa_e);
        extractBoundaryDataForContinousFunction(G0, w, bdry0_w, bdry0_wn);
    } else {
        std::cout << "WARNING!!!" << std::endl;
    }
    
    freeVector(zeros);
}

void Wi(CartesianGridAndControlPoints* G0, CartesianGridAndControlPoints* G1, double* phi_vec, double** w, double* bdry0_w, double* bdry0_wn)
{
    double* zeros = createVector(G1->M);
    
    computeDoubleLayerPotential(G0, phi_vec, zeros, w, bdry0_w, bdry0_wn, kappa_i);

    freeVector(zeros);
}

void We(CartesianGridAndControlPoints* G0, CartesianGridAndControlPoints* G1, int flag, double* phi_vec, double** w, double* bdry0_w, double* bdry0_wn, double* bdry1_w, double* bdry1_wn)
{
    double* zeros = createVector(G1->M);

    if (flag == 0) {
        computeDoubleLayerPotential(G0, phi_vec, zeros, w, bdry0_w, bdry0_wn, kappa_e);
        extractBoundaryDataForContinousFunction(G1, w, bdry1_w, bdry1_wn);
    } else if (flag == 1) {
        computeDoubleLayerPotential(G1, phi_vec, zeros, w, bdry1_w, bdry1_wn, kappa_e);
        extractBoundaryDataForContinousFunction(G0, w, bdry0_w, bdry0_wn);
    } else {
        std::cout << "WARNING!!!" << std::endl;
    }
    
    freeVector(zeros);
}
// ============================================================


void solveInterfacePDE(CartesianGridAndControlPoints* G0, CartesianGridAndControlPoints* G1, double** numerical_solution, const std::string& filename)
{   
    kappa_i = lambda * n_c;
    kappa_e = lambda;
    double* BIE_rhs = createVector(3 * G1->M); 
    double* zeros = createVector(G1->M);
    double* USELESS = createVector(G1->M);

    double* g = createVector(G1->M);
    double* j = createVector(G1->M);
    double* h = createVector(G1->M);
    for (int i = 0; i < G1->M; i++) {
        g[i] = Jump_g(G0->xy[i][0], G0->xy[i][1]);
        j[i] = Jump_j(G0->xy[i][0], G0->xy[i][1], G0->nxny[i][0], G0->nxny[i][1]);
        h[i] = U_e(G1->xy[i][0], G1->xy[i][1]);
    }

    double** Y_i = createMatrix(G1->I + 1, G1->J + 1); //
    double** Y_e = createMatrix(G1->I + 1, G1->J + 1); //
    double* bdry0_Y_i = createVector(G1->M); //
    double* bdry0_Y_e = createVector(G1->M); //
    double* bdry0_dn_Y_i = createVector(G1->M); //
    double* bdry0_dn_Y_e = createVector(G1->M); //
    double* bdry1_Y_e = createVector(G1->M); //

    double** W_0_e_g = createMatrix(G1->I + 1, G1->J + 1); //
    double* bdry0_W_0_e_g = createVector(G1->M); //
    double* bdry0_dn_W_0_e_g = createVector(G1->M); //
    double* bdry1_W_0_e_g = createVector(G1->M); //

    double** V_0_i_j = createMatrix(G1->I + 1, G1->J + 1); //
    double* bdry0_V_0_i_j = createVector(G1->M); //
    double* bdry0_dn_V_0_i_j = createVector(G1->M); //

    double** W_1_e_h = createMatrix(G1->I + 1, G1->J + 1); //
    double* bdry0_W_1_e_h = createVector(G1->M); //
    double* bdry0_dn_W_1_e_h = createVector(G1->M); // 
    double* bdry1_W_1_e_h = createVector(G1->M); //
    
    Yi(G0, G1, Y_i, bdry0_Y_i, bdry0_dn_Y_i);
    Ye(G0, G1, Y_e, bdry0_Y_e, bdry0_dn_Y_e, bdry1_Y_e, USELESS);

    We(G0, G1, 0, g, W_0_e_g, bdry0_W_0_e_g, bdry0_dn_W_0_e_g, bdry1_W_0_e_g, USELESS);
    Vi(G0, G1, j, V_0_i_j, bdry0_V_0_i_j, bdry0_dn_V_0_i_j);
    We(G0, G1, 1, h, W_1_e_h, bdry0_W_1_e_h, bdry0_dn_W_1_e_h, bdry1_W_1_e_h, USELESS);

    for (int i = 0; i < G1->M; i++) {
        BIE_rhs[i] = bdry0_W_0_e_g[i] + bdry0_V_0_i_j[i] + bdry0_W_1_e_h[i] + bdry0_Y_i[i] + bdry0_Y_e[i];
        BIE_rhs[i + G1->M] = -j[i] + bdry0_dn_W_0_e_g[i] + bdry0_dn_V_0_i_j[i] + bdry0_dn_W_1_e_h[i] + bdry0_dn_Y_i[i] + bdry0_dn_Y_e[i];
        BIE_rhs[i + 2 * G1->M] = h[i] - bdry1_W_0_e_g[i] - bdry1_W_1_e_h[i] - bdry1_Y_e[i];
    }
        
    double* BIE_sol = createVector(3 * G1->M); 
    int max_m = G1->M;
    int max_itr_num = 2 * G1->M;
    double tol = 1.e-8;
    int itr_num = 0;
    bool status = GMRES(G0, G1, BIE_rhs, BIE_sol, max_m, max_itr_num, tol, itr_num);

    // printVector(BIE_sol, 3 * G1->M);

    double* phi = createVector(G1->M);
    double* psi = createVector(G1->M);
    double* Phi = createVector(G1->M);

    for (int i = 0; i < G1->M; i++) {
        phi[i] = BIE_sol[i];
        psi[i] = BIE_sol[i + G1->M];
        Phi[i] = BIE_sol[i + 2 * G1->M];
    }

    double** W_0_i_phi = createMatrix(G1->I + 1, G1->J + 1);
    double** W_0_e_phi = createMatrix(G1->I + 1, G1->J + 1);

    double** V_0_i_psi = createMatrix(G1->I + 1, G1->J + 1);
    double** V_0_e_psi = createMatrix(G1->I + 1, G1->J + 1);

    double** V_1_e_Phi = createMatrix(G1->I + 1, G1->J + 1);

    double* USELESS1 = createVector(G1->M);
    double* USELESS2 = createVector(G1->M);
    double* USELESS3 = createVector(G1->M);
    double* USELESS4 = createVector(G1->M);

    Wi(G0, G1, phi, W_0_i_phi, USELESS1, USELESS2);
    We(G0, G1, 0, phi, W_0_e_phi, USELESS1, USELESS2, USELESS3, USELESS4);

    Vi(G0, G1, psi, V_0_i_psi, USELESS1, USELESS2);
    Ve(G0, G1, 0, psi, V_0_e_psi, USELESS1, USELESS2, USELESS3, USELESS4);

    Ve(G0, G1, 1, Phi, V_1_e_Phi, USELESS1, USELESS2, USELESS3, USELESS4);

    
    for (int i = 0; i < G1->I + 1; i++) {
        for (int j = 0; j < G1->J + 1; j++) {
            if (G1->interior[i][j]) {
                numerical_solution[i][j] = W_0_i_phi[i][j] - W_0_e_phi[i][j] + W_0_e_g[i][j] + V_0_i_psi[i][j] - V_0_e_psi[i][j] + V_0_i_j[i][j] + V_1_e_Phi[i][j] + W_1_e_h[i][j] + Y_i[i][j] + Y_e[i][j];
            } else {
                numerical_solution[i][j] = 0.;
            }
        }
    }

    std::cout << "The number of discrete points on the boundary is " << G1->M << "." << std::endl;
    std::cout << "Times of GMRES iteration is " << itr_num << "." << std::endl;

    saveMatrixToFile(const_cast<const double**>(numerical_solution), G0->I + 1, G0->J + 1, filename);


    freeVector(USELESS4);
    freeVector(USELESS3);
    freeVector(USELESS2);
    freeVector(USELESS1);

    freeMatrix(V_1_e_Phi, G1->I + 1);

    freeMatrix(V_0_e_psi, G1->I + 1);
    freeMatrix(V_0_i_psi, G1->I + 1);

    freeMatrix(W_0_e_phi, G1->I + 1);
    freeMatrix(W_0_i_phi, G1->I + 1);

    freeVector(Phi);
    freeVector(psi);
    freeVector(phi);

    freeVector(BIE_sol);

    freeVector(bdry1_W_1_e_h);
    freeVector(bdry0_dn_W_1_e_h);
    freeVector(bdry0_W_1_e_h);
    freeMatrix(W_1_e_h, G1->I + 1);

    freeVector(bdry0_dn_V_0_i_j);
    freeVector(bdry0_V_0_i_j);
    freeMatrix(V_0_i_j, G1->I + 1);

    freeVector(bdry1_W_0_e_g);
    freeVector(bdry0_dn_W_0_e_g);
    freeVector(bdry0_W_0_e_g);
    freeMatrix(W_0_e_g, G1->I + 1);

    freeVector(bdry1_Y_e);
    freeVector(bdry0_dn_Y_e);
    freeVector(bdry0_Y_e);
    freeVector(bdry0_dn_Y_i);
    freeVector(bdry0_Y_i);
    freeMatrix(Y_e, G1->I + 1);
    freeMatrix(Y_i, G1->I + 1);

    freeVector(h);
    freeVector(j);
    freeVector(g);
    freeVector(USELESS);
    freeVector(zeros);
    freeVector(BIE_rhs);
}
