#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>

#include "LinearAlgebra.h" 
#include "CartesianGridAndControlPoints.h"
#include "FastAlgorithm.h"
#include "Norm.h"
#include "Product.h"
#include "InterfaceBvpSolver.h"
#include "GMRES.h"

extern double kappa_i;
extern double kappa_e;

double makeQRdecomposition(double **H, double **R, double *cs, double *sn, double *b, int m, double atol)
{
    int m1 = m + 1;

    int ell = m - 1;
    for (int j = 0; j < m1; j++) {
        R[j][ell] = H[j][ell]; 
    }

    for (int k = 0; k < ell; k++) {
        double t = cs[k] * R[k][ell] + sn[k] * R[k + 1][ell];
        R[k + 1][ell] = - sn[k] * R[k][ell] + cs[k] * R[k + 1][ell]; 
        R[k][ell] = t;
    }

    if (fabs(R[ell + 1][ell]) > atol) {
        double x = R[ell][ell];
        double y = R[ell + 1][ell];
        double r = sqrt(x * x + y * y); 
        cs[ell] = x / r;
        sn[ell] = y / r; 

        int k = ell; 

        double t = cs[k] * R[k][ell] + sn[k] * R[k + 1][ell]; 
        R[k + 1][ell] = - sn[k] * R[k][ell] + cs[k] * R[k + 1][ell]; 
        R[k][ell] = t;

        t = cs[k] * b[k] + sn[k] * b[k + 1];
        b[k + 1] = - sn[k] * b[k] + cs[k] * b[k + 1]; 
        b[k] = t;
    }

    double res = b[m];

    return res; 
}

void computeMatrixVectorProduct(CartesianGridAndControlPoints* G0, CartesianGridAndControlPoints* G1, double* u, double* w) 
{
    double* phi = createVector(G1->M);
    double* psi = createVector(G1->M);
    double* Phi = createVector(G1->M);

    for (int i = 0; i < G1->M; i++) {
        phi[i] = u[i];
        psi[i] = u[i + G1->M];
        Phi[i] = u[i + 2 * G1->M];
    }

    double** W_0_i_phi = createMatrix(G1->I + 1, G1->J + 1);
    double** W_0_e_phi = createMatrix(G1->I + 1, G1->J + 1);
    double** V_0_i_psi = createMatrix(G1->I + 1, G1->J + 1);
    double** V_0_e_psi = createMatrix(G1->I + 1, G1->J + 1);
    double** V_1_e_Phi = createMatrix(G1->I + 1, G1->J + 1);

    double* bdry0_W_0_i_phi = createVector(G1->M);
    double* bdry0_W_0_e_phi = createVector(G1->M);
    double* bdry0_V_0_i_psi = createVector(G1->M);
    double* bdry0_V_0_e_psi = createVector(G1->M);
    double* bdry0_V_1_e_Phi = createVector(G1->M);

    double* bdry0_dn_W_0_i_phi = createVector(G1->M);
    double* bdry0_dn_W_0_e_phi = createVector(G1->M);
    double* bdry0_dn_V_0_i_psi = createVector(G1->M);
    double* bdry0_dn_V_0_e_psi = createVector(G1->M);
    double* bdry0_dn_V_1_e_Phi = createVector(G1->M);

    double* bdry1_W_0_e_phi = createVector(G1->M);
    double* bdry1_V_0_e_psi = createVector(G1->M);
    double* bdry1_V_1_e_Phi = createVector(G1->M);

    double* USELESS = createVector(G1->M);

    Wi(G0, G1, phi, W_0_i_phi, bdry0_W_0_i_phi, bdry0_dn_W_0_i_phi);
    We(G0, G1, 0, phi, W_0_e_phi, bdry0_W_0_e_phi, bdry0_dn_W_0_e_phi, bdry1_W_0_e_phi, USELESS);
    Vi(G0, G1, psi, V_0_i_psi, bdry0_V_0_i_psi, bdry0_dn_V_0_i_psi);
    Ve(G0, G1, 0, psi, V_0_e_psi, bdry0_V_0_e_psi, bdry0_dn_V_0_e_psi, bdry1_V_0_e_psi, USELESS);
    Ve(G0, G1, 1, Phi, V_1_e_Phi, bdry0_V_1_e_Phi, bdry0_dn_V_1_e_Phi, bdry1_V_1_e_Phi, USELESS);

    for (int i = 0; i < G1->M; i++) {
        w[i] = phi[i] - bdry0_W_0_i_phi[i] + bdry0_W_0_e_phi[i] - bdry0_V_0_i_psi[i] + bdry0_V_0_e_psi[i] - bdry0_V_1_e_Phi[i];
        w[i + G1->M] = psi[i] - bdry0_dn_W_0_i_phi[i] + bdry0_dn_W_0_e_phi[i] - bdry0_dn_V_0_i_psi[i] + bdry0_dn_V_0_e_psi[i] - bdry0_dn_V_1_e_Phi[i];
        w[i + 2 * G1->M] = - bdry1_W_0_e_phi[i] - bdry1_V_0_e_psi[i] + bdry1_V_1_e_Phi[i]; 
    }

    freeVector(USELESS);

    freeVector(bdry1_V_1_e_Phi);
    freeVector(bdry1_V_0_e_psi);
    freeVector(bdry1_W_0_e_phi);

    freeVector(bdry0_dn_V_1_e_Phi);
    freeVector(bdry0_dn_V_0_e_psi);
    freeVector(bdry0_dn_V_0_i_psi);
    freeVector(bdry0_dn_W_0_e_phi);
    freeVector(bdry0_dn_W_0_i_phi);

    freeVector(bdry0_V_1_e_Phi);
    freeVector(bdry0_V_0_e_psi);
    freeVector(bdry0_V_0_i_psi);
    freeVector(bdry0_W_0_e_phi);
    freeVector(bdry0_W_0_i_phi);

    freeMatrix(V_1_e_Phi, G1->I + 1);
    freeMatrix(V_0_e_psi, G1->I + 1);
    freeMatrix(V_0_i_psi, G1->I + 1);
    freeMatrix(W_0_e_phi, G1->I + 1);
    freeMatrix(W_0_i_phi, G1->I + 1);

    freeVector(Phi);
    freeVector(psi);
    freeVector(phi);
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
bool GMRES(CartesianGridAndControlPoints* G0, CartesianGridAndControlPoints* G1, const double *b, double *u, int max_m, int max_itr_num, double tol, int &itr_num)
{
    int n = 3 * G1->M;

    double atol = 1.0E-15;

    itr_num = 0; 

    double *r = new double[n];
    double *w = new double[n];
    for (int i = 0; i < n; i++) {
        r[i] = w[i] = 0.0;
    }

    computeMatrixVectorProduct(G0, G1, u, w);   // w = Au
    for (int i = 0; i < n; i++) {
        r[i] = b[i] - w[i]; 
    }   // r = b - Au
    // printVector(w, n);
    double denom = sqrt(n);
    double norm_r0 = computeVectorNorm(r, n);
    if ((norm_r0 / denom) < atol) {
        delete[] r; r = 0;
        delete[] w; w = 0;
        return true;
    }   // if r is small enough, u is the solution to Au = b

    int max_m1 = max_m + 1;
    int max_m2 = max_m + 2;

    double **V = new double*[max_m2]; 
    double **H = new double*[max_m2];
    double **R = new double*[max_m2]; 

    for (int i = 0; i < max_m2; i++) {
        V[i] = new double[n + 1]; 
        H[i] = new double[max_m1];
        R[i] = new double[max_m1];
    }

    for (int i = 0; i < max_m2; i++) {
        for (int j = 0; j <= n; j++) {
            V[i][j] = 0.0;
        }
        for (int j = 0; j < max_m1; j++) {
            H[i][j] = 0.0; 
            R[i][j] = 0.0; 
        }
    }

    double *cs = new double[max_m1]; 
    double *sn = new double[max_m1]; 
    for (int i = 0; i < max_m1; i++) {
        cs[i] = 1.0;
        sn[i] = 0.0; 
    }

    double *z = new double[max_m1]; 
    for (int i = 0; i < max_m1; i++) {
        z[i] = 0.0; 
    }

    int m = 0; 
    double beta = norm_r0; 
    double r_beta = 1.0 / beta;
    for (int j = 0; j < n; j++) {
        V[m][j] = r[j] * r_beta; 
    }

    double *c = new double [max_m2]; 
    for (int i = 1; i < max_m2; i++) {
        c[i] = 0.0; 
    }
    c[0] = beta; 

    int done = 0; 
    tol = tol * denom;

    double* temp = createVector(n);
    double* rr = createVector(n);
    while ((itr_num < max_itr_num) && (!done)) {
        // std::cout << "GMRES iteration with itr_num = " << itr_num ;

        computeMatrixVectorProduct(G0, G1, V[m], w);    // w_k = A v_k
        for (int j = 0; j <= m; j++) {
            H[j][m] = computeInnerProduct(V[j], w, n);
            for (int i = 0; i < n; i++) {
                w[i] -= V[j][i] * H[j][m];
            }
        }
        H[m + 1][m] = computeVectorNorm(w, n);

        if (fabs(H[m + 1][m]) > atol) {
            for (int i = 0; i < n; i++) {
                V[m + 1][i] = w[i] / H[m + 1][m];
            }
        } else {
            done = 1; 
        }

        double rz = makeQRdecomposition(H, R, cs, sn, c,  m + 1, atol); 

        itr_num++; 

        if (fabs(rz) < tol) {

            for (int i = m; i >= 0; i--) {
                double s = 0.0;
                for (int j = m; j > i; j--) {
                    s += R[i][j] * z[j]; 
                }
                z[i] = (c[i] - s) / R[i][i]; 
            }

            for (int i = 0; i < n; i++) {
                double sum = 0.0; 
                for (int j = 0; j <= m; j++) {
                    sum += V[j][i] * z[j];
                }
                u[i] += sum;
            }

            done = 1; 

        } else {

            if (m == max_m) {

                for (int i = m; i >= 0; i--) {
                    double s = 0.0;
                    for (int j = m; j > i; j--) {
                        s += R[i][j] * z[j]; 
                    }
                    z[i] = (c[i] - s) / R[i][i];
                }

                for (int i = 0; i < n; i++) {
                    double sum = 0.0; 
                    for (int j = 0; j <= m; j++) {
                        sum += V[j][i] * z[j];
                    }
                    u[i] += sum; 
                }

                computeMatrixVectorProduct(G0, G1, u, w); 
                for (int i = 0; i < n; i++) {
                    r[i] = b[i] - w[i]; 
                }
                norm_r0 = computeVectorNorm(r, n);
                if ((norm_r0 / sqrt(n)) < atol) {
                    done = 1; 
                } else {
                    m = 0; 
                    beta = norm_r0;
                    double r_beta = 1.0 / beta;
                    for (int j = 0; j < n; j++) {
                        V[m][j] = r[j] * r_beta; 
                    }
                    for (int i = 1; i < max_m2; i++) {
                        c[i] = 0.0; 
                    }
                    c[0] = beta; 
                }
            } else {
                m++; 
            }
        }

        computeMatrixVectorProduct(G0, G1, u, temp);   // w = Au
        for (int i = 0; i < n; i++) {
            rr[i] = b[i] - temp[i]; 
        }   // r = b - Au
        // std::cout << ", max-norm residual = " << computeMaxNorm(rr, n) << std::endl;
    }
    // freeVector(temp);
    // freeVector(rr);

    // computeMatrixVectorProduct(G, u, w);   // w = Au
    // for (int i = 0; i < n; i++) {
    //     r[i] = b[i] - w[i]; 
    // }   // r = b - Au
    // std::cout << "max-norm residual = " << computeMaxNorm(r, n) << std::endl;


    if (itr_num >= max_itr_num) {
        std::cout << "GMRES iteration failed." << std::endl;
        exit(EXIT_FAILURE);
        return false; 
    }

    delete[] cs; cs = 0; 
    delete[] sn; sn = 0; 

    for (int i = 0; i < max_m2; i++) {
        delete[] V[i]; V[i] = 0;
        delete[] H[i]; H[i] = 0; 
        delete[] R[i]; R[i] = 0;
    }
    delete[] H; H = 0;
    delete[] R; R = 0;
    delete[] V; V = 0;

    delete[] r; r = 0;
    delete[] w; w = 0;

    delete[] z; z = 0; 
    delete[] c; c = 0; 

    return true;
}
