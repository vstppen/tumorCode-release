#include <iostream>
#include <math.h>
#include <string.h>
#include <cstring>
#include <omp.h>    // OpenMP
#include <chrono>

#include "LinearAlgebra.h"
#include "CartesianGridAndControlPoints.h"
#include "PrimalDualAlgorithm.h"
#include "Const.h"
#include "CubicSpline.h"
#include "ExtractBoundaryPoints.h"
#include "NormalDerivative.h"
#include "Config.h"

#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseLU>

#include <vector>
#include <cmath> // For std::abs etc.

using Eigen::SparseMatrix;
using Eigen::Triplet;

double g_optimization(double x0, double y0) 
{   
    return 0;
}

double F_optimization(double x_min, double y_min, double dx, double dy, int I, int J, double** source, double x0, double y0)
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



void assembleSystemsSparse(CartesianGridAndControlPoints* G, double** source, SparseMatrix<double>& A, Eigen::VectorXd& b)
{
    int N = G->int_node_num;
    double dx = G->dx;
    double dy = G->dy;
    double dx2 = dx * dx;
    double dy2 = dy * dy;

    b = Eigen::VectorXd::Zero(N);
    A.resize(N, N);
    A.setZero();  

    for (int l = 0; l < N; l++) {
        int i = G->int_node_index[l][0];
        int j = G->int_node_index[l][1];

        double x = G->Grid_x[i];
        double y = G->Grid_y[j];

        b[l] = F_optimization(G->x_min, G->y_min, G->dx, G->dy, G->I, G->J, source, x, y);

        // ddx
        if (G->interior[i - 1][j] && G->interior[i + 1][j]) {
            A.coeffRef(l, l) += 2. / dx2;
            A.coeffRef(l, G->int_node_ij_to_l[i - 1][j]) = -1. / dx2;
            A.coeffRef(l, G->int_node_ij_to_l[i + 1][j]) = -1. / dx2;

        } else if (G->interior[i - 1][j] && !G->interior[i + 1][j]) {
            double x_b = G->findXIntersectByBisection(x, x + dx, y, 1e-10);
            double theta = (x_b - x) / dx;
            if (theta > 1e-6) {
                A.coeffRef(l, l) += 2. / (theta * dx2);
                A.coeffRef(l, G->int_node_ij_to_l[i - 1][j]) = - (2. * theta) / (theta * (1 + theta) * dx2);
                b[l] += 2. / (theta * (1 + theta) * dx2) * g_optimization(x_b, y);
            } else {
                A.coeffRef(l, l) = 1.;
                b[l] = g_optimization(x_b, y);
            }

        } else if (!G->interior[i - 1][j] && G->interior[i + 1][j]) {
            double x_b = G->findXIntersectByBisection(x - dx, x, y, 1e-10);
            double theta = (x - x_b) / dx;
            if (theta > 1e-6) {
                A.coeffRef(l, l) += 2. / (theta * dx2);
                A.coeffRef(l, G->int_node_ij_to_l[i + 1][j]) = - (2. * theta) / (theta * (1 + theta) * dx2);
                b[l] += 2. / (theta * (1 + theta) * dx2) * g_optimization(x_b, y);
            } else {
                A.coeffRef(l, l) = 1.;
                b[l] = g_optimization(x_b, y);
            }

        } else {
            double x_b_l = G->findXIntersectByBisection(x - dx, x, y, 1e-10);
            double theta_l = (x - x_b_l) / dx;
            double x_b_r = G->findXIntersectByBisection(x, x + dx, y, 1e-10);
            double theta_r = (x_b_r - x) / dx;
            if (theta_l < 1e-6) {
                A.coeffRef(l, l) = 1.;
                b[l] = g_optimization(x_b_l, y);
            } else if (theta_r < 1e-6) {
                A.coeffRef(l, l) = 1.;
                b[l] = g_optimization(x_b_r, y);
            } else {
                A.coeffRef(l, l) += 2. / (theta_l * theta_r * dx2);
                b[l] += 2. / (theta_l * (theta_l + theta_r) * dx2) * g_optimization(x_b_l, y);
                b[l] += 2. / (theta_r * (theta_l + theta_r) * dx2) * g_optimization(x_b_r, y);
            }
        }

        
        // ddy 
        if (G->interior[i][j - 1] && G->interior[i][j + 1]) {
            A.coeffRef(l, l) += 2. / dy2;
            A.coeffRef(l, G->int_node_ij_to_l[i][j - 1]) = -1. / dy2;
            A.coeffRef(l, G->int_node_ij_to_l[i][j + 1]) = -1. / dy2;

        } else if (G->interior[i][j - 1] && !G->interior[i][j + 1]) {
            double y_b = G->findYIntersectByBisection(x, y, y + dy, 1e-10);
            double theta = (y_b - y) / dy;
            if (theta > 1e-6) {
                A.coeffRef(l, l) += 2. / (theta * dy2);
                A.coeffRef(l, G->int_node_ij_to_l[i][j - 1]) = - (2. * theta) / (theta * (1 + theta) * dy2);
                b[l] += 2. / (theta * (1 + theta) * dy2) * g_optimization(x, y_b);
            } else {
                A.coeffRef(l, l) = 1.;
                b[l] = g_optimization(x, y_b);
            }

        } else if (!G->interior[i][j - 1] && G->interior[i][j + 1]) {
            double y_b = G->findYIntersectByBisection(x, y - dy, y, 1e-10);
            double theta = (y - y_b) / dy;
            if (theta > 1e-6) {
                A.coeffRef(l, l) += 2. / (theta * dy2);
                A.coeffRef(l, G->int_node_ij_to_l[i][j + 1]) = - (2. * theta) / (theta * (1 + theta) * dy2);
                b[l] += 2. / (theta * (1 + theta) * dy2) * g_optimization(x, y_b);
            } else {
                A.coeffRef(l, l) = 1.;
                b[l] = g_optimization(x, y_b);
            }

        } else {
            double y_b_l = G->findYIntersectByBisection(x, y - dy, y, 1e-10);
            double theta_l = (y - y_b_l) / dy;
            double y_b_r = G->findYIntersectByBisection(x, y, y + dy, 1e-10);
            double theta_r = (y_b_r - y) / dy;
            if (theta_l < 1e-6) {
                A.coeffRef(l, l) = 1.;
                b[l] = g_optimization(x, y_b_l);
            } else if (theta_r < 1e-6) {
                A.coeffRef(l, l) = 1.;
                b[l] = g_optimization(x, y_b_r);
            } else {
                A.coeffRef(l, l) += 2. / (theta_l * theta_r * dy2);
                b[l] += 2. / (theta_l * (theta_l + theta_r) * dy2) * g_optimization(x, y_b_l);
                b[l] += 2. / (theta_r * (theta_l + theta_r) * dy2) * g_optimization(x, y_b_r);
            }
        }
    }

    A.makeCompressed(); 
}


void selectSubmatrixAndVectorSparse(Eigen::SparseMatrix<double>& A, Eigen::VectorXd& b, int n, int* selected, 
                                    Eigen::SparseMatrix<double>& A_selected, Eigen::VectorXd& b_selected, int& n_selected, int*& index_map)
{
    n_selected = 0;
    for (int i = 0; i < n; ++i)
        if (selected[i]) ++n_selected;

    b_selected = Eigen::VectorXd::Zero(n_selected);
    index_map = new int[n];
    int idx = 0;
    for (int i = 0; i < n; ++i) {
        if (selected[i]) {
            index_map[i] = idx++;
        } else {
            index_map[i] = -1;
        }
    }

    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(A.nonZeros());  

    for (int k = 0; k < A.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(A, k); it; ++it) {
            int i = it.row();
            int j = it.col();
            if (selected[i] && selected[j]) {
                int row_new = index_map[i];
                int col_new = index_map[j];
                tripletList.emplace_back(row_new, col_new, it.value());
            }
        }
    }

    A_selected.resize(n_selected, n_selected);
    A_selected.setFromTriplets(tripletList.begin(), tripletList.end());
    A_selected.makeCompressed();  

    for (int i = 0; i < n; ++i) {
        if (selected[i]) {
            b_selected[index_map[i]] = b[i];
        }
    }
}


int PrimalDualAlgorithm(CartesianGridAndControlPoints* G, double** source, double** numerical_solution, double* numerical_v1, const std::string& filename_sol, const std::string& filename_v1, double**& int_bd_points)
{
    double itr_c = 1.;

    Eigen::SparseMatrix<double> A;
    Eigen::VectorXd b;
    assembleSystemsSparse(G, source, A, b);

    Eigen::VectorXd u = Eigen::VectorXd::Zero(G->int_node_num);
    Eigen::VectorXd lambda = Eigen::VectorXd::Zero(G->int_node_num);
    Eigen::VectorXd lambda_hat = Eigen::VectorXd::Zero(G->int_node_num);

    int* I = createIntVector(G->int_node_num);
    int* I_old = createIntVector(G->int_node_num);

    int n_selected = 0;
    int* index_map = nullptr;

    std::cout << "# {int_node} = " << G->int_node_num << std::endl;
    std::cout << "Running primal dual algorithm now... " << std::endl;

    for (int count = 0; count < G->int_node_num / 4; ++count) {

        for (int i = 0; i < G->int_node_num; ++i) {
            lambda_hat[i] = std::max(0.0, lambda[i] + itr_c * u[i]);
            I[i] = (std::abs(lambda_hat[i]) < 1e-9);
        }

        int difference = 0;
        for (int i = 0; i < G->int_node_num; ++i) {
            difference += (I[i] != I_old[i]);
        }

        // std::cout << "Running primal dual algorithm with iteration " << count << ", difference between I and I_old is " << difference << "." << std::endl;
        
        if (difference == 0) break;

        for (int i = 0; i < G->int_node_num; i++) {
            I_old[i] = I[i];
        }

        Eigen::SparseMatrix<double> A_selected;
        Eigen::VectorXd b_selected;
       
        if (index_map != nullptr) {
            delete[] index_map;
            index_map = nullptr;
        }

        selectSubmatrixAndVectorSparse(A, b, G->int_node_num, I, A_selected, b_selected, n_selected, index_map);

        A_selected.makeCompressed();

        Eigen::VectorXd u_I;

        Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::DiagonalPreconditioner<double>> solver;
        solver.compute(A_selected);
        if (solver.info() != Eigen::Success) {
            std::cerr << "BiCGSTAB decomposition failed!" << std::endl;
        }
        u_I = solver.solve(b_selected);
        if (solver.info() != Eigen::Success) {
            std::cerr << "BiCGSTAB solve failed!" << std::endl;
        }

        for (int i = 0; i < G->int_node_num; ++i) {
            u[i] = I[i] ? u_I[index_map[i]] : 0.0;
        }

        lambda = b - A * u;
        for (int i = 0; i < G->int_node_num; ++i) {
            lambda[i] = std::max(0.0, lambda[i]);
        }
    }

    // double** X_output = createMatrix(G->I + 1, G->J + 1);
    for (int i = 0; i < G->I + 1; i++) {
        for (int j = 0; j < G->J + 1; j++) {
            if (G->interior[i][j]) {
                numerical_solution[i][j] = - u[G->int_node_ij_to_l[i][j]];  // !!!
                // if (std::abs(numerical_solution[i][j]) > 1e-12) {
                //     X_output[i][j] = 1.;
                // } else {
                //     X_output[i][j] = 0.;
                // }
            } else {
                numerical_solution[i][j] = 0;
                // X_output[i][j] = -1.;
            }
        }
    }


    saveMatrixToFile(const_cast<const double**>(numerical_solution), G->I + 1, G->J + 1, filename_sol);
    
    for(int i = 0; i < G->n_ctrl; i++) {
        numerical_v1[i] = 0.;
    }
    
    getNormalDerivatives(G, g_optimization, numerical_solution, numerical_v1);

    for(int i = 0; i < G->n_ctrl; i++) {
        numerical_v1[i] = -numerical_v1[i];
    }
    
    saveVectorToFile(const_cast<const double*>(numerical_v1), G->n_ctrl, filename_v1);


    // Extract boundary points where u = 0
    Point* boundaryPoints = new Point[4 * ((G->I + 1) + (G->J + 1))];
    for (int i = 0; i < G->I + 1; i++) {
        for (int j = 0; j < G->J + 1; j++) {
            if (!G->interior[i][j] || G->irregular[i][j]) {
                numerical_solution[i][j] = -1.; 
            } 
        }
    }
    
    // saveMatrixToFile(const_cast<const double**>(numerical_solution), G->I + 1, G->J + 1, filename_sol);
    // saveMatrixToFile(const_cast<const double**>(X_output), G->I + 1, G->J + 1, "ttt.txt");
    // int numPoints = extractZeroLevelSet(X_output, G->I + 1, G->J + 1, G->Grid_x, G->Grid_y, boundaryPoints);
    int numPoints = extractZeroLevelSet(numerical_solution, G->I + 1, G->J + 1, G->Grid_x, G->Grid_y, boundaryPoints);
    std::cout << "# {bounadry points 0} = " << numPoints << std::endl;

    // Sort them counter-clockwise
    sortPointsCCW(boundaryPoints, numPoints);

    int_bd_points = new double*[numPoints];
    for (int i = 0; i < numPoints; ++i) {
        int_bd_points[i] = new double[2];
    }
    for (int i = 0; i < numPoints; i++) {
        int_bd_points[i][0] = boundaryPoints[i].x;
        int_bd_points[i][1] = boundaryPoints[i].y;
    }


    if (index_map != nullptr) {
        delete[] index_map;
        index_map = nullptr;
    }
    index_map = nullptr;
    freeIntVector(I_old);
    freeIntVector(I);
    delete[] boundaryPoints;

    return numPoints;
}