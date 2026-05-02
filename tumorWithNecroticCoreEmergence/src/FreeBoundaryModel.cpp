#include <iostream>
#include <iomanip>
#include <string>
#include <chrono>
#include <cmath>
#include <vector>
#include <cstring> // for memcpy
#include <gsl/gsl_sf_bessel.h>

#include "LinearAlgebra.h"
#include "PrimalDualAlgorithm.h"
#include "InterfaceBvpSolver.h"
#include "CartesianGridAndControlPoints0.h"
#include "CartesianGridAndControlPoints.h"
#include "FreeBoundaryModel.h"
#include "CleanFiles.h"
#include "CubicSpline.h"
#include "BoundaryInclusionCheck.h"
#include "PoissonSolver.h"
#include "MHelmholtzSolver.h"
#include "Config.h"

#include <sstream>
#include <cstdio>
#include <array>
#include <fstream>


void evolveCtrlPointsStrategy1(double** ctrl_points, int& n_ctrl, double** nxny, double* v, double time_step)
{
    for (int i = 0; i < n_ctrl; i++) {
        if (v[i] * time_step < 10) {
            ctrl_points[i][0] += v[i] * time_step * nxny[i][0];
            ctrl_points[i][1] += v[i] * time_step * nxny[i][1];
        } else {
            std::cout << "Invalid velocity!" << std::endl;
            exit(EXIT_FAILURE); 
        }
    }
}

void evolveCtrlPointsStrategy2(double** ctrl_points, int& n_ctrl, double** nxny, double* v, double time_step)
{
    evolveCtrlPointsStrategy1(ctrl_points, n_ctrl, nxny, v, time_step);

    double* M_x = createVector(n_ctrl);
    double* alpha_x = createVector(n_ctrl);
    double* beta_x = createVector(n_ctrl);
    double* M_y = createVector(n_ctrl);
    double* alpha_y = createVector(n_ctrl);
    double* beta_y = createVector(n_ctrl);

    periodicCubicSplineCurveInterpolation(ctrl_points, M_x, alpha_x, beta_x, M_y, alpha_y, beta_y, n_ctrl);

    // std::cout << "radius at next time step: " << computeTotalLength(n_ctrl, M_x, alpha_x, beta_x, M_y, alpha_y, beta_y) / M_2PI << std::endl;
    SelectPointsUniformly(n_ctrl, M_x, alpha_x, beta_x, M_y, alpha_y, beta_y, ctrl_points);
    
    freeVector(M_x);
    freeVector(alpha_x);
    freeVector(beta_x);
    freeVector(M_y);
    freeVector(alpha_y);
    freeVector(beta_y);
}



bool reparameterizeCtrlPoints(double** ctrl_points_0, int n_ctrl_0) {
    const int N = n_ctrl_0;
    const int K_max = std::min(10, N / 2);     // Maximum number of Fourier modes
    double best_error = 1.e8;              
    std::vector<double> best_a_x, best_b_x, best_a_y, best_b_y;

    // Step 1: Build a closed version of the input points (append the first point to the end)
    int N_closed = N + 1;
    std::vector<double> x(N_closed), y(N_closed), t(N_closed);
    for (int i = 0; i < N; ++i) {
        x[i] = ctrl_points_0[i][0];
        y[i] = ctrl_points_0[i][1];
        t[i] = 2 * M_PI * i / N;
    }
    x[N] = ctrl_points_0[0][0];
    y[N] = ctrl_points_0[0][1];
    t[N] = 2 * M_PI;

    // Internal lambda: compute error for given K, and also calculate Fourier coefficients
    auto compute_error_for_K = [&](int K, std::vector<double>& a_x, std::vector<double>& b_x,
                                   std::vector<double>& a_y, std::vector<double>& b_y) -> double {
        auto compute_coeffs = [&](const std::vector<double>& f, std::vector<double>& a, std::vector<double>& b) {
            a.resize(K + 1);
            b.resize(K + 1);
            for (int k = 0; k <= K; ++k) {
                double a_k = 0.0, b_k = 0.0;
                for (int i = 0; i < N_closed; ++i) {
                    a_k += f[i] * cos(k * t[i]);
                    b_k += f[i] * sin(k * t[i]);
                }
                a[k] = (k == 0) ? (a_k / N_closed) : (2 * a_k / N_closed);
                b[k] = (k == 0) ? 0.0 : (2 * b_k / N_closed);
            }
        };

        compute_coeffs(x, a_x, b_x);
        compute_coeffs(y, a_y, b_y);

        int Ns = 1000;
        double max_error = 0.0;

        for (int i = 0; i < N; ++i) {
            double px = ctrl_points_0[i][0];
            double py = ctrl_points_0[i][1];
            double min_dist = 1e9;

            for (int j = 0; j < Ns; ++j) {
                double tj = 2 * M_PI * j / Ns;

                double X = a_x[0], Y = a_y[0];
                for (int k = 1; k <= K; ++k) {
                    X += a_x[k] * cos(k * tj) + b_x[k] * sin(k * tj);
                    Y += a_y[k] * cos(k * tj) + b_y[k] * sin(k * tj);
                }

                double dx = X - px;
                double dy = Y - py;
                double d = std::sqrt(dx * dx + dy * dy);
                if (d < min_dist) min_dist = d;
            }

            max_error = std::max(min_dist, max_error);
        }

        return max_error;
    };

    // Step 2: Automatically choose the best K 
    int best_K = 1;
    std::vector<double> a_x, b_x, a_y, b_y;
    for (int K = 1; K <= K_max; ++K) {
        double err = compute_error_for_K(K, a_x, b_x, a_y, b_y);
        if (err < best_error) {
            best_K = K;
            best_error = err;
            best_a_x = a_x;
            best_b_x = b_x;
            best_a_y = a_y;
            best_b_y = b_y;
        }
    }

    // Print the selected K to stdout
    std::cout << "Selected Fourier order K = " << best_K << ", error = " << best_error << std::endl;

    if (best_error > 0.1) {
        std::cout << "Failed to reconstruct the inner boundary!" << std::endl;
        return 0;
    }

    // Step 3: Reconstruct N resampled points on the fitted curve
    for (int i = 0; i < N; ++i) {
        double ti = 2 * M_PI * i / N;
        double xi = best_a_x[0], yi = best_a_y[0];
        for (int k = 1; k <= best_K; ++k) {
            xi += best_a_x[k] * cos(k * ti) + best_b_x[k] * sin(k * ti);
            yi += best_a_y[k] * cos(k * ti) + best_b_y[k] * sin(k * ti);
        }
        ctrl_points_0[i][0] = xi;
        ctrl_points_0[i][1] = yi;
    }

    return 1;
}


void recordTime(const std::string& filename_time_log, double current_time) {
    std::ofstream ofs(filename_time_log, std::ios::app);
    if (!ofs.is_open()) {
        std::cerr << "Cannot open time log file: " << filename_time_log << std::endl;
        return;
    }
    ofs << current_time << "\n";
}


void testModelWithoutNecroticCore(double x_min_, double x_max_, double y_min_, double y_max_, int I_, int J_, double** ctrl_points_, int n_ctrl_, int M_, 
                                  double std_time_step, double t_end, double& current_time, const std::string& filename_time_log, int& current_step)
{  
    double** c = createMatrix(I_ + 1, J_ + 1);
    double** p = createMatrix(I_ + 1, J_ + 1);
    double* v = createVector(n_ctrl_);
    
    while (current_time + 1e-14 < t_end)
    {
        double dt = std_time_step;

        if (current_time + dt > t_end) dt = t_end - current_time;
        
        std::cout << "\ncurrent step = " << current_step << ", there is no necrotic core!" << std::endl;
        std::cout << "Computing the cell population density and limit pressure at step " << current_step << ", time = " << current_time << ", current delta_t = " << dt << std::endl;

        std::string filename_ctrlPts_ = "./result/points/ctrl_points_1_" + std::to_string(current_step) + ".txt";
        std::string filename_ctrl_nxny_ = "./result/points/ctrl_nxny_1_" + std::to_string(current_step) + ".txt";
        std::string filename_bdryPts_ = "./result/points/bdry_points_1_" + std::to_string(current_step) + ".txt";
        std::string filename_c_ = "./result/solutions/c_" + std::to_string(current_step) + ".txt";
        std::string filename_p_ = "./result/solutions/p_" + std::to_string(current_step) + ".txt";
        std::string filename_v_ = "./result/solutions/v_" + std::to_string(current_step) + ".txt";
        std::string filename_ctrlPts_0_ = "./result/points/ctrl_points_0_" + std::to_string(current_step) + ".txt";
        std::string filename_ctrl_0_nxny_ = "./result/points/ctrl_nxny_0_" + std::to_string(current_step) + ".txt";
        std::string filename_bdryPts_0_ = "./result/points/bdry_points_0_" + std::to_string(current_step) + ".txt";
        std::ofstream(filename_ctrlPts_0_).close();
        std::ofstream(filename_ctrl_0_nxny_).close();
        std::ofstream(filename_bdryPts_0_).close();
        
        CartesianGridAndControlPoints0* G = new CartesianGridAndControlPoints0(x_min_, x_max_, y_min_, y_max_, I_, J_, ctrl_points_, n_ctrl_, M_, filename_ctrlPts_, filename_ctrl_nxny_, filename_bdryPts_, 1);

        solveModifiedHelmholtz(G, c, filename_c_);

        for (int i = 0; i < G->I + 1; i++) {
            for (int j = 0; j < G->J + 1; j++) {
                c[i][j] = - model_G0 * (c[i][j] - c_bar);
            }
        }

        solvePoisson(G, c, p, v, filename_p_, filename_v_);

        bool fflag = false;

        for (int i = 0; i < G->I + 1; i++) {
            for (int j = 0; j < G->J + 1; j++) {
                if (G->interior[i][j] && !G->near_bdry[i][j] && p[i][j] < 1e-8) {
                    fflag = true;
                }
            }
        }

        if (fflag)
        {
            std::cout << "Necrotic core detected at time = "
                      << current_time << "\n";
            delete G;
            break;
        }

        evolveCtrlPointsStrategy2(ctrl_points_, n_ctrl_, G->ctrlPts_nxny, v, dt);

        delete G;

        current_time += dt;
        current_step++;

        recordTime(filename_time_log, current_time);

        if (current_time >= t_end) break;
    }

    freeVector(v);
    freeMatrix(p, I_ + 1);
    freeMatrix(c, I_ + 1);
}


void testModelWithNecroticCore(double x_min_, double x_max_, double y_min_, double y_max_, int I_, int J_, double** ctrl_points_0, int n_ctrl_0, double** ctrl_points_1, int n_ctrl_1, int M_, 
                               double std_time_step, double t_end, double& current_time, const std::string& filename_time_log, int& current_step)
{   
    double** c = createMatrix(I_ + 1, J_ + 1);
    double** p = createMatrix(I_ + 1, J_ + 1);
    double* v1 = createVector(n_ctrl_1);

    int numPoints_int = 1;
    double** int_bd_points = nullptr;
    int numPoints_int_star = 1;
    double** int_bd_points_star = nullptr;
    int numPoints_int_star_star = 1;
    double** int_bd_points_star_star = nullptr;

    double* M_x_R = createVector(n_ctrl_1);
    double* alpha_x_R = createVector(n_ctrl_1);
    double* beta_x_R = createVector(n_ctrl_1);
    double* M_y_R = createVector(n_ctrl_1);
    double* alpha_y_R = createVector(n_ctrl_1);
    double* beta_y_R = createVector(n_ctrl_1);
    double* M_x_r = createVector(n_ctrl_0);
    double* alpha_x_r = createVector(n_ctrl_0);
    double* beta_x_r = createVector(n_ctrl_0);
    double* M_y_r = createVector(n_ctrl_0);
    double* alpha_y_r = createVector(n_ctrl_0);
    double* beta_y_r = createVector(n_ctrl_0);

    double** ctrl_points_0_star = createMatrix(n_ctrl_0, 2);
    double** ctrl_points_0_star_star = createMatrix(n_ctrl_0, 2);
    

    while (current_time + 1e-14 < t_end)
    {
        double dt = std_time_step;

        if (current_time + dt > t_end) dt = t_end - current_time;

        std::cout << std::endl;
        std::cout << "current time = " << current_time << ", there is a necrotic core!" << std::endl;

        if (!isBoundaryInside(ctrl_points_0, n_ctrl_0, ctrl_points_1, n_ctrl_1)) {
            std::cout << "Boundary 0 is NOT inside Boundary 1 !\n";
            exit(EXIT_FAILURE);
        } 

        periodicCubicSplineCurveInterpolation(ctrl_points_1, M_x_R, alpha_x_R, beta_x_R, M_y_R, alpha_y_R, beta_y_R, n_ctrl_1);
        periodicCubicSplineCurveInterpolation(ctrl_points_0, M_x_r, alpha_x_r, beta_x_r, M_y_r, alpha_y_r, beta_y_r, n_ctrl_0);
       
        // std::cout << "radius R at time step " << current_step << ": "<< std::fixed << std::setprecision(8) << computeTotalLength(n_ctrl_1, M_x_R, alpha_x_R, beta_x_R, M_y_R, alpha_y_R, beta_y_R) / M_2PI << std::defaultfloat << std::endl;
        // std::cout << "radius r at time step " << current_step << ": " << std::fixed << std::setprecision(8) << computeTotalLength(n_ctrl_0, M_x_r, alpha_x_r, beta_x_r, M_y_r, alpha_y_r, beta_y_r) / M_2PI << std::defaultfloat << std::endl;

        std::cout << "Computing the cell population density and limit pressure at step " << current_step << ", time = " << current_time << ", current delta_t = " << dt << std::endl;
        std::string filename_ctrlPts_0_ = "./result/points/ctrl_points_0_" + std::to_string(current_step) + ".txt";
        std::string filename_ctrl_0_nxny_ = "./result/points/ctrl_nxny_0_" + std::to_string(current_step) + ".txt";
        std::string filename_bdryPts_0_ = "./result/points/bdry_points_0_" + std::to_string(current_step) + ".txt";
        std::string filename_ctrlPts_1_ = "./result/points/ctrl_points_1_" + std::to_string(current_step) + ".txt";
        std::string filename_ctrl_1_nxny_ = "./result/points/ctrl_nxny_1_" + std::to_string(current_step) + ".txt";
        std::string filename_bdryPts_1_ = "./result/points/bdry_points_1_" + std::to_string(current_step) + ".txt";
        std::string filename_c_ = "./result/solutions/c_" + std::to_string(current_step) + ".txt";
        std::string filename_p_ = "./result/solutions/p_" + std::to_string(current_step) + ".txt";
        std::string filename_v_ = "./result/solutions/v_" + std::to_string(current_step) + ".txt";
        
        CartesianGridAndControlPoints* G0 = new CartesianGridAndControlPoints(x_min_, x_max_, y_min_, y_max_, I_, J_, ctrl_points_0, n_ctrl_0, M_, filename_ctrlPts_0_, filename_ctrl_0_nxny_, filename_bdryPts_0_, 1);
        CartesianGridAndControlPoints* G1 = new CartesianGridAndControlPoints(x_min_, x_max_, y_min_, y_max_, I_, J_, ctrl_points_1, n_ctrl_1, M_, filename_ctrlPts_1_, filename_ctrl_1_nxny_, filename_bdryPts_1_, 1);
        // std::cout << G0->int_node_num << "   !!!" << std::endl;
        // step 1
        solveInterfacePDE(G0, G1, c, filename_c_);

        for (int i = 0; i < G1->I + 1; i++) {
            for (int j = 0; j < G1->J + 1; j++) {
                c[i][j] = - model_G0 * (c[i][j] - c_bar);
            }
        }

        if (int_bd_points != nullptr) {
            for (int i = 0; i < numPoints_int; ++i)
                delete[] int_bd_points[i];
            delete[] int_bd_points;
            int_bd_points = nullptr;
        }
        numPoints_int = PrimalDualAlgorithm(G1, c, p, v1, filename_p_, filename_v_, int_bd_points);
        
        if (numPoints_int < 10) {
            exit(EXIT_FAILURE);
        }

        double* M_x_star = createVector(numPoints_int);
        double* alpha_x_star = createVector(numPoints_int);
        double* beta_x_star = createVector(numPoints_int);
        double* M_y_star = createVector(numPoints_int);
        double* alpha_y_star = createVector(numPoints_int);
        double* beta_y_star = createVector(numPoints_int);

        double* M_x_star_ = createVector(n_ctrl_0);
        double* alpha_x_star_ = createVector(n_ctrl_0);
        double* beta_x_star_ = createVector(n_ctrl_0);
        double* M_y_star_ = createVector(n_ctrl_0);
        double* alpha_y_star_ = createVector(n_ctrl_0);
        double* beta_y_star_ = createVector(n_ctrl_0);

        periodicCubicSplineCurveInterpolation(int_bd_points, M_x_star, alpha_x_star, beta_x_star, M_y_star, alpha_y_star, beta_y_star, numPoints_int);
        SelectPointsUniformly(numPoints_int, M_x_star, alpha_x_star, beta_x_star, M_y_star, alpha_y_star, beta_y_star, ctrl_points_0_star, n_ctrl_0);

        reparameterizeCtrlPoints(ctrl_points_0_star, n_ctrl_0);

        periodicCubicSplineCurveInterpolation(ctrl_points_0_star, M_x_star_, alpha_x_star_, beta_x_star_, M_y_star_, alpha_y_star_, beta_y_star_, n_ctrl_0);
        SelectPointsUniformly(n_ctrl_0, M_x_star_, alpha_x_star_, beta_x_star_, M_y_star_, alpha_y_star_, beta_y_star_, ctrl_points_0_star);
        // step 1


        // step 2
        CartesianGridAndControlPoints* G0_star = new CartesianGridAndControlPoints(x_min_, x_max_, y_min_, y_max_, I_, J_, ctrl_points_0_star, n_ctrl_0, M_, filename_ctrlPts_0_, filename_ctrl_0_nxny_, filename_bdryPts_0_, 0);
        
        solveInterfacePDE(G0_star, G1, c, filename_c_);
        
        for (int i = 0; i < G1->I + 1; i++) {
            for (int j = 0; j < G1->J + 1; j++) {
                c[i][j] = - model_G0 * (c[i][j] - c_bar);
            }
        }

        if (int_bd_points_star != nullptr) {
            for (int i = 0; i < numPoints_int_star; ++i)
                delete[] int_bd_points_star[i];
            delete[] int_bd_points_star;
            int_bd_points_star = nullptr;
        }
        numPoints_int_star = PrimalDualAlgorithm(G1, c, p, v1, filename_p_, filename_v_, int_bd_points_star);

        if (numPoints_int_star < 10) {
            exit(EXIT_FAILURE);
        }

        evolveCtrlPointsStrategy2(ctrl_points_1, n_ctrl_1, G1->ctrlPts_nxny, v1, dt);
        // step 2


        // step 3
        CartesianGridAndControlPoints* G1_star = new CartesianGridAndControlPoints(x_min_, x_max_, y_min_, y_max_, I_, J_, ctrl_points_1, n_ctrl_1, M_, filename_ctrlPts_0_, filename_ctrl_0_nxny_, filename_bdryPts_0_, 0);
        
        solveInterfacePDE(G0_star, G1_star, c, filename_c_);
        
        for (int i = 0; i < G1_star->I + 1; i++) {
            for (int j = 0; j < G1_star->J + 1; j++) {
                c[i][j] = - model_G0 * (c[i][j] - c_bar);
            }
        }

        if (int_bd_points_star_star != nullptr) {
            for (int i = 0; i < numPoints_int_star_star; ++i)
                delete[] int_bd_points_star_star[i];
            delete[] int_bd_points_star_star;
            int_bd_points_star_star = nullptr;
        }
        numPoints_int_star_star = PrimalDualAlgorithm(G1_star, c, p, v1, filename_p_, filename_v_, int_bd_points_star_star);

        if (numPoints_int_star_star < 10) {
            exit(EXIT_FAILURE);
        }

        double* M_x_star_star = createVector(numPoints_int_star_star);
        double* alpha_x_star_star = createVector(numPoints_int_star_star);
        double* beta_x_star_star = createVector(numPoints_int_star_star);
        double* M_y_star_star = createVector(numPoints_int_star_star);
        double* alpha_y_star_star = createVector(numPoints_int_star_star);
        double* beta_y_star_star = createVector(numPoints_int_star_star);

        double* M_x_star_star_ = createVector(n_ctrl_0);
        double* alpha_x_star_star_ = createVector(n_ctrl_0);
        double* beta_x_star_star_ = createVector(n_ctrl_0);
        double* M_y_star_star_ = createVector(n_ctrl_0);
        double* alpha_y_star_star_ = createVector(n_ctrl_0);
        double* beta_y_star_star_ = createVector(n_ctrl_0);

        periodicCubicSplineCurveInterpolation(int_bd_points_star_star, M_x_star_star, alpha_x_star_star, beta_x_star_star, M_y_star_star, alpha_y_star_star, beta_y_star_star, numPoints_int_star_star);
        SelectPointsUniformly(numPoints_int_star_star, M_x_star_star, alpha_x_star_star, beta_x_star_star, M_y_star_star, alpha_y_star_star, beta_y_star_star, ctrl_points_0_star_star, n_ctrl_0);

        reparameterizeCtrlPoints(ctrl_points_0_star_star, n_ctrl_0);

        periodicCubicSplineCurveInterpolation(ctrl_points_0_star_star, M_x_star_star_, alpha_x_star_star_, beta_x_star_star_, M_y_star_star_, alpha_y_star_star_, beta_y_star_star_, n_ctrl_0);
        SelectPointsUniformly(n_ctrl_0, M_x_star_star_, alpha_x_star_star_, beta_x_star_star_, M_y_star_star_, alpha_y_star_star_, beta_y_star_star_, ctrl_points_0_star_star);
        // step 3


        // step 4
        for (int i = 0; i < n_ctrl_0; i++) {
            ctrl_points_0[i][0] = 0.5 * ctrl_points_0_star[i][0] + 0.5 * ctrl_points_0_star_star[i][0];
            ctrl_points_0[i][1] = 0.5 * ctrl_points_0_star[i][1] + 0.5 * ctrl_points_0_star_star[i][1];
        }
        // step 4

        current_time += dt;
        current_step++;

        recordTime(filename_time_log, current_time);

        if (current_time >= t_end) break;


        delete G1;
        delete G0;
        delete G1_star;
        delete G0_star;

        freeVector(M_x_star_);
        freeVector(alpha_x_star_);
        freeVector(beta_x_star_);
        freeVector(M_y_star_);
        freeVector(alpha_y_star_);
        freeVector(beta_y_star_);

        freeVector(M_x_star);
        freeVector(alpha_x_star);
        freeVector(beta_x_star);
        freeVector(M_y_star);
        freeVector(alpha_y_star);
        freeVector(beta_y_star);

        freeVector(M_x_star_star_);
        freeVector(alpha_x_star_star_);
        freeVector(beta_x_star_star_);
        freeVector(M_y_star_star_);
        freeVector(alpha_y_star_star_);
        freeVector(beta_y_star_star_);     
        
        freeVector(M_x_star_star);
        freeVector(alpha_x_star_star);
        freeVector(beta_x_star_star);
        freeVector(M_y_star_star);
        freeVector(alpha_y_star_star);
        freeVector(beta_y_star_star);

    }

    if (int_bd_points != nullptr) {
        for (int i = 0; i < numPoints_int; ++i)
            delete[] int_bd_points[i];
        delete[] int_bd_points;
        int_bd_points = nullptr;
    }

    if (int_bd_points_star != nullptr) {
        for (int i = 0; i < numPoints_int_star; ++i)
            delete[] int_bd_points_star[i];
        delete[] int_bd_points_star;
        int_bd_points_star = nullptr;
    }

    if (int_bd_points_star_star != nullptr) {
        for (int i = 0; i < numPoints_int_star_star; ++i)
            delete[] int_bd_points_star_star[i];
        delete[] int_bd_points_star_star;
        int_bd_points_star_star = nullptr;
    }

    freeVector(v1);
    freeMatrix(p, I_ + 1);
    freeMatrix(c, I_ + 1);

    freeMatrix(ctrl_points_0_star, n_ctrl_0);
    freeMatrix(ctrl_points_0_star_star, n_ctrl_0);

    freeVector(M_x_R);
    freeVector(alpha_x_R);
    freeVector(beta_x_R);
    freeVector(M_y_R);
    freeVector(alpha_y_R);
    freeVector(beta_y_R);
    freeVector(M_x_r);
    freeVector(alpha_x_r);
    freeVector(beta_x_r);
    freeVector(M_y_r);
    freeVector(alpha_y_r);
    freeVector(beta_y_r);
}


void testModel(double x_min_, double x_max_, double y_min_, double y_max_, int I_, int J_, double** ctrl_points_1, int n_ctrl_1, int M_, double std_time_step, double t_end, double d)
{   
    clear_directory("./result/points");
    clear_directory("./result/solutions");

    std::string filename_time_log = "./result/time_log.txt";

    std::ofstream os(filename_time_log, std::ios::trunc);
    os << "0\n";
    os.close();

    int current_step = 0;
    double current_time = 0.; 
    
    // Stage 1 ============================================================================
    std::cout << "\nBegin stage 1! " << std::endl;
    testModelWithoutNecroticCore(x_min_, x_max_, y_min_, y_max_, I_, J_, ctrl_points_1, n_ctrl_1, M_, std_time_step, t_end, current_time, filename_time_log , current_step);
    // ============================================================================ Stage 1


    // Stage 2 ============================================================================
    std::cout << "\nBegin stage 2! " << std::endl;

    int n_ctrl_0 = n_ctrl_1;
    double** ctrl_points_0 = createMatrix(n_ctrl_0, 2);
    bool flag = 0;
    while (flag == 0 && current_time + 1e-14 < t_end) {
        double time_step = std_time_step / d;
        if (current_time + time_step > t_end) time_step = t_end - current_time;

        std::cout << "\nComputing the cell population density and limit pressure at step " << current_step << ", time = " << current_time << ", current delta_t = " << time_step << std::endl;
        std::string filename_ctrlPts_ = "./result/points/ctrl_points_1_" + std::to_string(current_step) + ".txt";
        std::string filename_ctrl_nxny_ = "./result/points/ctrl_nxny_1_" + std::to_string(current_step) + ".txt";
        std::string filename_bdryPts_ = "./result/points/bdry_points_1_" + std::to_string(current_step) + ".txt";
        std::string filename_c_ = "./result/solutions/c_" + std::to_string(current_step) + ".txt";
        std::string filename_p_ = "./result/solutions/p_" + std::to_string(current_step) + ".txt";
        std::string filename_v_ = "./result/solutions/v_" + std::to_string(current_step) + ".txt";

        std::string filename_ctrlPts_0_ = "./result/points/ctrl_points_0_" + std::to_string(current_step) + ".txt";
        std::string filename_ctrl_0_nxny_ = "./result/points/ctrl_nxny_0_" + std::to_string(current_step) + ".txt";
        std::string filename_bdryPts_0_ = "./result/points/bdry_points_0_" + std::to_string(current_step) + ".txt";
        std::ofstream(filename_ctrlPts_0_).close();
        std::ofstream(filename_ctrl_0_nxny_).close();
        std::ofstream(filename_bdryPts_0_).close();
        
        CartesianGridAndControlPoints0* G = new CartesianGridAndControlPoints0(x_min_, x_max_, y_min_, y_max_, I_, J_, ctrl_points_1, n_ctrl_1, M_, filename_ctrlPts_, filename_ctrl_nxny_, filename_bdryPts_, 1);
        CartesianGridAndControlPoints* G1 = new CartesianGridAndControlPoints(x_min_, x_max_, y_min_, y_max_, I_, J_, ctrl_points_1, n_ctrl_1, M_, filename_ctrlPts_, filename_ctrl_nxny_, filename_bdryPts_, 1);
        double** c = createMatrix(I_ + 1, J_ + 1);
        double** p = createMatrix(I_ + 1, J_ + 1);
        int numPoints_int = 1;
        double** int_bd_points = nullptr;
        int numPoints_int_star = 1;
        double* v1 = createVector(n_ctrl_1);
        
        
        solveModifiedHelmholtz(G, c, filename_c_);


        for (int i = 0; i < G->I + 1; i++) {
            for (int j = 0; j < G->J + 1; j++) {
                c[i][j] = - model_G0 * (c[i][j] - c_bar);
            }
        }

        if (int_bd_points != nullptr) {
            for (int i = 0; i < numPoints_int; ++i)
                delete[] int_bd_points[i];
            delete[] int_bd_points;
            int_bd_points = nullptr;
        }
        numPoints_int = PrimalDualAlgorithm(G1, c, p, v1, filename_p_, filename_v_, int_bd_points);
        
        if (numPoints_int > 10) {
            double* M_x = createVector(numPoints_int);
            double* alpha_x = createVector(numPoints_int);
            double* beta_x = createVector(numPoints_int);
            double* M_y = createVector(numPoints_int);
            double* alpha_y = createVector(numPoints_int);
            double* beta_y = createVector(numPoints_int);

            double* M_x_ = createVector(n_ctrl_0);
            double* alpha_x_ = createVector(n_ctrl_0);
            double* beta_x_ = createVector(n_ctrl_0);
            double* M_y_ = createVector(n_ctrl_0);
            double* alpha_y_ = createVector(n_ctrl_0);
            double* beta_y_ = createVector(n_ctrl_0);

            periodicCubicSplineCurveInterpolation(int_bd_points, M_x, alpha_x, beta_x, M_y, alpha_y, beta_y, numPoints_int);
            SelectPointsUniformly(numPoints_int, M_x, alpha_x, beta_x, M_y, alpha_y, beta_y, ctrl_points_0, n_ctrl_0);

            flag = reparameterizeCtrlPoints(ctrl_points_0, n_ctrl_0);

            if (flag == 1) {
                periodicCubicSplineCurveInterpolation(ctrl_points_0, M_x_, alpha_x_, beta_x_, M_y_, alpha_y_, beta_y_, n_ctrl_0);
                SelectPointsUniformly(n_ctrl_0, M_x_, alpha_x_, beta_x_, M_y_, alpha_y_, beta_y_, ctrl_points_0);
            } 
            
            freeVector(M_x_);
            freeVector(alpha_x_);
            freeVector(beta_x_);
            freeVector(M_y_);
            freeVector(alpha_y_);
            freeVector(beta_y_);
            freeVector(M_x);
            freeVector(alpha_x);
            freeVector(beta_x);
            freeVector(M_y);
            freeVector(alpha_y);
            freeVector(beta_y);
        }
        
        evolveCtrlPointsStrategy2(ctrl_points_1, n_ctrl_1, G1->ctrlPts_nxny, v1, time_step);

        if (int_bd_points != nullptr) {
            for (int i = 0; i < numPoints_int; ++i)
                delete[] int_bd_points[i];
            delete[] int_bd_points;
            int_bd_points = nullptr;
        }
        freeVector(v1);    
        freeMatrix(p, I_ + 1);
        freeMatrix(c, I_ + 1);
        delete G1;
        delete G;

        current_time += time_step;
        current_step ++;

        recordTime(filename_time_log, current_time);

        if (current_time >= t_end) break;
    }
    
    // ============================================================================ Stage 2

    // saveMatrixToFile(const_cast<const double**>(ctrl_points_0), n_ctrl_0, 2, "./result/000.txt");
    // saveMatrixToFile(const_cast<const double**>(ctrl_points_1), n_ctrl_1, 2, "./result/111.txt");
    // Stage 3 ============================================================================
    std::cout << "\nBegin stage 3! " << std::endl;
    testModelWithNecroticCore(x_min_, x_max_, y_min_, y_max_, I_, J_, ctrl_points_0, n_ctrl_0, ctrl_points_1, n_ctrl_1, M_, std_time_step, t_end, current_time, filename_time_log, current_step);
    
    std::cout << "\nSaving Data at t_end." << std::endl;
    std::string filename_ctrlPts_0_ = "./result/points/ctrl_points_0_" + std::to_string(current_step) + ".txt";
    std::string filename_ctrl_0_nxny_ = "./result/points/ctrl_nxny_0_" + std::to_string(current_step) + ".txt";
    std::string filename_bdryPts_0_ = "./result/points/bdry_points_0_" + std::to_string(current_step) + ".txt";
    std::string filename_ctrlPts_1_ = "./result/points/ctrl_points_1_" + std::to_string(current_step) + ".txt";
    std::string filename_ctrl_1_nxny_ = "./result/points/ctrl_nxny_1_" + std::to_string(current_step) + ".txt";
    std::string filename_bdryPts_1_ = "./result/points/bdry_points_1_" + std::to_string(current_step) + ".txt";
    std::string filename_c_ = "./result/solutions/c_" + std::to_string(current_step) + ".txt";
    std::string filename_p_ = "./result/solutions/p_" + std::to_string(current_step) + ".txt";
    std::string filename_v_ = "./result/solutions/v_" + std::to_string(current_step) + ".txt";
    CartesianGridAndControlPoints* G0 = new CartesianGridAndControlPoints(x_min_, x_max_, y_min_, y_max_, I_, J_, ctrl_points_0, n_ctrl_0, M_, filename_ctrlPts_0_, filename_ctrl_0_nxny_, filename_bdryPts_0_, 1);
    CartesianGridAndControlPoints* G1 = new CartesianGridAndControlPoints(x_min_, x_max_, y_min_, y_max_, I_, J_, ctrl_points_1, n_ctrl_1, M_, filename_ctrlPts_1_, filename_ctrl_1_nxny_, filename_bdryPts_1_, 1);
    // ============================================================================ Stage 3


    freeMatrix(ctrl_points_0, n_ctrl_0);
}
