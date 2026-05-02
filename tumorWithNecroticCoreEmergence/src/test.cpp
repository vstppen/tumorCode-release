#include <iostream>
#include <string>
#include <chrono>
#include <omp.h>    // OpenMP
#include "LinearAlgebra.h"
#include "Domains.h"
#include "CartesianGridAndControlPoints.h"
#include "PrimalDualAlgorithm.h"
#include "InterfaceBvpSolver.h"
#include "FreeBoundaryModel.h"
#include "Config.h"

int main() 
{   
    // set the parameters
    lambda = 1;
    c_B = 10.;
    n_c = 1.e-3;
    model_G0 = 1.;
    c_bar = 0.5 * c_B;
    sqrt_lambda = sqrt(lambda);

    int n_ctrl_1 = 64;
    double** ctrl_points_1 = createMatrix(n_ctrl_1, 2);

    for (int i = 0; i < n_ctrl_1; i++) {
        double theta = i * 2 * M_PI / n_ctrl_1;
        double r = 2.2;
        // r = 2.2 + 0.2 * cos(4 * theta);
        ctrl_points_1[i][0] = 2.2 * cos(theta);
        ctrl_points_1[i][1] = 2.2 * sin(theta); 
    }

    double x_min = -5.;
    double x_max = 5.;
    double y_min = -5.;
    double y_max = 5.;
    int I = 512;
    int J = 512;
    int M = std::max(I, J);

    std::cout << "G_0 = " << model_G0 << ", c_B = " << c_B << ", lambda = " << lambda << ", n_c = " << n_c << "." << std::endl; 
    std::cout << "Cartesian grid size: " << I << " x " << J << "." << std::endl;
    std::cout << "Number of boundary samples: " << M << "." << std::endl; 
    
    auto start_ = std::chrono::steady_clock::now();

    testModel(x_min, x_max, y_min, y_max, I, J, ctrl_points_1, n_ctrl_1, M, 0.01, 0.2, 5.);

    auto end_ = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end_ - start_;
    std::cout << "Elapsed time: " << elapsed_seconds.count() << " seconds\n";

    freeMatrix(ctrl_points_1, n_ctrl_1);

    return 0;
}