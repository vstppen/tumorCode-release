#include <iostream>
#include <string>
#include <chrono>
#include <omp.h>    // OpenMP
#include "LinearAlgebra.h"
#include "Domains.h"
#include "CartesianGridAndControlPoints.h"
#include "PoissonSolver.h"
#include "MHelmholtzSolver.h"
#include "FreeBoundaryModel.h"

int main() 
{   
    // set the parameters
    lambda = 1.;
    update_sqrt_lambda();
    c_B = 10;
    G0 = 1;

    int domain_case = 1;
    int n_ctrl;
    double** ctrl_points;

    double x_min = -10;
    double x_max = 10;
    double y_min = -10;
    double y_max = 10;
    int I = 256;
    int J = 256;
    int M = 4 * std::max(I, J);

    switch (domain_case) {
        case 1:     // ellipse or disk
            n_ctrl = 32;
            ctrl_points = createMatrix(n_ctrl, 2);
            for (int i = 0; i < n_ctrl; i++) {
                double theta = i * 2 * M_PI / n_ctrl;
                ctrl_points[i][0] = 2.3 * cos(theta);
                ctrl_points[i][1] = 1.1 * sin(theta); 
            }
            break;
        case 2:     // kidney
            n_ctrl = n_kidney;
            ctrl_points = createMatrix(n_ctrl, 2);
            for (int i = 0; i < n_ctrl; i++) {
                ctrl_points[i][0] = kidney[i][0];
                ctrl_points[i][1] = kidney[i][1]; 
            }
            break;
        case 3:     // flower
            n_ctrl = n_flower;
            ctrl_points = createMatrix(n_ctrl, 2);
            for (int i = 0; i < n_ctrl; i++) {
                ctrl_points[i][0] = flower[i][0];
                ctrl_points[i][1] = flower[i][1]; 
            }
            break;
        case 4:     // perturbed disk
            n_ctrl = 128;
            ctrl_points = createMatrix(n_ctrl, 2);
            for (int i = 0; i < n_ctrl; i++) {
                double theta = i * 2 * M_PI / n_ctrl;
                double r = 0.8 + 0.02 * cos(12 * theta);
                ctrl_points[i][0] = r * cos(theta);
                ctrl_points[i][1] = r * sin(theta); 
            }
            break;
        default:
            printf("Invalid domain case!\n");
            break;
    }

    std::cout << "Begin free boundary model for tumor pressure coupled with in vitro nutrient regime now! " << std::endl;
    std::cout << "The initial domain is defined by 'domain case' " << domain_case << "." << std::endl;
    std::cout << "G_0 = " << G0 << ", c_B = " << c_B << ", lambda = " << lambda << "." << std::endl; 
    std::cout << "Cartesian grid size: " << I << " x " << J << "." << std::endl;
    std::cout << "Number of boundary samples: " << M << "." << std::endl << std::endl; 
    
    auto start_ = std::chrono::steady_clock::now();

    testModel(x_min, x_max, y_min, y_max, I, J, ctrl_points, n_ctrl, M, 0.01, 101);

    auto end_ = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end_ - start_;
    std::cout << "Elapsed time: " << elapsed_seconds.count() << " seconds\n";

    freeMatrix(ctrl_points, n_ctrl);

}