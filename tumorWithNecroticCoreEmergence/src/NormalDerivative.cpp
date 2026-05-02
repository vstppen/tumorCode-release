#include "LinearAlgebra.h"
#include "CartesianGridAndControlPoints.h"
#include "CartesianGridAndControlPoints0.h"
// #include <omp.h>    // OpenMP
#include <Const.h>
#include <LU.h>
#include <algorithm>
#include "NormalDerivative.h"
#include <vector>


void findInteriorPoints0(CartesianGridAndControlPoints0* G, double (*g)(double, double), int bdry_idx, double** u, double& x, double& y, double** coord, double* grid_u, int& count) 
{
    double xi = G->ctrl_points[bdry_idx][0];
    double eta = G->ctrl_points[bdry_idx][1];
    
    int i = static_cast<int>((xi - G->Grid_x[0]) / G->dx + 0.5); 
    int j = static_cast<int>((eta - G->Grid_y[0]) / G->dy + 0.5); 

    if (i >= 0 && i < G->I + 1 && j >= 0 && j < G->J + 1 && G->interior[i][j]) {
        i = i; j = j;
    } else if (i - 1 >= 0 && i - 1 < G->I + 1 && j >= 0 && j < G->J + 1 && G->interior[i - 1][j]) {
        i = i - 1; j = j;
    } else if (i + 1 >= 0 && i + 1 < G->I + 1 && j >= 0 && j < G->J + 1 && G->interior[i + 1][j]) {
        i = i + 1; j = j;
    } else if (i >= 0 && i < G->I + 1 && j - 1 >= 0 && j - 1 < G->J + 1 && G->interior[i][j - 1]) {
        i = i; j = j - 1;
    } else if (i >= 0 && i < G->I + 1 && j + 1 >= 0 && j + 1 < G->J + 1 && G->interior[i][j + 1]) {
        i = i; j = j + 1;
    } else {
        std::cout << "Warning about finding interior node!" << std::endl;
    }

    coord[0][0] = 0.0;  coord[0][1] = 0.0;  grid_u[0] = u[i][j];      // i, j
    // coord[0][0] = G->Grid_x[i];  coord[0][1] = G->Grid_y[j];  grid_u[0] = u[i][j];

    double xi0 = G->Grid_x[i]; 
    double eta0 = G->Grid_y[j];

    x = (xi - xi0) / G->dx;
    y = (eta - eta0) / G->dy;
    // x = xi;
    // y = eta;

    coord[1][0] = x;  coord[1][1] = y;  grid_u[1] = g(xi, eta);
    
    // std::cout << x << " " << y << " " << (xi - xi0) << " " << (eta - eta0) << " " << xi << " " << eta << std::endl;

    int directions[24][2] = {{-1, 0}, {1, 0}, {0, -1}, {0, 1}, {-1, -1}, {1, 1}, {-1, 1}, {1, -1}, 
                            {-2, 0}, {2, 0}, {0, -2}, {0, 2}, {-2, 1}, {-1, 2}, {-2, -1}, {-1, -2}, 
                            {1, -2}, {2, -1}, {2, 1}, {1, 2}, {-2, -2}, {-2, 2}, {2, -2}, {2, 2}}; 
    // std::cout << "flag1" << std::endl;
    count = 2;  
    for (int k = 0; k < 24; ++k) {
        int di = directions[k][0]; 
        int dj = directions[k][1];

        int ni = i + di;
        int nj = j + dj;
        // std::cout << ni << " " << nj << " " << std::endl;
        if (ni >= 0 && ni < G->I + 1 && nj >= 0 && nj < G->J + 1 && G->interior[ni][nj]) {
            coord[count][0] = di / 1.; 
            coord[count][1] = dj / 1.; 
            // coord[count][0] = G->Grid_x[ni]; 
            // coord[count][1] = G->Grid_y[nj]; 
            grid_u[count] = u[ni][nj];
            count++;
            // if (count >= 6) {
            //     break;
            // }
        }
        // std::cout << ni << " " << nj << " " << std::endl;
    }
    // std::cout << "flag2" << std::endl;
    // if (count < 6) {
    //     std::cout << "Warning: Less than 6 interior points found!" << std::endl;
    // }
    
}

void computeBoundaryValues0(double x, double y, double** coord, double* grid_u, int count, double& u, double& ux, double& uy) {
    // Arrays to store distances and indices of the points
    double* distances = createVector(count);
    int* indices = createIntVector(count);

    // Calculate distances from (x, y) to each point in coor
    for (int i = 0; i < count; i++) {
        distances[i] = (x - coord[i][0]) * (x - coord[i][0]) + (y - coord[i][1]) * (y - coord[i][1]);
        indices[i] = i;  // Store original index of each point
    }

    // Sort points by distance to (x, y)
    std::sort(indices, indices + count, [&](int i1, int i2) {
        return distances[i1] < distances[i2];
    });

    if (count > 10) {
        // std::cout << count << std::endl;
        count = 10;
    }

    double** A = createMatrix(count, 6);
    double** AT = createMatrix(6, count);
    double** ATA = createMatrix(6, 6);
    double* b = createVector(count);
    double* ATb = createVector(6);
    double* sol = createVector(6);
    for (int i = 0; i < count; i++) {
        int idx = indices[i];  
        A[i][0] = coord[i][0] * coord[i][0];
        A[i][1] = coord[i][0] * coord[i][1];
        A[i][2] = coord[i][1] * coord[i][1];
        A[i][3] = coord[i][0];
        A[i][4] = coord[i][1];
        A[i][5] = 1;
        b[i] = grid_u[i];
    }
 
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < count; j++) {
            AT[i][j] = A[j][i];
        }
    }

    MatrixMulMatrix_L(AT, A, ATA, 6, count, 6);
    MatrixMulVector_L(AT, b, ATb, 6, count);
    SolveLinearSystemByCG(ATA, ATb, sol, 6, 0);
    u = sol[0] * x * x + sol[1] * x * y + sol[2] * y * y + sol[3] * x + sol[4] * y + sol[5];
    ux = 2 * sol[0] * x + sol[1] * y + sol[3];
    uy = sol[1] * x + 2 * sol[2] * y + sol[4];
    
    freeMatrix(A, count);
    freeMatrix(AT, 6);
    freeMatrix(ATA, 6);
    freeVector(b);
    freeVector(ATb);
    freeVector(sol);

}

void extractboundaryData0(CartesianGridAndControlPoints0* G, double (*g)(double, double), double** w, double* ux, double* uy) 
{
    double x, y;
    double** coord = createMatrix(24, 2);
    double* grid_u = createVector(24);
    int count;
    double temp;

    for (int i = 0; i < G->n_ctrl; i++) {
        count = 0;
        findInteriorPoints0(G, g, i, w, x, y, coord, grid_u, count);

        // std::cout << i << std::endl;

        computeBoundaryValues0(x, y, coord, grid_u, count, temp, ux[i], uy[i]);
        ux[i] /= G->dx;
        uy[i] /= G->dy;
        // if (i == 103) {
        //     std::cout << std::endl;
        //     std::cout << coord[0][0] << " " << coord[0][1] << std::endl;
        //     std::cout << coord[1][0] << " " << coord[1][1] << std::endl;
        //     std::cout << coord[2][0] << " " << coord[2][1] << std::endl;
        //     std::cout << coord[3][0] << " " << coord[3][1] << std::endl;
        //     std::cout << coord[4][0] << " " << coord[4][1] << std::endl;
        //     std::cout << coord[5][0] << " " << coord[5][1] << std::endl;
        //     std::cout << grid_u[0] << std::endl;
        //     std::cout << grid_u[1] << std::endl;
        //     std::cout << grid_u[2] << std::endl;
        //     std::cout << grid_u[3] << std::endl;
        //     std::cout << grid_u[4] << std::endl;
        //     std::cout << grid_u[5] << std::endl;
        //     std::cout << x << "  " << y << " " << v[i] << std::endl;
        // }
    }

    freeVector(grid_u);
    freeMatrix(coord, 24);
}

void getNormalDerivatives0(CartesianGridAndControlPoints0* G, double (*g)(double, double), double** u, double* v)
{
    double* ux_bdry = createVector(G->n_ctrl);
    double* uy_bdry = createVector(G->n_ctrl);

    extractboundaryData0(G, g, u, ux_bdry, uy_bdry);

    double x, y, dx, dy, nx, ny;
    for (int i = 0; i < G->n_ctrl; i++) {
        v[i] = ux_bdry[i] * G->ctrlPts_nxny[i][0] + uy_bdry[i] * G->ctrlPts_nxny[i][1];
    }

    freeVector(uy_bdry);
    freeVector(ux_bdry);
}














void findInteriorPoints(CartesianGridAndControlPoints* G, double (*g)(double, double), int bdry_idx, double** u, double& x, double& y, double** coord, double* grid_u, int& count) 
{
    double xi = G->ctrl_points[bdry_idx][0];
    double eta = G->ctrl_points[bdry_idx][1];
    
    int i = static_cast<int>((xi - G->Grid_x[0]) / G->dx + 0.5); 
    int j = static_cast<int>((eta - G->Grid_y[0]) / G->dy + 0.5); 

    if (i >= 0 && i < G->I + 1 && j >= 0 && j < G->J + 1 && G->interior[i][j]) {
        i = i; j = j;
    } else if (i - 1 >= 0 && i - 1 < G->I + 1 && j >= 0 && j < G->J + 1 && G->interior[i - 1][j]) {
        i = i - 1; j = j;
    } else if (i + 1 >= 0 && i + 1 < G->I + 1 && j >= 0 && j < G->J + 1 && G->interior[i + 1][j]) {
        i = i + 1; j = j;
    } else if (i >= 0 && i < G->I + 1 && j - 1 >= 0 && j - 1 < G->J + 1 && G->interior[i][j - 1]) {
        i = i; j = j - 1;
    } else if (i >= 0 && i < G->I + 1 && j + 1 >= 0 && j + 1 < G->J + 1 && G->interior[i][j + 1]) {
        i = i; j = j + 1;
    } else {
        std::cout << "Warning about finding interior node!" << std::endl;
    }

    coord[0][0] = 0.0;  coord[0][1] = 0.0;  grid_u[0] = u[i][j];      // i, j
    // coord[0][0] = G->Grid_x[i];  coord[0][1] = G->Grid_y[j];  grid_u[0] = u[i][j];

    double xi0 = G->Grid_x[i]; 
    double eta0 = G->Grid_y[j];

    x = (xi - xi0) / G->dx;
    y = (eta - eta0) / G->dy;
    // x = xi;
    // y = eta;

    coord[1][0] = x;  coord[1][1] = y;  grid_u[1] = g(xi, eta);
    
    // std::cout << x << " " << y << " " << (xi - xi0) << " " << (eta - eta0) << " " << xi << " " << eta << std::endl;

    int directions[24][2] = {{-1, 0}, {1, 0}, {0, -1}, {0, 1}, {-1, -1}, {1, 1}, {-1, 1}, {1, -1}, 
                            {-2, 0}, {2, 0}, {0, -2}, {0, 2}, {-2, 1}, {-1, 2}, {-2, -1}, {-1, -2}, 
                            {1, -2}, {2, -1}, {2, 1}, {1, 2}, {-2, -2}, {-2, 2}, {2, -2}, {2, 2}}; 
    // std::cout << "flag1" << std::endl;
    count = 2;  
    for (int k = 0; k < 24; ++k) {
        int di = directions[k][0]; 
        int dj = directions[k][1];

        int ni = i + di;
        int nj = j + dj;
        // std::cout << ni << " " << nj << " " << std::endl;
        if (ni >= 0 && ni < G->I + 1 && nj >= 0 && nj < G->J + 1 && G->interior[ni][nj]) {
            coord[count][0] = di / 1.; 
            coord[count][1] = dj / 1.; 
            // coord[count][0] = G->Grid_x[ni]; 
            // coord[count][1] = G->Grid_y[nj]; 
            grid_u[count] = u[ni][nj];
            count++;
            // if (count >= 6) {
            //     break;
            // }
        }
        // std::cout << ni << " " << nj << " " << std::endl;
    }
    // std::cout << "flag2" << std::endl;
    // if (count < 6) {
    //     std::cout << "Warning: Less than 6 interior points found!" << std::endl;
    // }
    
}

void computeBoundaryValues(double x, double y, double** coord, double* grid_u, int count, double& u, double& ux, double& uy) {
    // Arrays to store distances and indices of the points
    double* distances = createVector(count);
    int* indices = createIntVector(count);

    // Calculate distances from (x, y) to each point in coor
    for (int i = 0; i < count; i++) {
        distances[i] = (x - coord[i][0]) * (x - coord[i][0]) + (y - coord[i][1]) * (y - coord[i][1]);
        indices[i] = i;  // Store original index of each point
    }

    // Sort points by distance to (x, y)
    std::sort(indices, indices + count, [&](int i1, int i2) {
        return distances[i1] < distances[i2];
    });

    if (count > 10) {
        // std::cout << count << std::endl;
        count = 10;
    }

    double** A = createMatrix(count, 6);
    double** AT = createMatrix(6, count);
    double** ATA = createMatrix(6, 6);
    double* b = createVector(count);
    double* ATb = createVector(6);
    double* sol = createVector(6);
    for (int i = 0; i < count; i++) {
        int idx = indices[i];  
        A[i][0] = coord[i][0] * coord[i][0];
        A[i][1] = coord[i][0] * coord[i][1];
        A[i][2] = coord[i][1] * coord[i][1];
        A[i][3] = coord[i][0];
        A[i][4] = coord[i][1];
        A[i][5] = 1;
        b[i] = grid_u[i];
    }
 
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < count; j++) {
            AT[i][j] = A[j][i];
        }
    }

    MatrixMulMatrix_L(AT, A, ATA, 6, count, 6);
    MatrixMulVector_L(AT, b, ATb, 6, count);
    SolveLinearSystemByCG(ATA, ATb, sol, 6, 0);
    u = sol[0] * x * x + sol[1] * x * y + sol[2] * y * y + sol[3] * x + sol[4] * y + sol[5];
    ux = 2 * sol[0] * x + sol[1] * y + sol[3];
    uy = sol[1] * x + 2 * sol[2] * y + sol[4];
    
    freeMatrix(A, count);
    freeMatrix(AT, 6);
    freeMatrix(ATA, 6);
    freeVector(b);
    freeVector(ATb);
    freeVector(sol);

}

void extractboundaryData(CartesianGridAndControlPoints* G, double (*g)(double, double), double** w, double* ux, double* uy) 
{
    double x, y;
    double** coord = createMatrix(24, 2);
    double* grid_u = createVector(24);
    int count;
    double temp;

    for (int i = 0; i < G->n_ctrl; i++) {
        count = 0;
        findInteriorPoints(G, g, i, w, x, y, coord, grid_u, count);

        computeBoundaryValues(x, y, coord, grid_u, count, temp, ux[i], uy[i]);
        ux[i] /= G->dx;
        uy[i] /= G->dy;
    }

    freeVector(grid_u);
    freeMatrix(coord, 24);
}

void getNormalDerivatives(CartesianGridAndControlPoints* G, double (*g)(double, double), double** u, double* v)
{
    double* ux_bdry = createVector(G->n_ctrl);
    double* uy_bdry = createVector(G->n_ctrl);

    extractboundaryData(G, g, u, ux_bdry, uy_bdry);

    double x, y, dx, dy, nx, ny;
    for (int i = 0; i < G->n_ctrl; i++) {
        v[i] = ux_bdry[i] * G->ctrlPts_nxny[i][0] + uy_bdry[i] * G->ctrlPts_nxny[i][1];
    }

    freeVector(uy_bdry);
    freeVector(ux_bdry);
}