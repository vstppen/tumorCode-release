#include "ExtractBoundaryPoints.h"
#include <cmath>
#include <cstdlib>
#include <vector>

int extractZeroLevelSet(double** u, int Nx, int Ny, double* x, double* y, Point* pts) {
    int count = 0;

    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            if (std::abs(u[i][j]) < 1.e-8) {
                bool hasPositiveNeighbor = false;

                if (i > 0     && u[i - 1][j] > 1.e-8) hasPositiveNeighbor = true;
                if (i < Nx-1  && u[i + 1][j] > 1.e-8) hasPositiveNeighbor = true;
                if (j > 0     && u[i][j - 1] > 1.e-8) hasPositiveNeighbor = true;
                if (j < Ny-1  && u[i][j + 1] > 1.e-8) hasPositiveNeighbor = true;

                if (hasPositiveNeighbor) {
                    pts[count].x = x[i];
                    pts[count].y = y[j];
                    count++;
                }
            }
        }
    }

    return count;
}

void sortPointsCCW(Point* pts, int N) {
    // Compute centroid of the points
    double cx = 0, cy = 0;
    for (int i = 0; i < N; ++i) {
        cx += pts[i].x;
        cy += pts[i].y;
    }
    cx /= N;
    cy /= N;

    // Bubble sort based on angle from centroid
    for (int i = 0; i < N - 1; ++i) {
        for (int j = 0; j < N - i - 1; ++j) {
            double angle1 = atan2(pts[j].y - cy, pts[j].x - cx);
            double angle2 = atan2(pts[j+1].y - cy, pts[j+1].x - cx);
            if (angle1 > angle2) {
                Point tmp = pts[j];
                pts[j] = pts[j+1];
                pts[j+1] = tmp;
            }
        }
    }
}