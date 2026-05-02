#ifndef EXTRACT_BOUNDARY_POINTS_H
#define EXTRACT_BOUNDARY_POINTS_H

struct Point {
    double x, y;
};

// Extracts intersection points where u=0 on a rectangular grid.
// Returns the number of points found, and stores them in outputPoints (caller must free).
int extractZeroLevelSet(double** u, int Nx, int Ny, double* x, double* y, Point* pts);
int extractZeroLevelSet(int** u, int Nx, int Ny, double* x, double* y, Point* pts);

// Sorts the points counterclockwise based on angle with respect to their centroid.
void sortPointsCCW(Point* pts, int N);

#endif // EXTRACT_BOUNDARY_POINTS_H
