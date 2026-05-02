#ifndef BOUNDARYINCLUSIONCHECK_H
#define BOUNDARYINCLUSIONCHECK_H

#include <cmath>
#include <algorithm>

inline bool isPointInPolygon(double x, double y, double** polygon, int n) {
    int cnt = 0;
    for (int i = 0; i < n; ++i) {
        double x1 = polygon[i][0], y1 = polygon[i][1];
        double x2 = polygon[(i + 1) % n][0], y2 = polygon[(i + 1) % n][1];

        // Skip horizontal edges (avoid ambiguous intersections)
        if (fabs(y1 - y2) < 1e-12) continue;

        // Ensure y1 <= y2
        if (y1 > y2) {
            std::swap(x1, x2);
            std::swap(y1, y2);
        }

        // Check if y is within the edge's vertical range
        if (y < y1 || y >= y2) continue;

        // Compute intersection x-coordinate
        double x_inters = (y - y1) * (x2 - x1) / (y2 - y1) + x1;

        if (x_inters >= x)
            cnt++;
    }
    return (cnt % 2 == 1);  // Inside if odd intersections
}

inline bool isBoundaryInside(double** ctrl_points_0, int n_ctrl_0,
                             double** ctrl_points_1, int n_ctrl_1) {
    for (int i = 0; i < n_ctrl_0; ++i) {
        double x = ctrl_points_0[i][0];
        double y = ctrl_points_0[i][1];
        if (!isPointInPolygon(x, y, ctrl_points_1, n_ctrl_1)) {
            return false;  
        }
    }
    return true;
}

#endif 
