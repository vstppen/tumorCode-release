#include <iostream>
#include <math.h>
#include <string.h>
#include <cstring>

#include "LinearAlgebra.h"
#include "CartesianGridAndControlPoints.h"
#include "TangentialDerivatives.h"
#include "BirkoffInterpolation.h"
#include "LU.h"
#include "QR.h"
#include "CubicSpline.h"

CartesianGridAndControlPoints::CartesianGridAndControlPoints(double x_min_, double x_max_, double y_min_, double y_max_, int I_, int J_, 
                                double** ctrl_points_, int n_ctrl_, int M_, 
                                const std::string& filename_ctrlPts_,
                                const std::string& filename_ctrl_nxny_, 
                                const std::string& filename_bdryPts_) 
{   
    x_min = x_min_;
    x_max = x_max_;
    y_min = y_min_;
    y_max = y_max_;
    I = I_;
    J = J_;
    
    n_ctrl = n_ctrl_;
    ctrl_points = createMatrix(n_ctrl, 2);
    for (int i = 0; i < n_ctrl; i++) {
        ctrl_points[i][0] = ctrl_points_[i][0];
        ctrl_points[i][1] = ctrl_points_[i][1];
    }

    M = M_;

    dx = (x_max - x_min) / I;
    dy = (y_max - y_min) / J;

    safetyCheck();

    Grid_x = createVector(I + 1);
    Grid_x[0] = x_min; 
    for (int i = 0; i < I; i++) {
        Grid_x[i + 1] = Grid_x[i] + dx; 
    }

    Grid_y = createVector(J + 1);
    Grid_y[0] = y_min; 
    for (int i = 0; i < J; i++) {
        Grid_y[i + 1] = Grid_y[i] + dy; 
    }

    M_x = createVector(n_ctrl);
    alpha_x = createVector(n_ctrl);
    beta_x = createVector(n_ctrl);
    M_y = createVector(n_ctrl);
    alpha_y = createVector(n_ctrl);
    beta_y = createVector(n_ctrl);

    periodicCubicSplineCurveInterpolation(ctrl_points, M_x, alpha_x, beta_x, M_y, alpha_y, beta_y, n_ctrl);

    generate_boundary_points(filename_ctrlPts_, filename_ctrl_nxny_, filename_bdryPts_);
    
    identifyInteriorGridNodes();
    identifyIrregularGridNodes();
    identifyNearBdryGridNodes();

    findIntersectionPoints();

}

CartesianGridAndControlPoints::~CartesianGridAndControlPoints()
{
    freeMatrix(ctrl_points, n_ctrl);
    freeVector(Grid_x);
    freeVector(Grid_y);

    freeVector(M_x);
    freeVector(alpha_x);
    freeVector(beta_x);
    freeVector(M_y);
    freeVector(alpha_y);
    freeVector(beta_y);

    freeVector(bdry_theta);
    freeMatrix(xy, M);
    freeMatrix(dxdy, M);
    freeMatrix(nxny, M);
    freeMatrix(ctrlPts_nxny, n_ctrl);
    freeMatrix(ddxddy, M);

    freeBoolMatrix(interior, I + 1);
    freeIntMatrix(int_node_index, int_node_num);

    freeBoolMatrix(irregular, I + 1);
    freeIntMatrix(irreg_node_index, irreg_node_num);

    freeBoolMatrix(near_bdry, I + 1);
    freeIntMatrix(near_node_index, near_node_num);
    freeIntVector(near_node_proj_to_bdry_index);

    freeIntMatrix(intersect_index, intersect_number);
    freeVector(intersect_theta);
    freeMatrix(intersect_coord, intersect_number);
    freeMatrix(intersect_tan, intersect_number);
    freeMatrix(intersect_nml, intersect_number);
}

inline bool isPowerOfTwo(int x) {
    return (x > 0) && ((x & (x - 1)) == 0);
}

void CartesianGridAndControlPoints::safetyCheck()
{   
    // 1. 
    if (I <= 8 or !isPowerOfTwo(I)) {
        std::cout <<"Failed to initialize the grid and the periodic cubic curve, since I is too small or it isn't a power of 2!" << std::endl;
        exit(EXIT_FAILURE);
    }

    if (J <= 8 or !isPowerOfTwo(J)) {
        std::cout <<"Failed to initialize the grid and the periodic cubic curve, since J is too small or it isn't a power of 2!" << std::endl;
        exit(EXIT_FAILURE);
    }

    if (M <= 8) {
        std::cout <<"Failed to initialize the grid and the periodic cubic curve, since M is too small!" << std::endl;
        exit(EXIT_FAILURE);
    }

    // 2.
    for (int i = 0; i < n_ctrl; i++) {
        if (ctrl_points[i][0] <= x_min + 2 * dx or ctrl_points[i][0] >= x_max - 2 * dx or ctrl_points[i][1] <= y_min + 2 * dy or ctrl_points[i][1] >= y_max - 2 * dy) {
            std::cout <<"Failed to initialize the grid and the periodic cubic curve, since some control points are out of the grid range!" << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    // 3. 
    double length = 0;
    for (int i = 0; i < n_ctrl - 1; i++) {
        length += sqrt((ctrl_points[i][0] - ctrl_points[i + 1][0]) * (ctrl_points[i][0] - ctrl_points[i + 1][0]) 
                    + (ctrl_points[i][1] - ctrl_points[i + 1][1]) * (ctrl_points[i][1] - ctrl_points[i + 1][1]));
    }
    if (length < 4e-1) {
        std::cout <<"Failed to initialize the grid and the periodic cubic curve, since the control points are too concentrated!" << std::endl;
        exit(EXIT_FAILURE);
    }

    // 4.
    double sum = 0.0;
    for (int i = 0; i < n_ctrl; ++i) {
        int next_i = (i + 1) % n_ctrl; // Wrap around for the last point
        double xi = ctrl_points[i][0];
        double yi = ctrl_points[i][1];
        double xj = ctrl_points[next_i][0];
        double yj = ctrl_points[next_i][1];
        sum += (xi * yj) - (xj * yi); // Shoelace formula component
    }

    if (sum <= 0.0) { // Includes degenerate cases (sum == 0)
        std::cout <<"Failed to initialize the grid and the periodic cubic curve, since the control points are not aligned counterclockwise!" << std::endl;
        exit(EXIT_FAILURE);
    }

}

void CartesianGridAndControlPoints::getPoint1(double t, double &x, double &y)
{
    x = periodicCubicSplineGetS(M_x, alpha_x, beta_x, n_ctrl, t);
    y = periodicCubicSplineGetS(M_y, alpha_y, beta_y, n_ctrl, t);
}

void CartesianGridAndControlPoints::getPoint2(double theta, double &x, double &y, double &tx, double &ty, double &nx, double &ny)
{
    x = periodicCubicSplineGetS(M_x, alpha_x, beta_x, n_ctrl, theta);
    y = periodicCubicSplineGetS(M_y, alpha_y, beta_y, n_ctrl, theta);

    tx = periodicCubicSplineGetDS(M_x, alpha_x, beta_x, n_ctrl, theta);
    ty = periodicCubicSplineGetDS(M_y, alpha_y, beta_y, n_ctrl, theta);

    double r_len = 1.0 / sqrt(tx * tx + ty * ty); 
    nx = ty * r_len;  ny = - tx * r_len; 
}

void CartesianGridAndControlPoints::getPoint3(double theta, double &x, double &y, double &dx, double &dy, double &ddx, double &ddy)
{
    x = periodicCubicSplineGetS(M_x, alpha_x, beta_x, n_ctrl, theta);
    y = periodicCubicSplineGetS(M_y, alpha_y, beta_y, n_ctrl, theta);

    dx = periodicCubicSplineGetDS(M_x, alpha_x, beta_x, n_ctrl, theta);
    dy = periodicCubicSplineGetDS(M_y, alpha_y, beta_y, n_ctrl, theta);

    ddx = periodicCubicSplineGetDDS(M_x, alpha_x, beta_x, n_ctrl, theta);
    ddy = periodicCubicSplineGetDDS(M_y, alpha_y, beta_y, n_ctrl, theta);
}

void CartesianGridAndControlPoints::generate_boundary_points(const std::string& filename_ctrlPts, const std::string& filename_ctrl_nxny, const std::string& filename_bdryPts) 
{
    bdry_delta = M_2PI / M;

    bdry_theta = createVector(M);  
    xy = createMatrix(M, 2);
    dxdy = createMatrix(M, 2);
    ctrlPts_nxny = createMatrix(n_ctrl, 2);
    nxny = createMatrix(M, 2);
    ddxddy = createMatrix(M, 2);

    double a, b, c, d;
    for (int i = 0; i < n_ctrl; i++) {
        getPoint2(i * M_2PI / n_ctrl, a, b, c, d, ctrlPts_nxny[i][0], ctrlPts_nxny[i][1]);
    }

    for (int i = 0; i < M; i++) {
        bdry_theta[i] = i * bdry_delta;
        getPoint2(bdry_theta[i], xy[i][0], xy[i][1], dxdy[i][0], dxdy[i][1], nxny[i][0], nxny[i][1]);
        getPoint3(bdry_theta[i], xy[i][0], xy[i][1], dxdy[i][0], dxdy[i][1], ddxddy[i][0], ddxddy[i][1]);
    }

    saveMatrixToFile(const_cast<const double**>(ctrl_points), n_ctrl, 2, filename_ctrlPts);
    saveMatrixToFile(const_cast<const double**>(ctrlPts_nxny), n_ctrl, 2, filename_ctrl_nxny);
    saveMatrixToFile(const_cast<const double**>(xy), M, 2, filename_bdryPts);
}

int CartesianGridAndControlPoints::findClosestIndex(double p, double q)
{
    double dist = 10000.;
    double dist_temp;
    int i_temp;
    for (int i = 0; i < M; i++) {
        dist_temp = (p - xy[i][0]) * (p - xy[i][0]) + (q - xy[i][1]) * (q - xy[i][1]);
        if (dist_temp < dist) {
            i_temp = i;
            dist = dist_temp;
        }
    }

    return i_temp;
}

void CartesianGridAndControlPoints::findClosestPointRecursive(double p, double q, double t_start, double t_end, 
                                                              double &t0, double &x0, double &y0) 
{
    double tol = EPSILON12;
    if (fabs(t_end - t_start) < tol) { 
        t0 = (t_start + t_end) / 2.0;
        getPoint1(t0, x0, y0);
        return;
    }

    double t_mid1 = t_start + (t_end - t_start) / 3.0;
    double t_mid2 = t_start + 2.0 * (t_end - t_start) / 3.0;
    
    double x1, y1, x2, y2;
    getPoint1(t_mid1, x1, y1);
    getPoint1(t_mid2, x2, y2);
    
    double dist1 = (p - x1) * (p - x1) + (q - y1) * (q - y1);
    double dist2 = (p - x2) * (p - x2) + (q - y2) * (q - y2);
    
    if (dist1 < dist2)
        findClosestPointRecursive(p, q, t_start, t_mid2, t0, x0, y0);
    else
        findClosestPointRecursive(p, q, t_mid1, t_end, t0, x0, y0);
}

void CartesianGridAndControlPoints::findClosestPoint(double p, double q, double &t0, double &x0, double &y0)
{
    int i0 = findClosestIndex(p, q);
    double t_start = bdry_theta[i0] - 2. / 3. * bdry_delta;
    double t_end = bdry_theta[i0] + 2. / 3. * bdry_delta;
    findClosestPointRecursive(p, q, t_start, t_end, t0, x0, y0);
}

double CartesianGridAndControlPoints::computeParameter(double x, double y) 
{
    double result, x0, y0;
    findClosestPoint(x, y, result, x0, y0);
    return result;
}

bool CartesianGridAndControlPoints::interiorOrNot(double p, double q, double x, double y, double nx, double ny) {
    double s = (p - x) * nx + (q - y) * ny;

    if (s > 0.0) {
        return false;
    } else {
        return true;
    }
}

bool CartesianGridAndControlPoints::interiorOrNot(double p, double q) {
    int i0 = findClosestIndex(p, q);
    double t = bdry_theta[i0];
    double x, y, tx, ty, nx, ny, s;

    getPoint2(t, x, y, tx, ty, nx, ny);
    s = (p - x) * nx + (q - y) * ny;

    if (s > 0.0) {
        return false;
    } else {
        return true;
    }
}

void CartesianGridAndControlPoints::identifyInteriorGridNodes()
{   
    interior = createBoolMatrix(I + 1, J + 1);
    int_node_num = 0;
    for (int i = 0; i < I + 1; i++) {
        for (int j = 0; j < J + 1; j++) {
            interior[i][j] = interiorOrNot(Grid_x[i], Grid_y[j]);
            if (interior[i][j] == true) {
                int_node_num++;
            }
        }
    }

    int count = 0;
    int_node_index = createIntMatrix(int_node_num, 2);
    for (int i = 1; i < I; i++) {
        for (int j = 1; j < J; j++) {
            if (interior[i][j]) {
                int_node_index[count][0] = i; 
                int_node_index[count][1] = j; 
                count++;
            }
        }
    }

    if (I == 64) {
        double** tt = createMatrix(int_node_num, 2);
        for (int i = 0; i < int_node_num; i++) {
            tt[i][0] = Grid_x[int_node_index[i][0]];
            tt[i][1] = Grid_y[int_node_index[i][1]];
        }
        saveMatrixToFile(const_cast<const double**>(tt), int_node_num, 2, "./result/points/int_points.txt");
        freeMatrix(tt, int_node_num);
    }
}

/**
 *****************************************************************************
 * irregular: the bool pointer we desired                                    *
 * irreg_node_num: the number of irregular points which we desired           *
 * irreg_node_index: the indices of irregular points which we desired        *
 *****************************************************************************
 **/
void CartesianGridAndControlPoints::identifyIrregularGridNodes()
{
    irregular = createBoolMatrix(I + 1, J + 1);
    int count = 0; 
    for (int i = 1; i < I; i++) {
        for (int j = 1; j < J; j++) {
            if ((interior[i][j] != interior[i + 1][j]) ||
                (interior[i][j] != interior[i - 1][j]) ||
                (interior[i][j] != interior[i][j + 1]) ||
                (interior[i][j] != interior[i][j - 1])) {
                irregular[i][j] = true; 
                count++;
            }
        }
      }

    irreg_node_num = count; 
    irreg_node_index = createIntMatrix(irreg_node_num, 2);

    count = 0;
    for (int i = 1; i < I; i++) {
        for (int j = 1; j < J; j++) {
            if (irregular[i][j]) {
                irreg_node_index[count][0] = i; 
                irreg_node_index[count][1] = j; 
                count++;
            }
        }
    }
}

void CartesianGridAndControlPoints::identifyNearBdryGridNodes()
{
    near_bdry = createBoolMatrix(I + 1, J + 1);
    int count = 0; 
    for (int i = 3; i < I - 2; i++) {
        for (int j = 3; j < J - 2; j++) {
            if ((interior[i][j] != interior[i + 1][j]) ||
                (interior[i][j] != interior[i + 2][j]) ||
                (interior[i][j] != interior[i - 1][j]) ||
                (interior[i][j] != interior[i - 2][j]) ||
                (interior[i][j] != interior[i][j + 1]) ||
                (interior[i][j] != interior[i][j + 2]) ||
                (interior[i][j] != interior[i][j - 1]) ||
                (interior[i][j] != interior[i][j - 2]) ) {
                near_bdry[i][j] = true; 
                count++;
            }
        }
      }

    near_node_num = count; 
    near_node_index = createIntMatrix(near_node_num, 2);

    count = 0;
    for (int i = 3; i < I - 2; i++) {
        for (int j = 3; j < J - 2; j++) {
            if (near_bdry[i][j]) {
                near_node_index[count][0] = i; 
                near_node_index[count][1] = j; 
                count++;
            }
        }
    }

    near_node_proj_to_bdry_index = createIntVector(near_node_num);
    double temp_theta = 0.;
    for (int i = 0; i < near_node_num; i++) {
        temp_theta = computeParameter(Grid_x[near_node_index[i][0]], Grid_y[near_node_index[i][1]]);
        near_node_proj_to_bdry_index[i] = round(temp_theta / bdry_delta);
        if (near_node_proj_to_bdry_index[i] > M - 1) {
            near_node_proj_to_bdry_index[i] -= M;
        }
        // std::cout << near_node_proj_to_bdry_index[i] << " ";
    }

}

// find the x-coordinate of intersection points between (a, q) and (b, q)
double CartesianGridAndControlPoints::findXIntersectByBisection(double a, double b, double q, double tol)
{
    int ia = findClosestIndex(a, q);
    int ib = findClosestIndex(b, q);
    double c = 0.0;

    if (ia == ib) {

        double t = bdry_theta[ia];
        double x, y, tx, ty, nx, ny, s;
        getPoint2(t, x, y, tx, ty, nx, ny);
        
        bool s_l = interiorOrNot(a, q, x, y, nx, ny);
        bool s_r = interiorOrNot(b, q, x, y, nx, ny);

        if (s_l == s_r) {
            std::cout << "The resolution of the rectangular grid is not enough! X " 
            << a << " " << b << "  " << q << "  " << s_l << "  " <<  s_r << "  " 
            << findClosestIndex(a, q) << "  " << findClosestIndex(b, q) << std::endl;
            exit(EXIT_FAILURE);
        }

        bool done = false; 
        do {
            c = 0.5 * (a + b); 
            bool s_c = interiorOrNot(c, q, x, y, nx, ny);
            if (fabs(c - a) < tol) {
                done = true;
            } else {
                if (s_c == s_l) {
                    a = c;  s_l = s_c; 
                } else {
                    b = c;  s_r = s_c; 
                }
                if (fabs(b - a) < tol) {
                    done = true;
                }
            }
        } while (!done); 

    } else {
    
        bool s_l = interiorOrNot(a, q);
        bool s_r = interiorOrNot(b, q);

        if (s_l == s_r) {
            std::cout << "The resolution of the rectangular grid is not enough! X " 
            << a << " " << b << "  " << q << "  " << s_l << "  " <<  s_r << "  " 
            << findClosestIndex(a, q) << "  " << findClosestIndex(b, q) << std::endl;
            exit(EXIT_FAILURE);
        }

        bool done = false; 
        do {
            c = 0.5 * (a + b); 
            bool s_c = interiorOrNot(c, q);
            if (fabs(c - a) < tol) {
                done = true;
            } else {
                if (s_c == s_l) {
                    a = c;  s_l = s_c; 
                } else {
                    b = c;  s_r = s_c; 
                }
                if (fabs(b - a) < tol) {
                    done = true;
                }
            }
        } while (!done); 

    }

    return c; 
}

// find the y-coordinate of intersection points between (p, a) and (p, b)
double CartesianGridAndControlPoints::findYIntersectByBisection(double p, double a, double b, double tol)
{
    int ia = findClosestIndex(p, a);
    int ib = findClosestIndex(p, b);
    double c = 0.0;

    if (ia == ib) {

        double t = bdry_theta[ia];
        double x, y, tx, ty, nx, ny, s;
        getPoint2(t, x, y, tx, ty, nx, ny);
        
        bool s_l = interiorOrNot(p, a, x, y, nx, ny);
        bool s_u = interiorOrNot(p, b, x, y, nx, ny);

        if (s_l == s_u) {
            std::cout << "The resolution of the rectangular grid is not enough! Y " 
            << p << " " << a << "  " << b << "  " << s_l << "  " <<  s_u << "  " 
            << findClosestIndex(p, a) << "  " << findClosestIndex(p, b) << std::endl;
            exit(EXIT_FAILURE);
        }

        bool done = false; 
        do {
            c = 0.5 * (a + b); 
            bool s_c = interiorOrNot(p, c, x, y, nx, ny);
            if (fabs(c - a) < tol) {
                done = true;
            } else {
                if (s_c == s_l) {
                    a = c;  s_l = s_c; 
                } else {
                    b = c;  s_u = s_c; 
                }
                if (fabs(b - a) < tol) {
                    done = true;
                }
            }
        } while (!done); 

    } else {
    
        bool s_l = interiorOrNot(p, a);
        bool s_u = interiorOrNot(p, b);

        if (s_l == s_u) {
            std::cout << "The resolution of the rectangular grid is not enough! Y " 
            << p << " " << a << "  " << b << "  " << s_l << "  " <<  s_u << "  " 
            << findClosestIndex(p, a) << "  " << findClosestIndex(p, b) << std::endl;
            exit(EXIT_FAILURE);
        }

        bool done = false; 
        do {
            c = 0.5 * (a + b); 
            bool s_c = interiorOrNot(p, c);
            if (fabs(c - a) < tol) {
                done = true;
            } else {
                if (s_c == s_l) {
                    a = c;  s_l = s_c; 
                } else {
                    b = c;  s_u = s_c; 
                }
                if (fabs(b - a) < tol) {
                    done = true;
                }
            }
        } while (!done); 

    }

    return c; 
}

/**
 *****************************************************************************************************
 * intersect_num: the number of intersection points on vertical and horizontal direction             *
 * intersect_idx: the x-idx and y-idx of intersection points                                         *
 * intersect_theta: the parameter of intersection points                                             *
 * intersect_coord: the coordinates (x, y) of intersection points                                    * 
 * intersect_tan: (dx, dy) of intersection points                                                    *
 * intersect_nml: (nx, ny) of intersection points                                                    *
 *****************************************************************************************************
 **/
void CartesianGridAndControlPoints::findIntersectionPoints()
{
    double tol = 1.0E-12; 

    int count[2] = {0, 0};

    int m2 = 2 * irreg_node_num; 

    int (*idx_h)[2] = new int[m2][2];
    int (*idx_v)[2] = new int[m2][2]; 
    for (int ell = 0; ell < m2; ell++) {
        idx_h[ell][0] = idx_h[ell][1] = - 1;
        idx_v[ell][0] = idx_v[ell][1] = - 1; 
    }

    double *theta_h = new double[m2]; 
    double *theta_v = new double[m2]; 
    for (int ell = 0; ell < m2; ell++) {
        theta_h[ell] = 0.0; 
        theta_v[ell] = 0.0; 
    }

    double (*tan_h)[2] = new double[m2][2]; 
    double (*tan_v)[2] = new double[m2][2];
    for (int ell = 0; ell < m2; ell++) {
        tan_h[ell][0] = 0.0;  tan_h[ell][1] = 0.0;
        tan_v[ell][0] = 0.0;  tan_v[ell][1] = 0.0; 
    }

    double (*nml_h)[2] = new double[m2][2];
    double (*nml_v)[2] = new double[m2][2]; 
    for (int ell = 0; ell < m2; ell++) {
        nml_h[ell][0] = 0.0;  nml_h[ell][1] = 0.0; 
        nml_v[ell][0] = 0.0;  nml_v[ell][1] = 0.0; 
    }

    double (*coord_h)[2] = new double[m2][2]; 
    double (*coord_v)[2] = new double[m2][2]; 
    for (int ell = 0; ell < m2; ell++) {
        coord_h[ell][0] = Grid_x[0];  coord_h[ell][1] = Grid_y[0]; 
        coord_v[ell][0] = Grid_x[0];  coord_v[ell][1] = Grid_y[0];
    }

    double p1, q1, tx1, ty1, nx1, ny1;
    double p2, q2, tx2, ty2, nx2, ny2; 

    for (int ell = 0; ell < irreg_node_num; ell++) {

        int i = irreg_node_index[ell][0];
        int j = irreg_node_index[ell][1];

        // check the edge on the right
        if (interior[i][j] != interior[i + 1][j]) {

        p1 = findXIntersectByBisection(Grid_x[i], Grid_x[i + 1], Grid_y[j], tol); 

        double s = computeParameter(p1, Grid_y[j]);

        getPoint2(s, p1, q1, tx1, ty1, nx1, ny1);

        idx_h[count[0]][0] = i; 
        idx_h[count[0]][1] = j;
        coord_h[count[0]][0] = p1;
        coord_h[count[0]][1] = q1;
        tan_h[count[0]][0] = tx1; 
        tan_h[count[0]][1] = ty1;
        nml_h[count[0]][0] = nx1; 
        nml_h[count[0]][1] = ny1;
        theta_h[count[0]] = s; 
        count[0]++; 
        }

        // check the edge above
        if (interior[i][j] != interior[i][j + 1]) {

        q2 = findYIntersectByBisection(Grid_x[i], Grid_y[j], Grid_y[j + 1], tol);

        double s = computeParameter(Grid_x[i], q2);

        getPoint2(s, p2, q2, tx2, ty2, nx2, ny2); 

        idx_v[count[1]][0] = i; 
        idx_v[count[1]][1] = j; 
        coord_v[count[1]][0] = p2; 
        coord_v[count[1]][1] = q2;
        tan_v[count[1]][0] = tx2;
        tan_v[count[1]][1] = ty2;
        nml_v[count[1]][0] = nx2;
        nml_v[count[1]][1] = ny2;
        theta_v[count[1]] = s; 
        count[1]++; 
        }
    }

    /////////////////////////////////////////////////////////////////////////////
    intersect_num[0] = count[0];
    intersect_num[1] = count[1];
    intersect_number = count[0] + count[1];

    intersect_index = createIntMatrix(intersect_number, 2);
    for (int i = 0; i < count[0]; i++) {
        intersect_index[i][0] = idx_h[i][0];
        intersect_index[i][1] = idx_h[i][1]; 
    }
    for (int j = 0, i = count[0]; j < count[1]; i++, j++) {
        intersect_index[i][0] = idx_v[j][0]; 
        intersect_index[i][1] = idx_v[j][1];
    }

    intersect_coord = createMatrix(intersect_number, 2);
    for (int i = 0; i < count[0]; i++) {
        intersect_coord[i][0] = coord_h[i][0]; 
        intersect_coord[i][1] = coord_h[i][1];
    }
    for (int j = 0, i = count[0]; j < count[1]; i++, j++) {
        intersect_coord[i][0] = coord_v[j][0];
        intersect_coord[i][1] = coord_v[j][1];
    }

    intersect_nml = createMatrix(intersect_number, 2);
    for (int i = 0; i < count[0]; i++) {
        intersect_nml[i][0] = nml_h[i][0];
        intersect_nml[i][1] = nml_h[i][1];
    }
    for (int j = 0, i = count[0]; j < count[1]; i++, j++) {
        intersect_nml[i][0] = nml_v[j][0];
        intersect_nml[i][1] = nml_v[j][1];
    }

    intersect_tan = createMatrix(intersect_number, 2); 
    for (int i = 0; i < count[0]; i++) {
        intersect_tan[i][0] = tan_h[i][0];
        intersect_tan[i][1] = tan_h[i][1];
    }
    for (int j = 0, i = count[0]; j < count[1]; i++, j++) {
        intersect_tan[i][0] = tan_v[j][0]; 
        intersect_tan[i][1] = tan_v[j][1];
    }

    intersect_theta = createVector(intersect_number);
    for (int i = 0; i < count[0]; i++) {
        intersect_theta[i] = theta_h[i];
    }
    for (int j = 0, i = count[0]; j < count[1]; i++, j++) {
        intersect_theta[i] = theta_v[j];
    }


    /////////////////////////////////////////////////////////////////////////////
    // free pointers

    delete[] theta_h; 
    delete[] theta_v; 

    delete[] coord_h; 
    delete[] coord_v;

    delete[] nml_h;
    delete[] nml_v; 

    delete[] tan_h;
    delete[] tan_v;

    delete[] idx_h;
    delete[] idx_v; 
}

bool CartesianGridAndControlPoints::isDuplicate(double** coor, int count, double x, double y, double epsilon = 1e-10) 
{
    for (int i = 0; i < count; ++i) {
        if (fabs(coor[i][0] - x) < epsilon && fabs(coor[i][1] - y) < epsilon) {
            return true;  
        }
    }
    return false;
}

void CartesianGridAndControlPoints::findGoodPoints(int i, int j, double** coor, double* u, double** solution, double (*g)(double, double), bool selfOrNot, int& count)
{   
    if (i >= 0 && i <= I && j >= 0 && j <= J) {
        if (selfOrNot & interior[i][j] & !near_bdry[i][j]) {
            if (!isDuplicate(coor, count, Grid_x[i], Grid_y[j])) {  // 如果坐标不重复
                coor[count][0] = Grid_x[i];
                coor[count][1] = Grid_y[j];
                u[count] = solution[i][j];
                count++;
            }
        }

        if (i + 1 <= I && interior[i + 1][j] & !near_bdry[i + 1][j]) {
            if (!isDuplicate(coor, count, Grid_x[i + 1], Grid_y[j])) {
                coor[count][0] = Grid_x[i + 1];
                coor[count][1] = Grid_y[j];
                u[count] = solution[i + 1][j];
                count++;
            }
        }

        if (i + 1 <= I && interior[i][j] & !interior[i + 1][j]) {
            double new_x = findXIntersectByBisection(Grid_x[i], Grid_x[i + 1], Grid_y[j], EPSILON12);
            if (!isDuplicate(coor, count, new_x, Grid_y[j])) {
                coor[count][0] = new_x;
                coor[count][1] = Grid_y[j];
                u[count] = g(coor[count][0], coor[count][1]);
                count++;
            }
        }

        if (i - 1 >= 0 && interior[i - 1][j] & !near_bdry[i - 1][j]) {
            if (!isDuplicate(coor, count, Grid_x[i - 1], Grid_y[j])) {
                coor[count][0] = Grid_x[i - 1];
                coor[count][1] = Grid_y[j];
                u[count] = solution[i - 1][j];
                count++;
            }
        }

        if (i - 1 >= 0 && interior[i][j] & !interior[i - 1][j]) {
            double new_x = findXIntersectByBisection(Grid_x[i - 1], Grid_x[i], Grid_y[j], EPSILON12);
            if (!isDuplicate(coor, count, new_x, Grid_y[j])) {
                coor[count][0] = new_x;
                coor[count][1] = Grid_y[j];
                u[count] = g(coor[count][0], coor[count][1]);
                count++;
            }
        }

        if (j + 1 <= J && interior[i][j + 1] & !near_bdry[i][j + 1]) {
            if (!isDuplicate(coor, count, Grid_x[i], Grid_y[j + 1])) {
                coor[count][0] = Grid_x[i];
                coor[count][1] = Grid_y[j + 1];
                u[count] = solution[i][j + 1];
                count++;
            }
        }

        if (j + 1 <= J && interior[i][j] & !interior[i][j + 1]) {
            double new_y = findYIntersectByBisection(Grid_x[i], Grid_y[j], Grid_y[j + 1], EPSILON12);
            if (!isDuplicate(coor, count, Grid_x[i], new_y)) {
                coor[count][0] = Grid_x[i];
                coor[count][1] = new_y;
                u[count] = g(coor[count][0], coor[count][1]);
                count++;
            }
        }

        if (j - 1 >= 0 && interior[i][j - 1] & !near_bdry[i][j - 1]) {
            if (!isDuplicate(coor, count, Grid_x[i], Grid_y[j - 1])) {
                coor[count][0] = Grid_x[i];
                coor[count][1] = Grid_y[j - 1];
                u[count] = solution[i][j - 1];
                count++;
            }
        }

        if (j - 1 >= 0 &&interior[i][j] & !interior[i][j - 1]) {
            double new_y = findYIntersectByBisection(Grid_x[i], Grid_y[j - 1], Grid_y[j], EPSILON12);
            if (!isDuplicate(coor, count, Grid_x[i], new_y)) {
                coor[count][0] = Grid_x[i];
                coor[count][1] = new_y;
                u[count] = g(coor[count][0], coor[count][1]);
                count++;
            }
        }
    }
}

void CartesianGridAndControlPoints::computeDDXandDN(double xi, double eta, double ddx[2], double dn[2])
{
    double theta = computeParameter(xi, eta); 

    double delta = 1.0E-5; 
    double delta2 = delta + delta;

    double s1 = theta - delta; 
    double s2 = theta + delta;

    double x1, y1, tx1, ty1, nx1, ny1; 
    double x2, y2, tx2, ty2, nx2, ny2;

    getPoint2(s1, x1, y1, tx1, ty1, nx1, ny1);
    getPoint2(s2, x2, y2, tx2, ty2, nx2, ny2);

    double x0, y0, dx0, dy0, ddx0, ddy0;

    getPoint3(theta, x0, y0, dx0, dy0, ddx0, ddy0); 

    ddx[0] = ddx0;  ddx[1] = ddy0;

    dn[0] = (nx2 - nx1) / delta2;
    dn[1] = (ny2 - ny1) / delta2; 
}

void CartesianGridAndControlPoints::computeJumps1(double x, double y, double tx, double ty, double nx, double ny, double psi, double d_phi, double du[2]) 
{
    double mat[2][2], rhs[2];

    mat[0][0] = tx;  mat[0][1] = ty;
    mat[1][0] = nx;  mat[1][1] = ny; 

    rhs[0] = d_phi;  rhs[1] = psi; 

    double det = mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0]; 

    du[0] = (mat[1][1] * rhs[0] - mat[0][1] * rhs[1]) / det; 
    du[1] = (mat[0][0] * rhs[1] - mat[1][0] * rhs[0]) / det; 
}

void CartesianGridAndControlPoints::computeJumps2(double x, double y, double s, double tx, double ty, double nx, 
                                            double ny, double d_psi, double dd_phi, const double du[2], double ddu[3])
{
  double ddX[2], dn[2]; 

  computeDDXandDN(x, y, ddX, dn);

  double ddx = ddX[0]; 
  double ddy = ddX[1]; 

  double dnx = dn[0]; 
  double dny = dn[1];

  double txtx = tx * tx;
  double txty = tx * ty;
  double tyty = ty * ty;

  double nxtx = nx * tx; 
  double nxty = nx * ty; 

  double nytx = ny * tx;
  double nyty = ny * ty;

  double mat[3][3]; 

  mat[0][0] = txtx; 
  mat[0][1] = tyty;
  mat[0][2] = txty + txty;

  mat[1][0] = nxtx;
  mat[1][1] = nyty; 
  mat[1][2] = nxty + nytx;

  mat[2][0] = 1.0; 
  mat[2][1] = 1.0; 
  mat[2][2] = 0.0; 

  double rhs[3]; 

  rhs[0] = dd_phi - ddx * du[0] - ddy * du[1];
  rhs[1] = d_psi - dnx * du[0] - dny * du[1]; 
  rhs[2] = 0.0;

  bool status = solveByQRdecomposition<3>(mat, rhs, ddu, 3); 
}

void CartesianGridAndControlPoints::computeJumps(const double *phi_vec, const double *psi_vec, int K,
                  double x0, double y0, double s, double tx, double ty,
                  double nx, double ny, double &j0, double j1[2], double j2[3]) 
{
  double psi, d_psi; 
  double phi, d_phi, dd_phi;

  computeFirstDerivative(bdry_delta, bdry_theta, psi_vec, K, s, psi, d_psi); 

  computeFirstAndSecondDerivatives(bdry_delta, bdry_theta, phi_vec, K, s, phi, d_phi, dd_phi); 

  double jmp1[2], jmp2[3];

  computeJumps1(x0, y0, tx, ty, nx, ny, psi, d_phi, jmp1); 

  computeJumps2(x0, y0, s, tx, ty, nx, ny, d_psi, dd_phi, jmp1, jmp2); 

  j0 = phi; 

  j1[0] = jmp1[0];
  j1[1] = jmp1[1]; 

  j2[0] = jmp2[0]; 
  j2[1] = jmp2[1]; 
  j2[2] = jmp2[2]; 
}

void CartesianGridAndControlPoints::computeJumps2(double x, double y, double s, double tx, double ty,
                                            double nx, double ny, double f, double ddu[3])
{
  double ddX[2], dn[2]; 

  computeDDXandDN(x, y, ddX, dn); 

  double ddx = ddX[0]; 
  double ddy = ddX[1]; 

  double dnx = dn[0]; 
  double dny = dn[1]; 

  double txtx = tx * tx; 
  double txty = tx * ty;
  double tyty = ty * ty; 

  double nxtx = nx * tx; 
  double nxty = nx * ty; 

  double nytx = ny * tx; 
  double nyty = ny * ty;

  double mat[3][3];

  mat[0][0] = txtx;
  mat[0][1] = tyty; 
  mat[0][2] = txty + txty;

  mat[1][0] = nxtx;
  mat[1][1] = nyty;
  mat[1][2] = nxty + nytx;

  mat[2][0] = 1.0;
  mat[2][1] = 1.0; 
  mat[2][2] = 0.0; 

  double rhs[3];

  rhs[0] = 0.0;
  rhs[1] = 0.0;
  rhs[2] = f; 

  bool status = solveByQRdecomposition<3>(mat, rhs, ddu, 3);
}

void CartesianGridAndControlPoints::computeJumps(double x0, double y0, double s, double tx, double ty, double nx, double ny, double f, double jmp2[3])
{
  computeJumps2(x0, y0, s, tx, ty, nx, ny, f, jmp2);
}

void CartesianGridAndControlPoints::makeCorrection(double (*F)(double, double), double **b)
{

  double d1, d2, c; 

  double jmp2[3];

  int M0 = intersect_num[0]; 
  int M1 = intersect_num[1]; 

  for (int ell = 0; ell < M0; ell++) {

    int i = intersect_index[ell][0];
    int j = intersect_index[ell][1];

    double p = intersect_coord[ell][0];
    double q = intersect_coord[ell][1];

    double s = intersect_theta[ell]; 
    double nx = intersect_nml[ell][0]; 
    double ny = intersect_nml[ell][1]; 
    double tx = intersect_tan[ell][0];
    double ty = intersect_tan[ell][1]; 

    double f = F(p, q); 
    
    // if (ell == 1) {
    //     std::cout << p << " " << Grid_y[j] << " " << s << " " << tx << " " << ty << " " << nx << " " << ny << " " << f << std::endl;
    // }

    // if (ell == 1) {
    //     std::cout << jmp2[0] << " " << jmp2[1] << " " << jmp2[2] << " " << std::endl;
    // }
    computeJumps(p, Grid_y[j], s, tx, ty, nx, ny, f, jmp2);
    // if (ell == 1) {
    //     std::cout << jmp2[0] << " " << jmp2[1] << " " << jmp2[2] << " " << std::endl;
    // }

    d1 = Grid_x[i + 1] - p; 
    d2 = 0.5 * d1 * d1;
    c = jmp2[0] * d2; 
    if (interior[i + 1][j]) {
      b[i][j] -= c; 
    } else {
      b[i][j] += c;
    }

    d1 = Grid_x[i] - p; 
    d2 = 0.5 * d1 * d1; 
    c = jmp2[0] * d2;
    if (interior[i][j]) {
      b[i + 1][j] -= c; 
    } else {
      b[i + 1][j] += c; 
    }
  }

  for (int ell = M0; ell < M0 + M1; ell++) {

    int i = intersect_index[ell][0]; 
    int j = intersect_index[ell][1];

    double p = intersect_coord[ell][0]; 
    double q = intersect_coord[ell][1]; 

    double s = intersect_theta[ell]; 
    double nx = intersect_nml[ell][0];
    double ny = intersect_nml[ell][1]; 
    double tx = intersect_tan[ell][0];
    double ty = intersect_tan[ell][1];

    double f = F(p, q); 

    // if (ell == M0 + 1) {
    //     std::cout << jmp2[0] << " " << jmp2[1] << " " << jmp2[2] << " " << std::endl;
    // }
    computeJumps(Grid_x[i], q, s, tx, ty, nx, ny, f, jmp2);
    // if (ell == M0 + 1) {
    //     std::cout << jmp2[0] << " " << jmp2[1] << " " << jmp2[2] << " " << std::endl;
    // }

    d1 = Grid_y[j + 1] - q; 
    d2 = 0.5 * d1 * d1; 
    c = jmp2[1] * d2;
    if (interior[i][j + 1]) {
      b[i][j] -= c;
    } else {
      b[i][j] += c;
    }

    d1 = Grid_y[j] - q; 
    d2 = 0.5 * d1 * d1;
    c = jmp2[1] * d2; 
    if (interior[i][j]) {
      b[i][j + 1] -= c;
    } else {
      b[i][j + 1] += c; 
    }
  }
}

void CartesianGridAndControlPoints::makeCorrection(double (*F)(double, double, double, double, int, int, double**, double, double), double** source, double **b)
{

  double d1, d2, c; 

  double jmp2[3];

  int M0 = intersect_num[0]; 
  int M1 = intersect_num[1]; 

  for (int ell = 0; ell < M0; ell++) {

    int i = intersect_index[ell][0];
    int j = intersect_index[ell][1];

    double p = intersect_coord[ell][0];
    double q = intersect_coord[ell][1];

    double s = intersect_theta[ell]; 
    double nx = intersect_nml[ell][0]; 
    double ny = intersect_nml[ell][1]; 
    double tx = intersect_tan[ell][0];
    double ty = intersect_tan[ell][1]; 

    double f = F(x_min, y_min, dx, dy, I, J, source, p, q); 
    
    // if (ell == 1) {
    //     std::cout << p << " " << Grid_y[j] << " " << s << " " << tx << " " << ty << " " << nx << " " << ny << " " << f << std::endl;
    // }

    // if (ell == 1) {
    //     std::cout << jmp2[0] << " " << jmp2[1] << " " << jmp2[2] << " " << std::endl;
    // }
    computeJumps(p, Grid_y[j], s, tx, ty, nx, ny, f, jmp2);
    // if (ell == 1) {
    //     std::cout << jmp2[0] << " " << jmp2[1] << " " << jmp2[2] << " " << std::endl;
    // }

    d1 = Grid_x[i + 1] - p; 
    d2 = 0.5 * d1 * d1;
    c = jmp2[0] * d2; 
    if (interior[i + 1][j]) {
      b[i][j] -= c; 
    } else {
      b[i][j] += c;
    }

    d1 = Grid_x[i] - p; 
    d2 = 0.5 * d1 * d1; 
    c = jmp2[0] * d2;
    if (interior[i][j]) {
      b[i + 1][j] -= c; 
    } else {
      b[i + 1][j] += c; 
    }
  }

  for (int ell = M0; ell < M0 + M1; ell++) {

    int i = intersect_index[ell][0]; 
    int j = intersect_index[ell][1];

    double p = intersect_coord[ell][0]; 
    double q = intersect_coord[ell][1]; 

    double s = intersect_theta[ell]; 
    double nx = intersect_nml[ell][0];
    double ny = intersect_nml[ell][1]; 
    double tx = intersect_tan[ell][0];
    double ty = intersect_tan[ell][1];

    double f = F(x_min, y_min, dx, dy, I, J, source, p, q); 

    // if (ell == M0 + 1) {
    //     std::cout << jmp2[0] << " " << jmp2[1] << " " << jmp2[2] << " " << std::endl;
    // }
    computeJumps(Grid_x[i], q, s, tx, ty, nx, ny, f, jmp2);
    // if (ell == M0 + 1) {
    //     std::cout << jmp2[0] << " " << jmp2[1] << " " << jmp2[2] << " " << std::endl;
    // }

    d1 = Grid_y[j + 1] - q; 
    d2 = 0.5 * d1 * d1; 
    c = jmp2[1] * d2;
    if (interior[i][j + 1]) {
      b[i][j] -= c;
    } else {
      b[i][j] += c;
    }

    d1 = Grid_y[j] - q; 
    d2 = 0.5 * d1 * d1;
    c = jmp2[1] * d2; 
    if (interior[i][j]) {
      b[i][j + 1] -= c;
    } else {
      b[i][j + 1] += c; 
    }
  }
}

bool CartesianGridAndControlPoints::computeBoundaryValues(double **u, double xi, double eta, const double jump[6], double bdry_u[12])
{
  double x0 = Grid_x[0];
  double y0 = Grid_y[0]; 

  int i = static_cast<int>((xi - x0) / dx + 0.5); 
  int j = static_cast<int>((eta - y0) / dy + 0.5); 

  double xi0 = x0 + i * dx; 
  double eta0 = y0 + j * dy;

  int dir[2];

  if (xi > xi0) {
    dir[0] = 1; 
  } else {
    dir[0] = - 1; 
  }
  if (eta > eta0) {
    dir[1] = 1;
  } else {
    dir[1] = - 1; 
  }

  double coord[6][2];

  coord[0][0] = 0.0;      coord[0][1] = 0.0; 

  coord[1][0] = 1.0;      coord[1][1] = 0.0; 
  coord[2][0] = - 1.0;    coord[2][1] = 0.0;

  coord[3][0] = 0.0;      coord[3][1] = 1.0;
  coord[4][0] = 0.0;      coord[4][1] = - 1.0;

  coord[5][0] = dir[0];   coord[5][1] = dir[1]; 

  bool inside[6];

  inside[0] = interior[i][j]; 

  inside[1] = interior[i + 1][j];
  inside[2] = interior[i - 1][j];

  inside[3] = interior[i][j + 1]; 
  inside[4] = interior[i][j - 1]; 

  inside[5] = interior[i + dir[0]][j + dir[1]];

  double grid_u[6]; 

  grid_u[0] = u[i][j];
  grid_u[1] = u[i + 1][j];
  grid_u[2] = u[i - 1][j];
  grid_u[3] = u[i][j + 1];
  grid_u[4] = u[i][j - 1];
  grid_u[5] = u[i + dir[0]][j + dir[1]]; 

  double jp[6]; 

  jp[0] = jump[0];

  jp[1] = jump[1] * dx;
  jp[2] = jump[2] * dy;

  double h0h0 = dx * dx; 
  double h1h1 = dy * dy;
  double h0h1 = dx * dy;

  jp[3] = jump[3] * h0h0;
  jp[4] = jump[4] * h1h1; 
  jp[5] = jump[5] * h0h1; 

  double q[2]; 

  q[0] = (xi - xi0) / dx;
  q[1] = (eta - eta0) / dy;

  extern bool computeBoundaryValues(const double q[2], const double jp[6],
                                    const double coord[6][2], const bool s[6], const double grid_u[6], double bdry_u[12]);

  bool status = computeBoundaryValues(q, jp, coord, inside, grid_u, bdry_u);

  if (status == true) {
    bdry_u[1] /= dx;
    bdry_u[2] /= dy; 

    bdry_u[3] /= h0h0; 
    bdry_u[4] /= h1h1;
    bdry_u[5] /= h0h1;

    bdry_u[7] /= dx; 
    bdry_u[8] /= dy;

    bdry_u[9] /= h0h0;
    bdry_u[10] /= h1h1; 
    bdry_u[11] /= h0h1;
  }

  return status;
}


void CartesianGridAndControlPoints::extractDirichletBoundaryData(double **w, double (*F)(double, double), double *bdry_w_vec, int K) 
{
  double j2[3], jump_w[6], bdry_w[12];
  
  double tx, ty, nx, ny;

  for (int k = 0; k < K; k++) {

    double xi = xy[k][0]; 
    double eta = xy[k][1]; 

    double f = F(xi, eta);

    getPoint2(bdry_theta[k], xi, eta, tx, ty, nx, ny); 

    computeJumps(xi, eta, bdry_theta[k], tx, ty, nx, ny, f, j2); 

    jump_w[0] = 0.0;    jump_w[1] = 0.0;     jump_w[2] = 0.0;
    jump_w[3] = j2[0];  jump_w[4] = j2[1];   jump_w[5] = j2[2];

    computeBoundaryValues(w, xi, eta, jump_w, bdry_w);

    bdry_w_vec[k] = bdry_w[0];
  }
}

void CartesianGridAndControlPoints::extractNeumannBoundaryData(double **w, double (*F)(double, double), double *bdry_wn_vec, int K)
{
    double j2[3], jump_w[6], bdry_w[12];

    double tx, ty, nx, ny; 

    for (int k = 0; k < K; k++) {

        double xi = xy[k][0];
        double eta = xy[k][1]; 

        double f = F(xi, eta); 

        getPoint2(bdry_theta[k], xi, eta, tx, ty, nx, ny); 

        computeJumps(xi, eta, bdry_theta[k], tx, ty, nx, ny, f, j2);

        jump_w[0] = 0.0;    jump_w[1] = 0.0;     jump_w[2] = 0.0;
        jump_w[3] = j2[0];  jump_w[4] = j2[1];   jump_w[5] = j2[2]; 

        computeBoundaryValues(w, xi, eta, jump_w, bdry_w);

        double nx = nxny[k][0]; 
        double ny = nxny[k][1]; 

        bdry_wn_vec[k] = bdry_w[1] * nx + bdry_w[2] * ny; 
    }
}

void CartesianGridAndControlPoints::extractDirichletBoundaryData(double **w, double (*F)(double, double, double, double, int, int, double**, double, double), double** source, double *bdry_w_vec, int K) 
{
  double j2[3], jump_w[6], bdry_w[12];
  
  double tx, ty, nx, ny;

  for (int k = 0; k < K; k++) {

    double xi = xy[k][0]; 
    double eta = xy[k][1]; 

    double f = F(x_min, y_min, dx, dy, I, J, source, xi, eta);

    getPoint2(bdry_theta[k], xi, eta, tx, ty, nx, ny); 

    computeJumps(xi, eta, bdry_theta[k], tx, ty, nx, ny, f, j2); 

    jump_w[0] = 0.0;    jump_w[1] = 0.0;     jump_w[2] = 0.0;
    jump_w[3] = j2[0];  jump_w[4] = j2[1];   jump_w[5] = j2[2];

    computeBoundaryValues(w, xi, eta, jump_w, bdry_w);

    bdry_w_vec[k] = bdry_w[0];
  }
}

void CartesianGridAndControlPoints::extractNeumannBoundaryData(double **w, double (*F)(double, double, double, double, int, int, double**, double, double), double** source, double *bdry_wn_vec, int K)
{
    double j2[3], jump_w[6], bdry_w[12];

    double tx, ty, nx, ny; 

    for (int k = 0; k < K; k++) {

        double xi = xy[k][0];
        double eta = xy[k][1]; 

        double f = F(x_min, y_min, dx, dy, I, J, source, xi, eta); 

        getPoint2(bdry_theta[k], xi, eta, tx, ty, nx, ny); 

        computeJumps(xi, eta, bdry_theta[k], tx, ty, nx, ny, f, j2);

        jump_w[0] = 0.0;    jump_w[1] = 0.0;     jump_w[2] = 0.0;
        jump_w[3] = j2[0];  jump_w[4] = j2[1];   jump_w[5] = j2[2]; 

        computeBoundaryValues(w, xi, eta, jump_w, bdry_w);

        double nx = nxny[k][0]; 
        double ny = nxny[k][1]; 

        bdry_wn_vec[k] = bdry_w[1] * nx + bdry_w[2] * ny; 
    }
}