# ifndef CUBICSPLINE_H
# define CUBICSPLINE_H

# include "LinearAlgebra.h"
# include "Const.h"

inline void periodicCubicSplineInterpolation(double* f, double* M, double* alpha, double* beta, int n)
{
    double h = M_2PI / n;
    double h2 = h * h;

    double** A = createMatrix(n, n);
    double* g = createVector(n);

    A[0][0] = 2. / 3.;
    A[0][1] = 1. / 6.;
    A[0][n - 1] = 1. / 6.;
    for (int i = 1; i < n - 1; i++) {
        A[i][i - 1] = 1. / 6.;
        A[i][i] = 2. / 3.;
        A[i][i + 1] = 1. / 6.;
    }
    A[n - 1][0] = 1. / 6.;
    A[n - 1][n - 2] = 1. / 6.;
    A[n - 1][n - 1] = 2. / 3.;

    g[0] =  (f[n - 1] - 2. * f[0] + f[1]) / h2;
    for (int i = 1; i < n - 1; i++) {
        g[i] = (f[i - 1] - 2. * f[i] + f[i + 1]) / h2;
    }
    g[n - 1] = (f[n - 2] - 2. * f[n - 1] + f[0]) / h2;

    // printMatrix(A, n, n);
    SolveLinearSystemByCG(A, g, M, n, 0);

    for (int i = 0; i < n; i++) {
        alpha[i] = f[i] - 1. / 6. * h2 * M[i];
        beta[i] = f[(i + 1) % n] - 1. / 6. * h2 * M[(i + 1) % n];
    }

    // printVector(M, n);
    // printVector(alpha, n);
    // printVector(beta, n);

    freeVector(g);
    freeMatrix(A, n);
}

inline void periodicCubicSplineCurveInterpolation(double** coor, double* M_x, double* alpha_x, double* beta_x, double* M_y, double* alpha_y, double* beta_y, int n)
{
    double* x = createVector(n);
    double* y = createVector(n);

    for (int i = 0; i < n; i++) {
        x[i] = coor[i][0];
        y[i] = coor[i][1];
    }

    periodicCubicSplineInterpolation(x, M_x, alpha_x, beta_x, n);
    periodicCubicSplineInterpolation(y, M_y, alpha_y, beta_y, n);

    freeVector(x);
    freeVector(y);
}


inline double periodicCubicSplineGetS(double* M, double* alpha, double* beta, int n, double x)
{   
    double h = M_2PI / n;
    x = std::fmod(x, M_2PI);
    if (x < 0) {
        x += 2 * M_PI;
    }

    double s;
    
    for (int i = 0; i < n; i++) {
        double x_start = i * h;
        double x_end = (i + 1) * h;
        if (x > x_start - EPSILON12 & x < x_end) {
            s = M[i] * (x_end - x) * (x_end - x) * (x_end - x) / (6. * h);
            s += M[(i + 1) % n] * (x - x_start) * (x - x_start) * (x - x_start) / (6. * h);
            s += alpha[i] * (x_end - x) / h + beta[i] * (x - x_start) / h;
        }
    }

    return s;
}

inline double periodicCubicSplineGetDS(double* M, double* alpha, double* beta, int n, double x)
{   
    double h = M_2PI / n;
    x = std::fmod(x, M_2PI);
    if (x < 0) {
        x += 2 * M_PI;
    }

    double s;
    
    for (int i = 0; i < n; i++) {
        double x_start = i * h;
        double x_end = (i + 1) * h;
        if (x > x_start - EPSILON12 & x < x_end) {
            s = - M[i] * (x_end - x) * (x_end - x) / (2. * h);
            s += M[(i + 1) % n] * (x - x_start) * (x - x_start) / (2. * h);
            s += (beta[i] - alpha[i]) / h;
        }
    }

    return s;
}

inline double periodicCubicSplineGetDDS(double* M, double* alpha, double* beta, int n, double x)
{   
    double h = M_2PI / n;
    x = std::fmod(x, M_2PI);
    if (x < 0) {
        x += 2 * M_PI;
    }
    
    double s;
    
    for (int i = 0; i < n; i++) {
        double x_start = i * h;
        double x_end = (i + 1) * h;
        if (x > x_start - EPSILON12 & x < x_end) {
            s = M[i] * (x_end - x) / h + M[(i + 1) % n] * (x - x_start) / h;
        }
    }

    return s;
}


inline double computeTotalLength(int n, double* M_x, double* alpha_x, double* beta_x, double* M_y, double* alpha_y, double* beta_y)
{
    double h = M_2PI / n;

    double result = 0.;

    double a, b, t1, t2, t3, w1, w2, w3, f1, f2, f3;
    for (int i = 0; i < n; i++) {
        a = i * h;
        b = (i + 1) * h;
        t1 = (a + b) / 2. - (b - a) / 2. * sqrt(0.6);
        t2 = (a + b) / 2.;
        t3 = (a + b) / 2. + (b - a) / 2. * sqrt(0.6);
        w1 = 5 * (b - a) / 18.;
        w2 = 4 * (b - a) / 9.;
        w3 = w1;
        f1 = sqrt(pow(periodicCubicSplineGetDS(M_x, alpha_x, beta_x, n, t1), 2) + pow(periodicCubicSplineGetDS(M_y, alpha_y, beta_y, n, t1), 2));
        f2 = sqrt(pow(periodicCubicSplineGetDS(M_x, alpha_x, beta_x, n, t2), 2) + pow(periodicCubicSplineGetDS(M_y, alpha_y, beta_y, n, t2), 2));
        f3 = sqrt(pow(periodicCubicSplineGetDS(M_x, alpha_x, beta_x, n, t3), 2) + pow(periodicCubicSplineGetDS(M_y, alpha_y, beta_y, n, t3), 2));
        result += w1 * f1 + w2 * f2 + w3 * f3;
    }

    return result;
}

inline double computeLength(int n, double* M_x, double* alpha_x, double* beta_x, double* M_y, double* alpha_y, double* beta_y, double t)
{
    double h = M_2PI / n;

    double result = 0.;
    int idx = static_cast<int>(t / h);

    double a, b, t1, t2, t3, w1, w2, w3, f1, f2, f3;
    for (int i = 0; i < idx; i++) {
        a = i * h;
        b = (i + 1) * h;
        t1 = (a + b) / 2. - (b - a) / 2. * sqrt(0.6);
        t2 = (a + b) / 2.;
        t3 = (a + b) / 2. + (b - a) / 2. * sqrt(0.6);
        w1 = 5 * (b - a) / 18.;
        w2 = 4 * (b - a) / 9.;
        w3 = w1;
        f1 = sqrt(pow(periodicCubicSplineGetDS(M_x, alpha_x, beta_x, n, t1), 2) + pow(periodicCubicSplineGetDS(M_y, alpha_y, beta_y, n, t1), 2));
        f2 = sqrt(pow(periodicCubicSplineGetDS(M_x, alpha_x, beta_x, n, t2), 2) + pow(periodicCubicSplineGetDS(M_y, alpha_y, beta_y, n, t2), 2));
        f3 = sqrt(pow(periodicCubicSplineGetDS(M_x, alpha_x, beta_x, n, t3), 2) + pow(periodicCubicSplineGetDS(M_y, alpha_y, beta_y, n, t3), 2));
        result += w1 * f1 + w2 * f2 + w3 * f3;
    }

    a = idx * h;
    b = t;
    t1 = (a + b) / 2. - (b - a) / 2. * sqrt(0.6);
    t2 = (a + b) / 2.;
    t3 = (a + b) / 2. + (b - a) / 2. * sqrt(0.6);
    w1 = 5 * (b - a) / 18.;
    w2 = 4 * (b - a) / 9.;
    w3 = w1;
    f1 = sqrt(pow(periodicCubicSplineGetDS(M_x, alpha_x, beta_x, n, t1), 2) + pow(periodicCubicSplineGetDS(M_y, alpha_y, beta_y, n, t1), 2));
    f2 = sqrt(pow(periodicCubicSplineGetDS(M_x, alpha_x, beta_x, n, t2), 2) + pow(periodicCubicSplineGetDS(M_y, alpha_y, beta_y, n, t2), 2));
    f3 = sqrt(pow(periodicCubicSplineGetDS(M_x, alpha_x, beta_x, n, t3), 2) + pow(periodicCubicSplineGetDS(M_y, alpha_y, beta_y, n, t3), 2));
    result += w1 * f1 + w2 * f2 + w3 * f3;

    return result;
}

inline double findParameterByArcLength(int n, double* M_x, double* alpha_x, double* beta_x, double* M_y, double* alpha_y, double* beta_y, double target_length) {
    const double eps = 1e-12;
    double low = 0.0, high = 2 * M_PI;
    while (high - low > eps) {
        double mid = 0.5 * (low + high);
        double current_length = computeLength(n, M_x, alpha_x, beta_x, M_y, alpha_y, beta_y, mid);
        if (current_length < target_length) {
            low = mid;
        } else {
            high = mid;
        }
    }
    return 0.5 * (low + high);
}

inline void SelectPointsUniformly(int n, double* M_x, double* alpha_x, double* beta_x, double* M_y, double* alpha_y, double* beta_y, double** new_points) {
    const double total_length = computeTotalLength(n, M_x, alpha_x, beta_x, M_y, alpha_y, beta_y);
    const double segment_length = total_length / n;

    for (int i = 0; i < n; i++) {
        const double target_length = i * segment_length;
        const double t = findParameterByArcLength(n, M_x, alpha_x, beta_x, M_y, alpha_y, beta_y, target_length);
        new_points[i][0] = periodicCubicSplineGetS(M_x, alpha_x, beta_x, n, t);
        new_points[i][1] = periodicCubicSplineGetS(M_y, alpha_y, beta_y, n, t);
    }
}

inline void SelectPointsUniformly(int n, double* M_x, double* alpha_x, double* beta_x, double* M_y, double* alpha_y, double* beta_y, double** new_points, int n_new) {
    const double total_length = computeTotalLength(n, M_x, alpha_x, beta_x, M_y, alpha_y, beta_y);
    const double segment_length = total_length / n_new;

    for (int i = 0; i < n_new; i++) {
        const double target_length = i * segment_length;
        const double t = findParameterByArcLength(n, M_x, alpha_x, beta_x, M_y, alpha_y, beta_y, target_length);
        new_points[i][0] = periodicCubicSplineGetS(M_x, alpha_x, beta_x, n, t);
        new_points[i][1] = periodicCubicSplineGetS(M_y, alpha_y, beta_y, n, t);
    }
}

# endif