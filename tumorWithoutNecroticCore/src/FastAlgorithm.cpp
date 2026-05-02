# include <cmath> 
# include <cassert> 

# include "ThomasAlgorithm.h"
# include "LinearAlgebra.h"
# include "FastAlgorithm.h"
# include "LU.h"

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline int Log2(int m) 
{
  int k = 0; 
  while (m > 1) {
    m >>= 1;  // binary operator
    k++; 
  }
  return k;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline int Power2(int m)
{
  assert(m >= 0); 

  int k = 1; 
  for (int i = 0; i < m; i++) { 
    // k *= 2;
    k <<= 1;
  }

  return k;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline bool checkIfPowerOf2(int n)
{
  if (n < 1) {
    return false; 
  }

  int p = n;
  int k = 0;
  while (p > 1) {
    p >>= 1;
    k++; 
  }
  //int k = Log2(n);

  int m = 1;
  for (int i = 0; i < k; i++) {
    m <<= 1; 
  }
  //int m = Power2(k);

  if (m == n) {
    return true; 
  } else {
    return false;
  }
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void computeFastFourierTransform(const double f_r[], 
                                        const double z_r[], 
                                        const double z_i[],
                                        double c_r[], 
                                        double c_i[],
                                        int n)
{
  const int m = Log2(n); 

  for (int i = 0; i < n; i++) {
    c_r[i] = f_r[i];
    c_i[i] = 0.0;
  }

  int j = 0;
  int i2 = n >> 1;
  for (int i = 0; i < n - 1; i++) {
    if (i < j) {
      double tr = c_r[i];
      double ti = c_i[i]; 
      c_r[i] = c_r[j]; 
      c_i[i] = c_i[j]; 
      c_r[j] = tr;
      c_i[j] = ti;
    }
    int k = i2;
    while (k <= j) {
      j -= k; 
      k >>= 1; 
    }
    j += k; 
  }

  int q = 1;
  int p = q << 1; 
  int s = n >> 1; 

  {
    for (int i = 0; i < n; i += p) {
      int j0 = i, j1 = i + q; 
      double tmp_r = c_r[j1];
      c_r[j1] = c_r[j0] - tmp_r;
      c_r[j0] += tmp_r; 
      j0++;  j1++; 
      for (int r = 1, t = s; r < q; r++, j0++, j1++, t += s) {
        double tmp_r = c_r[j1] * z_r[t]; 
        double tmp_i = c_r[j1] * z_i[t];
        c_r[j1] = c_r[j0] - tmp_r; 
        c_r[j0] += tmp_r;
      }
    }
    q = p;  p = q << 1;  s = s >> 1;
  }

  for (int k = 1; k < m; k++) {
    for (int i = 0; i < n; i += p) {
      int j0 = i, j1 = i + q; 
      double tmp_r = c_r[j1]; 
      double tmp_i = c_i[j1]; 
      c_r[j1] = c_r[j0] - tmp_r; 
      c_i[j1] = c_i[j0] - tmp_i; 
      c_r[j0] += tmp_r; 
      c_i[j0] += tmp_i; 
      j0++;  j1++;
      for (int r = 1, t = s; r < q; r++, j0++, j1++, t += s) {
        double tmp_r = c_r[j1] * z_r[t] - c_i[j1] * z_i[t]; 
        double tmp_i = c_r[j1] * z_i[t] + c_i[j1] * z_r[t]; 
        c_r[j1] = c_r[j0] - tmp_r; 
        c_i[j1] = c_i[j0] - tmp_i;
        c_r[j0] += tmp_r; 
        c_i[j0] += tmp_i;
      }
    }
    q = p;  p = q << 1;  s = s >> 1; 
  }
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void computeFastFourierTransform(const double f[], double c_r[], 
                                        double c_i[], int n) 
{
  const int n2 = n >> 1;

  const double h = 2.0 * M_PI / n; 

  double *f_r = new double[n]; 
  double *f_i = new double[n]; 

  double *z_r = new double[n2]; 
  double *z_i = new double[n2]; 

  double x = 0.0;
  for (int i = 0; i < n; i++) {
    f_r[i] = f[i]; 
    f_i[i] = 0.0;
    x += h; 
  }

  const double alpha = h;

  const double unit_r = cos(alpha);
  const double unit_i = - sin(alpha); 

  double tr = 1.0, ti = 0.0; 
  z_r[0] = 1.0;  z_i[0] = 0.0;
  for (int i = 1; i < n2; i++) {
    double sr = tr * unit_r - ti * unit_i; 
    double si = ti * unit_r + tr * unit_i; 
    z_r[i] = sr;  z_i[i] = si; 
    tr = sr;  ti = si;
  }

  computeFastFourierTransform(f_r, z_r, z_i, c_r, c_i, n); 

  for (int i = 0, j = n2; i < n2; i++, j++) {
    c_r[i] /= n;
    c_i[i] /= n;
    c_r[j] /= n; 
    c_i[j] /= n; 
  }

  delete[] z_r; 
  delete[] z_i; 

  delete[] f_r; 
  delete[] f_i; 
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void computeFastSineTransform(const double u[], double c[], int m)
{
  int n = 0;
  bool is_even = true; 
  if (m%2 == 0) {
    n = m << 1; 
  } else {
    is_even = false;
    n = (m + 1) << 1;
  }
  int k = Log2(n);

  double *f_r = new double[n];
  double *c_r = new double[n]; 
  double *c_i = new double[n];

  if (is_even) {
    for (int i = 0; i < m; i++) {
      f_r[i] = u[i];
    }
    f_r[0] = 0.0;
    for (int i = m; i < n; i++) {
      f_r[i] = 0.0; 
    }
  } else {
    int j = 0; 
    f_r[j++] = 0.0;
    for (int i = 0; i < m; i++, j++) {
      f_r[j] = u[i]; 
    }
    for (int i = j; i < n; i++) {
      f_r[i] = 0.0; 
    }
  }

  computeFastFourierTransform(f_r, c_r, c_i, n);

  const double factor = 2.0 * sqrt(n); 

  if (is_even) {
    for (int i = 1; i < m; i++) {
      c[i] = c_i[i] * factor; 
    }
  } else {
    for (int i = 0, j = 1; i < m; i++, j++) {
      c[i] = c_i[j] * factor;
    }
  }

  delete[] f_r;
  delete[] c_r; 
  delete[] c_i; 
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void computeFastCosineTransform(const double u[], double c[], int m)
{
  int m1 = m - 1; 
  int n = m1 << 1; 
  int k = Log2(n); 

  double *f_r = new double[n];
  double *c_r = new double[n];
  double *c_i = new double[n]; 

  int m12 = m1 >> 1;

  for (int i = 1; i < m1; i++) {
    f_r[i] = u[i] / m12; 
  }
  f_r[m1] = u[m1] / m1;
  f_r[0] = u[0] / m1; 

  for (int i = m; i < n; i++) {
    f_r[i] = 0.0;
  }

  computeFastFourierTransform(f_r, c_r, c_i, n);

  for (int i = 0; i < m; i++) {
    c[i] = c_r[i] * n;
  }

  const double gamma = 0.5 * sqrt(2.0); 

  c[m1] *= gamma; 
  c[0] *= gamma;

  delete[] f_r;
  delete[] c_r; 
  delete[] c_i; 
}

inline void computeFastInverseCosineTransform(const double u[], 
                                              double c[], int m)
{
  int m1 = m - 1; 
  int n = m1 << 1;
  int k = Log2(n); 

  double *f_r = new double[n];
  double *c_r = new double[n]; 
  double *c_i = new double[n];

  const double gamma = 0.5 * sqrt(2.0); 

  for (int i = 0; i < m; i++) {
    f_r[i] = u[i];
  }
  f_r[0] *= gamma;
  f_r[m1] *= gamma;
  for (int i = m; i < n; i++) {
    f_r[i] = 0.0;
  }

  computeFastFourierTransform(f_r, c_r, c_i, n);

  for (int i = 0; i < m; i++) {
    c[i] = c_r[i] * n; 
  }

  delete[] f_r;
  delete[] c_r; 
  delete[] c_i;
}

/**
 *****************************************************************************
 * The fast solver below solves the Dirichlet BVP of the Poisson equation:   *                                                              *
 *                 - \laplace u = f(x, y)                         *
 * on a box with the five-point node-centered finite difference method.      *
 *                                                                           *
 * Note: Currently, we require the grid has equal spacings along both the    *
 * horizontal and the vertical directions.                                   *
 *****************************************************************************
 **/
void FastPoissonSolver(double **unknown, const int sz[2])
{
  int size[2]; 

  size[0] = sz[0] - 1;
  size[1] = sz[1] - 1;

  {
    double *dst = new double[size[1]]; 
    double *src = new double[size[1]];
    for (int i = 1; i < size[0]; i++) {
      for (int j = 1; j < size[1]; j++) {
        src[j] = unknown[i][j]; 
      }
      src[0] = 0.0;

      computeFastSineTransform(src, dst, size[1]); 

      for (int j = 1; j < size[1]; j++) {
        unknown[i][j] = dst[j];
      }
    }
    delete[] src;
    delete[] dst;
  }

  {
    double *src = new double[size[0]]; 
    double *dst = new double[size[0]];

    double *alpha = new double[size[0] - 1];
    double *gamma = new double[size[0] - 1];
    double *beta = new double[size[0] - 1];

    for (int j = 1; j < size[1]; j++) {

      for (int i = 1; i < size[0]; i++) {
        src[i] = unknown[i][j];
      }

      double lambda_j = 2.0 * (1.0 - cos(j * M_PI / size[1]));

      const double diag = 2.0 + lambda_j; 
      for (int i = 0; i < size[0] - 1; i++) {
        beta[i] = gamma[i] = - 1.0; 
        alpha[i] = diag;
      }

      LU(alpha, beta, gamma, src + 1, size[0] - 1, dst + 1);

      for (int i = 1; i < size[0]; i++) {
        unknown[i][j] = dst[i];
      }
    }

    delete[] gamma;
    delete[] alpha; 
    delete[] beta; 

    delete[] src; 
    delete[] dst;
  }

  {
    double *dst = new double[size[1]]; 
    double *src = new double[size[1]];
    for (int i = 1; i < size[0]; i++) {
      for (int j = 1; j < size[1]; j++) {
        src[j] = unknown[i][j]; 
      }
      src[0] = 0.0;
      computeFastSineTransform(src, dst, size[1]); 
      for (int j = 1; j < size[1]; j++) {
        unknown[i][j] = dst[j]; 
      }
    }
    delete[] src; 
    delete[] dst;
  }
}

/**
 *****************************************************************************
 * This one solves the Dirichlet BVP of the anisotropic diffusion equation   *
 *         - \sigma[0] u_xx - \sigma[1] u_yy = f(x, y)                       *
 * on a box with the five-point node-centered finite difference method.      *
 *                                                                           *
 * Note: Currently, we require the grid has equal spacings along both the    *
 * horizontal and the vertical directions.                                   *
 *****************************************************************************
 **/
void FastPoissonSolver(const double sigma[2], double **unknown,
                       const int sz[2]) // node-centered
{
  int size[2]; 

  size[0] = sz[0] - 1;
  size[1] = sz[1] - 1;

  {
    /* Transform the rows */
    double *dst = new double[size[1]];
    double *src = new double[size[1]];
    for (int i = 1; i < size[0]; i++) {
      for (int j = 1; j < size[1]; j++) {
        src[j] = unknown[i][j]; 
      }
      src[0] = 0.0; 
      computeFastSineTransform(src, dst, size[1]); 
      for (int j = 1; j < size[1]; j++) {
        unknown[i][j] = dst[j];
      }
    }
    delete[] src;
    delete[] dst; 
  }

  {
    double *src = new double[size[0]]; 
    double *dst = new double[size[0]]; 

    double *alpha = new double[size[0] - 1];
    double *gamma = new double[size[0] - 1]; 
    double *beta = new double[size[0] - 1];

    for (int j = 1; j < size[1]; j++) {
      for (int i = 1; i < size[0]; i++) {
        src[i] = unknown[i][j]; 
      }

      double lambda_j = 2.0 * (1.0 - cos(j * M_PI / size[1])); 

//#ifdef USE_TRIDIAGONAL_SOLVER

      const double diag = 2.0 * sigma[0] + lambda_j * sigma[1];

      for (int i = 0; i < size[0] - 1; i++) {
        beta[i] = gamma[i] = - sigma[0]; 
        alpha[i] = diag; 
      }

      // Use Tridiagonal System Solver.
      LU(alpha, beta, gamma, src + 1, size[0] - 1, dst + 1); 

//#else

      /* double theta_i = M_PI / size[0];

      const double c = lambda_j * sigma[1];

      computeFastSineTransform(src, dst, size[0]);

      for (int i = 0; i < size[0]; i++) {
        double lambda_i = 2.0 * (1.0 - cos(i * theta_i));
        dst[i] /= lambda_i * sigma[0] + c;
      }

      computeFastInverseSineTransform(dst, dst, size[0]); */

//#endif

      for (int i = 1; i < size[0]; i++) {
        unknown[i][j] = dst[i];
      }
    }

    delete[] gamma;
    delete[] alpha;
    delete[] beta; 

    delete[] src;
    delete[] dst;
  }

  {
    /* Transform the rows */
    double *dst = new double[size[1]]; 
    double *src = new double[size[1]];
    for (int i = 1; i < size[0]; i++) {
      for (int j = 1; j < size[1]; j++) {
        src[j] = unknown[i][j];
      }
      src[0] = 0.0; 
      computeFastSineTransform(src, dst, size[1]); 
      for (int j = 1; j < size[1]; j++) {
        unknown[i][j] = dst[j]; 
      }
    }
    delete[] src;
    delete[] dst;
  }
}


/**
 *****************************************************************************
 * The fast solver below solves the Dirichlet BVP of the modified Helmholtz  *
 * equation:                                                                 *
 *       - \sigma[0] u_xx - \sigma[1] u_yy + \kappa u = f(x, y)              *
 * on a box with the five-point node-centered finite difference method.      *
 *                                                                           *
 * Note: Currently, we require the grid has equal spacings along both the    *
 * horizontal and the vertical directions.                                   *
 *****************************************************************************
 **/
void FastMHelmholtzSolver(const double sigma[2], double **unknown, const int sz[2], double eta) // node-centered
{
  int size[2]; 

  size[0] = sz[0] - 1;
  size[1] = sz[1] - 1; 

  {
    /* Transform the rows */
    double *dst = new double[size[1]]; 
    double *src = new double[size[1]];
    for (int i = 1; i < size[0]; i++) {
      for (int j = 1; j < size[1]; j++) {
        src[j] = unknown[i][j];
      }
      src[0] = 0.0;
      computeFastSineTransform(src, dst, size[1]);
      for (int j = 1; j < size[1]; j++) {
        unknown[i][j] = dst[j]; 
      }
    }
    delete[] src;
    delete[] dst; 
  }

  {
    double *src = new double[size[0]];
    double *dst = new double[size[0]];

    double *alpha = new double[size[0] - 1];
    double *gamma = new double[size[0] - 1]; 
    double *beta = new double[size[0] - 1];

    for (int j = 1; j < size[1]; j++) {
      for (int i = 1; i < size[0]; i++) {
        src[i] = unknown[i][j]; 
      }

      double lambda_j = 2.0 * (1.0 - cos(j * M_PI / size[1])); 

      const double diag = 2.0 * sigma[0] + lambda_j * sigma[1] + eta;

      for (int i = 0; i < size[0] - 1; i++) {
        beta[i] = gamma[i] = - sigma[0];
        alpha[i] = diag; 
      }

      LU(alpha, beta, gamma, src + 1, size[0] - 1, dst + 1);

      for (int i = 1; i < size[0]; i++) {
        unknown[i][j] = dst[i];
      }
    }

    delete[] gamma; 
    delete[] alpha; 
    delete[] beta; 

    delete[] src; 
    delete[] dst;
  }

  {
    /* Transform the rows */
    double *dst = new double[size[1]];
    double *src = new double[size[1]]; 
    for (int i = 1; i < size[0]; i++) {
      for (int j = 1; j < size[1]; j++) {
        src[j] = unknown[i][j]; 
      }
      src[0] = 0.0;
      computeFastSineTransform(src, dst, size[1]); 
      for (int j = 1; j < size[1]; j++) {
        unknown[i][j] = dst[j]; 
      }
    }
    delete[] src;
    delete[] dst; 
  }
}


/**
 *****************************************************************************
 * The fast solver below solves the Dirichlet BVP of the modified Helmholtz  *
 * equation:                                                                 *
 *      \lambda u - u_xx - u_yy = f(x, y)                                    *
 * on a box with the five-point node-centered finite difference method.      *
 *                                                                           *
 * Note: Currently, we require the grid has equal spacings along both the    *
 * horizontal and the vertical directions.                                   *
 *****************************************************************************
 **/
void FastMHelmholtzSolver(double **unknown, const int sz[2], double eta) 
{
  double sigma[2] = {1., 1.};
  FastMHelmholtzSolver(sigma, unknown, sz, eta);
}