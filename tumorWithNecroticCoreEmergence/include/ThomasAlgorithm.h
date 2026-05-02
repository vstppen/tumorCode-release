# ifndef THOMASALGORITHM_H
# define THOMASALGORITHM_H

#include "TangentialDerivatives.h"
#include "Const.h"

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
int LU(const double a[], const double b[], const double c[], const double rhs[], int m, double u[]) 
{
  double *beta = new double[m];
  double *alpha = new double[m];
  double *gamma = new double[m];

  for (int i = 0; i < m; i++) {
    beta[i] = b[i]; 
  }

  alpha[0] = a[0];
  gamma[0] = c[0] / alpha[0]; 
  for (int i = 1; i < m - 1; i++) {
    alpha[i] = a[i] - beta[i] * gamma[i - 1]; 
    if (fabs(alpha[i]) < EPSILON) {
      delete[] alpha; 
      delete[] gamma;
      delete[] beta;
      return (-1); 
    }
    gamma[i] = c[i] / alpha[i];
  }
  alpha[m - 1] = a[m - 1] - beta[m - 1] * gamma[m - 2];
  if (fabs(alpha[m - 1]) < EPSILON) {
    delete[] alpha; 
    delete[] gamma;
    delete[] beta;
    return (-1); 
  }

  double *z = new double[m];

  z[0] = rhs[0] / alpha[0];
  for (int i = 1; i < m; i++) {
    z[i] = (rhs[i] - beta[i] * z[i - 1]) / alpha[i]; 
  }

  u[m - 1] = z[m - 1];
  for (int i = m - 2; i >= 0; i--) {
    u[i] = z[i] - gamma[i] * u[i + 1];
  }

  delete[] z;
  delete[] alpha;
  delete[] gamma; 
  delete[] beta;

  return 0; 
}

# endif