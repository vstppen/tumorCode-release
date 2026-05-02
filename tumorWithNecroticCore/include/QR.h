// Copyright (c) Wenjun Ying, Department of Mathematics and Institute
// of Natural Sciences, Shanghai Jiao Tong University, Shanghai 200240.

// This file is NOT a free software. Without the author's permission,
// please do not copy and/or distribute it to others. 
 
#ifndef __QR_h_IS_INCLUDED__
#define __QR_h_IS_INCLUDED__

#include <cmath> 
#include <iostream> 

#include "Householder.h"

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template <int N>
bool solveByQRdecomposition(double A[N][N], const double rhs[N], 
                            double u[N], int n) 
{
  double b[N]; 
  for (int i = 0; i < n; i++) {
    b[i] = rhs[i]; 
  }

  double v[N];
  double x[N];
  double y[N]; 

  double R[N][N]; 

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      R[i][j] = A[i][j];
    }
  }

  for (int k = 0, m = n; k < n - 1; k++, m--) {
    for (int i = 0; i < m; i++) {
      x[i] = R[i + k][k]; 
    }
    bool status = computeHouseholderVector(x, m, y, v); 
    if (status) {
      for (int i = 0; i < m; i++) {
        R[i + k][k] = y[i]; 
      }

      for (int j = k + 1; j < n; j++) {
        for (int i = 0; i < m; i++) {
          x[i] = R[i + k][j];
        }
        makeHouseholderTransform(v, m, x, y); 
        for (int i = 0; i < m; i++) {
          R[i + k][j] = y[i];
        }
      }

      makeHouseholderTransform(v, m, b + k, y); 
      for (int i = 0; i < m; i++) {
        b[k + i] = y[i];
      }

    } else {
      return false; 
    }
  }

  bool status = true;
  if (fabs(R[n - 1][n - 1]) > EPSILON12) {
    for (int i = n - 1; i >= 0; i--) {
      double sum = 0.0; 
      for (int j = i + 1; j < n; j++) {
        sum += R[i][j] * u[j]; 
      }
      u[i] = (b[i] - sum) / R[i][i]; 
    }
  } else {
    status = false;
  }

  return status;
}

#endif
