// Copyright (c) Wenjun Ying, Department of Mathematics and Institute
// of Natural Sciences, Shanghai Jiao Tong University, Shanghai 200240.

// This file is NOT a free software. Without the author's permission,
// please do not copy and/or distribute it to others. 
 
#ifndef __Householder_h_IS_INCLUDED__
#define __Householder_h_IS_INCLUDED__ 

#include "LinearAlgebra.h" 
#include "Norm.h" 
#include "Product.h"

/**
 *****************************************************************************
 *****************************************************************************
 **/
// y = (|x|, 0, ..., 0) or (-|x|, 0, ..., 0) and v = (x - y) / |x - y|
inline bool computeHouseholderVector(const double *x, int m, double *y, double *v)  
{
  for (int i = 0; i < m; i++) {
    y[i] = 0.0;
  }

  double norm_x = computeL2Norm(x, m); 
  if (norm_x / sqrt(m) < EPSILON12) {
    for (int i = 0; i < m; i++) {
      y[i] = x[i]; 
      v[i] = 0.0; 
    }
    return false; 
  }

  if (x[0] > EPSILON12) {
    y[0] = - norm_x;
  } else {
    y[0] = norm_x;
  }

  for (int i = 0; i < m; i++) {
    v[i] = x[i];
  }
  v[0] -= y[0];

  double norm_v = sqrt(norm_x * norm_x - x[0] * x[0] + v[0] * v[0]); 
  if (norm_v < EPSILON12) {
    for (int i = 0; i < m; i++) {
      v[i] = 0.0;
    }
  } else {
    double r_norm_v = 1.0 / norm_v;
    for (int i = 0; i < m; i++) {
      v[i] *= r_norm_v;
    }
  }

  return true; 
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
// b = (|a|, 0, ..., 0) or (-|a|, 0, ..., 0) and w = (a - b) / |a - b| 
inline bool computeHouseholderVector(const double a[], double w[], int n) 
{
  double len = computeL2Norm(a, n); 
  if (len / sqrt(n) < EPSILON12) {
    for (int i = 0; i < n; i++) {
      w[i] = 0.0; 
    }
    return false;
  }

  if (a[0] > 0.0) {
    w[0] = a[0] + len;
  } else {
    w[0] = a[0] - len; 
  }
  for (int i = 1; i < n; i++) {
    w[i] = a[i]; 
  }

  double len2 = computeL2Norm(w, n);
  if (len2 > EPSILON) {
    for (int i = 0; i < n; i++) {
      w[i] /= len2; 
    }
  }

  return true; 
}

/**
 *****************************************************************************
 * w: the Householder vector                                                 *
 * n: the dimension of the vectors                                           *
 * u: the vector before reflection                                           *
 * v: the vector after reflection                                            *
 * v = (I - 2 w w^T) * u                                                     *
 *****************************************************************************
 **/
inline void makeHouseholderTransform(const double w[], const double u[],
                                     int n, double v[])
{
  double p = computeInnerProduct(w, u, n); 
  for (int i = 0; i < n; i++) {
    double c = p * w[i]; 
    v[i] = u[i] - (c + c);
  }
}

/**
 *****************************************************************************
 * w: the Householder vector                                                 *
 * n: the dimension of the vectors                                           *
 * u: the vector before reflection                                           *
 * v: the vector after reflection                                            *
 * v = (I - 2 w w^T) * u                                                     *
 *****************************************************************************
 **/
inline void makeHouseholderTransform(const double w[], int n, 
                                     const double u[], double v[]) 
{
  double p = computeInnerProduct(w, u, n);
  for (int i = 0; i < n; i++) {
    double c = p * w[i]; 
    v[i] = u[i] - (c + c); 
  }
}

#endif
