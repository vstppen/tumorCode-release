// Copyright (c) Wenjun Ying, Department of Mathematics and Institute
// of Natural Sciences, Shanghai Jiao Tong University, Shanghai 200240.

// This file is NOT a free software. Without the author's permission,
// please do not copy and/or distribute it to others. 
 
#ifndef __Product_h_IS_INCLUDED__
#define __Product_h_IS_INCLUDED__ 

/**
 *****************************************************************************
 *****************************************************************************
 **/
template <typename T>
void CrossProduct(const T u[3], const T v[3], T w[3])
{
  w[0] = u[1] * v[2] - u[2] * v[1];
  w[1] = u[2] * v[0] - u[0] * v[2];
  w[2] = u[0] * v[1] - u[1] * v[0]; 
}

/**
 *****************************************************************************
 *****************************************************************************
 **/
template <typename T>
void computeCrossProduct(const T u[3], const T v[3], T w[3])
{
  w[0] = u[1] * v[2] - u[2] * v[1]; 
  w[1] = u[2] * v[0] - u[0] * v[2];
  w[2] = u[0] * v[1] - u[1] * v[0];
}

/**
 *****************************************************************************
 *****************************************************************************
 **/
template <typename T>
T computeInnerProduct(const T a[], const T b[], int n) 
{
  T r = 0;
  for (int i = 0; i < n; i++) {
    r += a[i] * b[i]; 
  }
  return r;
}

/**
 *****************************************************************************
 *****************************************************************************
 **/
inline void multiply(double **A, double *x, int n, double *b)
{
  for (int i = 0; i < n; i++) {
    double sum = 0.0;
    for (int j = 0; j < n; j++) {
      sum += A[i][j] * x[j]; 
    }
    b[i] = sum; 
  }
}

/**
 *****************************************************************************
 *****************************************************************************
 **/
template <typename T>
void computeMatrixVectorProduct(T **const A, const T u[], T w[], int n) 
{
  for (int i = 0; i < n; i++) {
    T s = 0.0; 
    for (int j = 0; j < n; j++) {
      s += A[i][j] * u[j]; 
    }
    w[i] = s; 
  }
}

/**
 *****************************************************************************
 *****************************************************************************
 **/
inline double computeInnerProduct(double **const v,
                                  double **const w, 
                                  int m, int n) 
{
  double sum = 0.0;
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      sum += v[i][j] * w[i][j]; 
    }
  }
  return sum; 
}

#endif 
