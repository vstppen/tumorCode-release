// Copyright (c) Wenjun Ying, Department of Mathematics and Institute
// of Natural Sciences, Shanghai Jiao Tong University, Shanghai 200240.

// This file is NOT a free software. Without the author's permission,
// please do not copy and/or distribute it to others. 
 
#ifndef __Norm_h_IS_INCLUDED__
#define __Norm_h_IS_INCLUDED__ 

#include <cmath>

#include "Const.h"

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template <typename T>
T computeMaxNorm(const T u[], int n)
{
  T r = fabs(u[0]); 
  for (int i = 1; i < n; i++) {
    T e = fabs(u[i]); 
    if (e > r) {
      r = e;
    }
  }
  return r; 
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline double computeEuclidNorm(const double x[], int n)
{
  double sum = 0.0; 
  for (int i = 0; i < n; i++) {
    sum += x[i] * x[i]; 
  }
  return sqrt(sum);
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline double computeVectorNorm(const double x[], int n) 
{
  double sum = 0.0; 
  for (int i = 0; i < n; i++) {
    sum += x[i] * x[i];
  }
  return sqrt(sum); 
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline double computeL2Norm(const double x[], int n)
{
  double sum = 0.0; 
  for (int i = 0; i < n; i++) {
    sum += x[i] * x[i];
  }
  return sqrt(sum);
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline double computeScaledL2Norm(const double x[], int n) 
{
  double sum = 0.0;
  for (int i = 0; i < n; i++) {
    sum += x[i] * x[i]; 
  }
  return sqrt(sum / n);
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline double norm(const double x[], int n)
{
  double sum = 0.0;
  for (int i = 0; i < n; i++) {
    sum += x[i] * x[i];
  }
  return sqrt(sum);
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline int normalize(double x[], int n) 
{
  double sum = 0.0;
  for (int i = 0; i < n; i++) {
    sum += x[i] * x[i];
  }
  double len = sqrt(sum);
  if (len < EPSILON15) {
    return (-1); 
  }

  for (int i = 0; i < n; i++) {
    x[i] /= len; 
  }

  return 0;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline double computeL2Norm(double **const v, int m, int n)
{
  double sum = 0.0; 
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      sum += v[i][j] * v[i][j];
    }
  }
  double r = sqrt(sum);
  return r;
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline double computeScaledL2Norm(double **const v, int m, int n) 
{
  double sum = 0.0; 
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      sum += v[i][j] * v[i][j];
    }
  }
  double r = sqrt(sum / (m * n)); 
  return r; 
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline double computeMaxNorm(double **u, int m, int n)
{
  double max_u = fabs(u[0][0]); 
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      double v = fabs(u[i][j]);
      if (max_u < v) {
        max_u = v;
      }
    }
  }
  return max_u;
}

#endif 
