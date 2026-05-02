#ifndef __LU_h_IS_INCLUDED__
#define __LU_h_IS_INCLUDED__

#include <cmath> 
#include <iostream>
#include "Const.h"

template <int N>
int LUfactorize(double A[N][N], int p[N])
{
  int i, j, k; 
  double m, ajj, lij; 

  // initialize the permutation vector (p).
  for (j = 0; j < N; j++) {
    p[j] = j;
  }

  for (j = 0; j < N; j++) {

    // swap rows to put row with biggest |aij| (i>=j).

    k = j; 
    m = A[p[j]][j]; 
    for (i = j + 1; i < N; i++) {
      if (fabs(A[p[i]][j]) > fabs(m)) {
        m = A[p[i]][j]; 
        k = i;
      }
    }

    if (k != j) {
      int temp = p[j]; 
      p[j] = p[k]; 
      p[k] = temp;
    }

    ajj = A[p[j]][j]; 
    if (fabs(ajj) < EPSILON15) {
      printf("failure in the LU factorization.\n"); 
      return (-1);
    }

    for (i = j + 1; i < N; i++) {
      lij = A[p[i]][j] / ajj; 
      A[p[i]][j] = lij;
      for (k = j + 1; k < N; k++) {
        A[p[i]][k] -= lij * A[p[j]][k]; 
      }
    }
  }

  return 0;

} // LU factorization is A[P[i], j].

/**
 *****************************************************************************
 *****************************************************************************
 **/
template <int N>
void LUsolve_internal(double A[N][N], const int p[N],
                      const double b[N], double x[N])
{
  int i, j, n; 
  double rowsum;

  // forward substitution.
  x[0] = b[p[0]]; 
  for (i = 1; i < N; i++) {
    x[i] = b[p[i]];
    rowsum = 0.0; 
    for (j = 0; j < i; j++) {
      rowsum += A[p[i]][j] * x[j];
    }
    x[i] -= rowsum; 
  }

  n = N; 

  // back substitution.
  x[n - 1] = x[n - 1] / A[p[n - 1]][n - 1]; 
  for (i = n - 2; i >= 0; i--) {
    rowsum = 0.0;
    for (j = n - 1; j > i; j--) {
      rowsum += A[p[i]][j] * x[j]; 
    }
    x[i] = (x[i] - rowsum) / A[p[i]][i]; 
  }
}

/**
 *****************************************************************************
 *****************************************************************************
 **/
template <int N>
bool LUsolve(double A[N][N], const double b[N], double x[N])
{
  double B[N][N]; 
  int i, j, p[N], status; 

  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      B[i][j] = A[i][j];
    }
  }

  status = LUfactorize(B, p);

  if (status != 0) {
    printf("Failed in the LU factorization.\n"); 
    return false;
  }

  LUsolve_internal(B, p, b, x);

  return true;
}

# endif