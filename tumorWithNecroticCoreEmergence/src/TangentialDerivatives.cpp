// #include "LinearAlgebra.h"
#include "QR.h"
#include "TangentialDerivatives.h"

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void computeFirstDerivative(double delta,
                            const double *theta, 
                            const double *psi, int K,
                            double s, double &w0, double &w1) 
{
  int k = static_cast<int>(0.5 + s / delta); 

  int stencil[5];

  stencil[0] = k - 2; 
  stencil[1] = k - 1;
  stencil[2] = k;
  stencil[3] = k + 1; 
  stencil[4] = k + 2; 

  if ((k < 2) || (k >= (K - 2))) {
    if (k == 0) {
      stencil[0] = K - 2;
      stencil[1] = K - 1; 
      stencil[2] = 0;
      stencil[3] = 1; 
      stencil[4] = 2;
    } else if (k == 1) {
      stencil[0] = K - 1;
      stencil[1] = 0;
      stencil[2] = 1;
      stencil[3] = 2; 
      stencil[4] = 3;
    } else if (k == (K - 2)) {
      stencil[0] = K - 4;
      stencil[1] = K - 3;
      stencil[2] = K - 2; 
      stencil[3] = K - 1;
      stencil[4] = 0;
    } else if (k == (K - 1)) {
      stencil[0] = K - 3; 
      stencil[1] = K - 2;
      stencil[2] = K - 1;
      stencil[3] = 0;
      stencil[4] = 1;
    } else if (k == K) {
      stencil[0] = K - 2;
      stencil[1] = K - 1; 
      stencil[2] = 0; 
      stencil[3] = 1; 
      stencil[4] = 2; 
    }
  }

  if (s > theta[k]) {
    for (int i = 0; i < 4; i++) {
      stencil[i] = stencil[i + 1];
    }
  }

  double t[5];
  for (int i = 0; i < 4; i++) {
    t[i] = theta[stencil[i]]; 
  }
  t[4] = s; 

  int k_p = 0;
  int k_n = 0;
  for (int i = 0; i < 5; i++) {
    if (t[i] < M_PI_4) {
      k_p++; 
    }
    if (t[i] > (M_2PI - M_PI_4)) {
      k_n++; 
    }
  }
  if ((k_p > 0) && (k_n > 0)) {
    for (int i = 0; i < 5; i++) {
      if (t[i] > (M_2PI - M_PI_4)) {
        t[i] -= M_2PI; 
      }
    }
  }

  for (int i = 0; i < 4; i++) {
    t[i] -= t[4];
  }

  double mat[4][4];

  for (int i = 0; i < 4; i++) {
    mat[i][0] = 1.0;
    for (int j = 1; j < 4; j++) {
      mat[i][j] = mat[i][j - 1] * t[i];
    }
  }

  double a[4], b[4]; 

  for (int i = 0; i < 4; i++) {
    b[i] = psi[stencil[i]];
  }

  bool status = solveByQRdecomposition<4>(mat, b, a, 4); 

  w0 = a[0];  w1 = a[1]; 
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void computeFirstAndSecondDerivatives(double delta, 
                                      const double *theta,
                                      const double *phi, int K, 
                                      double s, double &w0, double &w1,
                                      double &w2)
{
  int k = static_cast<int>(0.5 + s / delta); 

  int stencil[5]; 

  stencil[0] = k - 2;
  stencil[1] = k - 1;
  stencil[2] = k;
  stencil[3] = k + 1;
  stencil[4] = k + 2;

  if ((k < 2) || (k >= (K - 2))) {
    if (k == 0) {
      stencil[0] = K - 2; 
      stencil[1] = K - 1;
      stencil[2] = 0;
      stencil[3] = 1;
      stencil[4] = 2;
    } else if (k == 1) {
      stencil[0] = K - 1;
      stencil[1] = 0; 
      stencil[2] = 1; 
      stencil[3] = 2; 
      stencil[4] = 3; 
    } else if (k == (K - 2)) {
      stencil[0] = K - 4;
      stencil[1] = K - 3;
      stencil[2] = K - 2; 
      stencil[3] = K - 1;
      stencil[4] = 0; 
    } else if (k == (K - 1)) {
      stencil[0] = K - 3; 
      stencil[1] = K - 2;
      stencil[2] = K - 1;
      stencil[3] = 0;
      stencil[4] = 1; 
    } else if (k == K) {
      stencil[0] = K - 2; 
      stencil[1] = K - 1;
      stencil[2] = 0; 
      stencil[3] = 1; 
      stencil[4] = 2; 
    }
  }

  double t[5];
  for (int i = 0; i < 4; i++) {
    t[i] = theta[stencil[i]];
  }
  t[4] = s; 

  int k_p = 0; 
  int k_n = 0;
  for (int i = 0; i < 5; i++) {
    if (t[i] < M_PI_4) {
      k_p++; 
    }
    if (t[i] > (M_2PI - M_PI_4)) {
      k_n++;
    }
  }
  if ((k_p > 0) && (k_n > 0)) {
    for (int i = 0; i < 5; i++) {
      if (t[i] > (M_2PI - M_PI_4)) {
        t[i] -= M_2PI; 
      }
    }
  }

  for (int i = 0; i < 4; i++) {
    t[i] -= t[4]; 
  }

  double mat[4][4]; 

  for (int i = 0; i < 4; i++) {
    mat[i][0] = 1.0; 
    for (int j = 1; j < 4; j++) {
      mat[i][j] = mat[i][j - 1] * t[i]; 
    }
  }

  double a[4], b[4];

  for (int i = 0; i < 4; i++) {
    b[i] = phi[stencil[i]];
  }

  bool status = solveByQRdecomposition<4>(mat, b, a, 4); 

  w0 = a[0];  w1 = a[1];  w2 = a[2] + a[2]; 
}