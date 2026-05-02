#include <cstdlib>
#include <cstdio>
#include <cmath> 

#include "LinearAlgebra.h"
#include "BirkoffInterpolation.h"
#include "LU.h"

bool computeBoundaryValues(const double p[2], const double jump[6], 
                           const double coord[6][2], const bool interior[6],
                           const double grid_u[6], double bdry_u[12]) 
{
  // jump[0]: u
  // jump[1]: du_dx
  // jump[2]: du_dy
  // jump[3]: d2u_dx2
  // jump[4]: d2u_dy2
  // jump[5]: d2u_dxdy

  int i, j, k; 

  double dx, dy;

  double mat[6][6]; 

  double u[6], b[6];

  double c[6], sum; 

  c[0] = 1.0; 

  for (i = 0; i < 6; i++) {

    b[i] = grid_u[i];

    dx = coord[i][0] - p[0];
    dy = coord[i][1] - p[1]; 

    c[1] = dx; 
    c[2] = dy; 

    c[3] = 0.5 * dx * dx; 
    c[4] = 0.5 * dy * dy; 

    c[5] = dx * dy;

    if (interior[i] == false) {
      sum = jump[0]; 
      for (j = 1; j < 6; j++) {
        sum += c[j] * jump[j]; 
      }
      b[i] += sum;
    }

    for (j = 0; j < 6; j++) {
      mat[i][j] = c[j];
    }
  }

  bool status = LUsolve(mat, b, u); 

  if (status == true) {
    for (i = 0, k = 6; i < 6; i++, k++) {
      bdry_u[k] = u[i] - jump[i];
      bdry_u[i] = u[i]; 
    }
  }

  return status;
}

/**
 *****************************************************************************
 *****************************************************************************
 **/
bool computeBoundaryValues(const double p[2],
                           const double coord[6][2], 
                           const double grid_u[6], 
                           double bdry_u[6]) 
{
  int i, j, k; 

  double dx, dy;

  double mat[6][6]; 

  double u[6], b[6]; 

  double c[6], sum;

  c[0] = 1.0;

  for (i = 0; i < 6; i++) {

    b[i] = grid_u[i];

    dx = coord[i][0] - p[0]; 
    dy = coord[i][1] - p[1]; 

    c[1] = dx; 
    c[2] = dy;

    c[3] = 0.5 * dx * dx;
    c[4] = 0.5 * dy * dy; 

    c[5] = dx * dy; 

    for (j = 0; j < 6; j++) {
      mat[i][j] = c[j];
    }
  }

  bool status = LUsolve(mat, b, u); 

  if (status == true) {
    for (i = 0; i < 6; i++, k++) {
      bdry_u[i] = u[i]; 
    }
  }

  return status; 
}
