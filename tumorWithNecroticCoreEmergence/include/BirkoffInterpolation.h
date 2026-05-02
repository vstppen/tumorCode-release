# ifndef BIRKOFFINTERPOLATION_H
# define BIRKOFFINTERPOLATION_H

bool computeBoundaryValues(const double p[2], const double jump[6], const double coord[6][2], const bool interior[6], const double grid_u[6], double bdry_u[12]); 
bool computeBoundaryValues(const double p[2], const double coord[6][2], const double grid_u[6], double bdry_u[6]); 

# endif