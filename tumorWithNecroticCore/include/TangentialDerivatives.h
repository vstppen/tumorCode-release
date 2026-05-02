# ifndef TANGENTIALDERIVATIVES_H
# define TANGENTIALDERIVATIVES_H

void computeFirstDerivative(double delta, const double *theta, const double *psi, int K, 
                            double s, double &w0, double &w1);

void computeFirstAndSecondDerivatives(double delta, const double *theta, const double *phi, 
                            int K, double s, double &w0, double &w1, double &w2);


# endif