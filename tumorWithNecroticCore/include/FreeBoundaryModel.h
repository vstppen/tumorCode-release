#ifndef FREEBOUNDARYMODEL_H
#define FREEBOUNDARYMODEL_H

class CartesianGridAndControlPoints;

extern double model_G0;

void testModel1(double x_min_, double x_max_, double y_min_, double y_max_, int I_, int J_, double** ctrl_points_0, int n_ctrl_0, double** ctrl_points_1, int n_ctrl_1, int M_, double time_step, int steps);
void testModel3(double x_min_, double x_max_, double y_min_, double y_max_, int I_, int J_, double** ctrl_points_0, int n_ctrl_0, double** ctrl_points_1, int n_ctrl_1, int M_, double time_step, int steps);

#endif