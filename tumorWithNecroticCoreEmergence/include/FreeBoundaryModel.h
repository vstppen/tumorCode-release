#ifndef FREEBOUNDARYMODEL_H
#define FREEBOUNDARYMODEL_H

class CartesianGridAndControlPoints0;
class CartesianGridAndControlPoints;

void testModel(double x_min_, double x_max_, double y_min_, double y_max_, int I_, int J_, double** ctrl_points_1, int n_ctrl_1, int M_, double std_time_step, double t_end, double d);

#endif