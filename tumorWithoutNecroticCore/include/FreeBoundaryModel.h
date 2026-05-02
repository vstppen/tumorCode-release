#ifndef FREEBOUNDARYMODEL
#define FREEBOUNDARYMODEL

void evolveCtrlPointsStrategy1(double** ctrl_points, int& n_ctrl, double** nxny, double* v, double time_step);
void evolveCtrlPointsStrategy2(double** ctrl_points, int& n_ctrl, double** nxny, double* v, double time_step);
void testModel(double x_min_, double x_max_, double y_min_, double y_max_, int I_, int J_, double** ctrl_points_, int n_ctrl_, int M_, double time_step, int steps);

#endif