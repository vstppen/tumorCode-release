# ifndef CARTESIANGRIDANDCONTROLPOINTS0_H
# define CARTESIANGRIDANDCONTROLPOINTS0_H

class CartesianGridAndControlPoints0 {
    public:
        CartesianGridAndControlPoints0(double x_min_, double x_max_, double y_min_, double y_max_, int I_, int J_, 
                                    double** ctrl_points_, int n_ctrl_, 
                                    int M_, 
                                    const std::string& filename_ctrlPts_, 
                                    const std::string& filename_ctrl_nxny_,
                                    const std::string& filename_bdryPts_, 
                                    bool recordOrNot);
        ~CartesianGridAndControlPoints0();  

        double x_min;
        double x_max;
        double y_min;
        double y_max;
        int I;
        int J;

        int n_ctrl;
        double** ctrl_points;
        int M;

        double dx;
        double dy;
        double* Grid_x;
        double* Grid_y;

        double* M_x;
        double* alpha_x;
        double* beta_x;
        double* M_y;
        double* alpha_y;
        double* beta_y;

        double bdry_delta;  
        double* bdry_theta;  

        double** xy;
        double** dxdy;
        double** nxny;
        double** ctrlPts_nxny;
        double** ddxddy;

        bool** interior;
        int int_node_num;
        int** int_node_index;

        bool** irregular;
        int irreg_node_num;
        int** irreg_node_index;

        bool** near_bdry;
        int near_node_num;
        int** near_node_index;
        int* near_node_proj_to_bdry_index;

        int intersect_num[2];
        int intersect_number;
        int** intersect_index;
        double* intersect_theta;
        double** intersect_coord; 
        double** intersect_tan; 
        double** intersect_nml;

    public:
        void safetyCheck();
        void getPoint1(double t, double &x, double &y);
        void getPoint2(double theta, double &x, double &y, double &tx, double &ty, double &nx, double &ny);
        void getPoint3(double theta, double &x, double &y, double &dx, double &dy, double &ddx, double &ddy);
        
        void generate_boundary_points(const std::string& filename_ctrlPts, const std::string& filename_ctrl_nxny, const std::string& filename_bdryPts, bool recordOrNot); 

        int findClosestIndex(double p, double q);
        void findClosestPointRecursive(double p, double q, double t_start, double t_end, double &t0, double &x0, double &y0); 
        void findClosestPoint(double p, double q, double &t0, double &x0, double &y0);
        double computeParameter(double x, double y);

        bool interiorOrNot(double p, double q, double x, double y, double nx, double ny);
        bool interiorOrNot(double p, double q);
        void identifyInteriorGridNodes();
        void identifyIrregularGridNodes();
        void identifyNearBdryGridNodes();

        double findXIntersectByBisection(double a, double b, double q, double tol);
        double findYIntersectByBisection(double p, double a, double b, double tol); 
        void findIntersectionPoints();

        bool isDuplicate(double** coor, int count, double x, double y, double epsilon);
        void findGoodPoints(int i, int j, double** coor, double* u, double** solution, double (*g)(double, double), bool selfOrNot, int& count);

        void computeDDXandDN(double xi, double eta, double ddx[2], double dn[2]);
        void computeJumps1(double x, double y, double tx, double ty, double nx, double ny, double psi, double d_phi, double du[2]); 
        void computeJumps2(double x, double y, double s, double tx, double ty, double nx, 
                                            double ny, double d_psi, double dd_phi, const double du[2], double ddu[3]);
        void computeJumps(const double *phi_vec, const double *psi_vec, int K, double x0, double y0, double s, double tx, double ty, double nx, double ny, double &j0, double j1[2], double j2[3]);
        void computeJumps2(double x, double y, double s, double tx, double ty, double nx, double ny, double f, double ddu[3]);
        void computeJumps(double x0, double y0, double s, double tx, double ty, double nx, double ny, double f, double jmp2[3]);
        
        void makeCorrection(double (*F)(double, double), double **b);

        void makeCorrection(double (*F)(double, double, double, double, int, int, double**, double, double), double** source, double **b);

        bool computeBoundaryValues(double **u, double xi, double eta, const double jump[6], double bdry_u[12]);

        void extractDirichletBoundaryData(double **w, double (*F)(double, double), double *bdry_w_vec, int K); 
        void extractNeumannBoundaryData(double **w, double (*F)(double, double), double *bdry_wn_vec, int K);

        void extractDirichletBoundaryData(double **w, double (*F)(double, double, double, double, int, int, double**, double, double), double** source, double *bdry_w_vec, int K); 
        void extractNeumannBoundaryData(double **w, double (*F)(double, double, double, double, int, int, double**, double, double), double** source, double *bdry_wn_vec, int K);

};


# endif