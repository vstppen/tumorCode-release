# ifndef CARTESIANGRIDANDCONTROLPOINTS_H
# define CARTESIANGRIDANDCONTROLPOINTS_H

// !!!0
class CartesianGridAndControlPoints {
    public:
        CartesianGridAndControlPoints(double x_min_, double x_max_, double y_min_, double y_max_, int I_, int J_, 
                                    double** ctrl_points_, int n_ctrl_, 
                                    int M_, 
                                    const std::string& filename_ctrlPts_, 
                                    const std::string& filename_ctrl_nxny_,
                                    const std::string& filename_bdryPts_,
                                    bool recordOrNot_);
        ~CartesianGridAndControlPoints();  

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
        int** int_node_ij_to_l;

        bool** irregular;
        int irreg_node_num;
        int** irreg_node_index;

        int intersect_num[2];
        int intersect_number;
        int** intersect_index;
        double* intersect_theta;
        double** intersect_coord; 
        double** intersect_tan; 
        double** intersect_nml;

    public:
        void getPoint1(double t, double &x, double &y);
        void getPoint2(double theta, double &x, double &y, double &tx, double &ty, double &nx, double &ny);
        void getPoint3(double theta, double &x, double &y, double &dx, double &dy, double &ddx, double &ddy);

        void generate_boundary_points(const std::string& filename_ctrlPts, const std::string& filename_ctrl_nxny, const std::string& filename_bdryPts, bool recordOrNot_); 

        int findClosestIndex(double p, double q);
        void findClosestPointRecursive(double p, double q, double t_start, double t_end, double &t0, double &x0, double &y0); 
        void findClosestPoint(double p, double q, double &t0, double &x0, double &y0);
        double computeParameter(double x, double y);

        bool interiorOrNot(double p, double q, double x, double y, double nx, double ny);
        bool interiorOrNot(double p, double q);
        void identifyInteriorGridNodes();
        void identifyIrregularGridNodes();

        double findXIntersectByBisection(double a, double b, double q, double tol);
        double findYIntersectByBisection(double p, double a, double b, double tol); 
        void findIntersectionPoints();

        void computeDDXandDN(double xi, double eta, double ddx[2], double dn[2]);
        
        void computeJumps1(double x, double y, double tx, double ty, double nx, double ny, double psi, double d_phi, double du[2]); 
        
        void computeJumps2(double x, double y, double s, double tx, double ty, double nx, double ny, double d_psi, double dd_phi, const double du[2], double ddu[3], double lambda, double phi);
        void computeJumps(const double *phi_vec, const double *psi_vec, int K, double x0, double y0, double s, double tx, double ty, double nx, double ny, double &j0, double j1[2], double j2[3], double lambda);
        
        void computeJumps2(double x, double y, double s, double tx, double ty, double nx, double ny, double f, double ddu[3], double lambda, bool interiorOrNot);
        void computeJumps(double x0, double y0, double s, double tx, double ty, double nx, double ny, double f, double jmp2[3], double lambda, bool interiorOrNot);

        void makeCorrection(const double *phi, const double *psi, int K, double **b, double lambda);
        void makeCorrection(double (*F)(double, double), double **b, double lambda, bool interiorOrNot);
        
        bool computeBoundaryValues(double **u, double xi, double eta, const double jump[6], double bdry_u[12]);
        
        void extractDirichletBoundaryData(double **u, const double *phi_vec, const double *psi_vec, double *bdry_u_vec, int K, double lambda);
        void extractDirichletBoundaryData(double **w, double (*F)(double, double), double *bdry_w_vec, int K, double lambda, bool interiorOrNot); 
        void extractNeumannBoundaryData(double **u, const double *phi_vec, const double *psi_vec, double *bdry_un_vec, int K, double lambda);
        void extractNeumannBoundaryData(double **w, double (*F)(double, double), double *bdry_wn_vec, int K, double lambda, bool interiorOrNot);

};


# endif