# ifndef FASTALGORITHM_H
# define FASTALGORITHM_H

inline int Log2(int m);
inline int Power2(int m);
inline bool checkIfPowerOf2(int n);

inline void computeFastFourierTransform(const double f_r[], const double z_r[], const double z_i[], double c_r[], double c_i[], int n);
inline void computeFastFourierTransform(const double f[], double c_r[], double c_i[], int n);
inline void computeFastSineTransform(const double u[], double c[], int m);
inline void computeFastCosineTransform(const double u[], double c[], int m);
inline void computeFastInverseCosineTransform(const double u[], double c[], int m);

void FastPoissonSolver(double **unknown, const int sz[2]);
void FastPoissonSolver(const double sigma[2], double **unknown, const int sz[2]);
void FastMHelmholtzSolver(const double sigma[2], double **unknown, const int sz[2], double eta);
void FastMHelmholtzSolver(double **unknown, const int sz[2], double lambda); 
void FastMHelmholtzSolver2(double **unknown, const int sz[2], double eta);   // ！！！

# endif