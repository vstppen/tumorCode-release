#ifndef __LinearAlgebra_h_IS_INCLUDED__
#define __LinearAlgebra_h_IS_INCLUDED__ 

#include <iostream> 
#include <cstdlib>
#include <cmath> 
#include <limits.h>
#include <ctime>    

double randDouble(void);
double randDouble(double a, double b);
double Random(double a, double b);

int* createIntVector(int size);
double* createVector(int size);
void freeIntVector(int*& vector);
void freeVector(double*& vector);
double** createMatrix(int rows, int cols);
void freeMatrix(double**& matrix, int rows);
bool** createBoolMatrix(int rows, int cols);
void freeBoolMatrix(bool**& matrix, int rows);
int** createIntMatrix(int rows, int cols);
void freeIntMatrix(int**& matrix, int rows);

double computeL2Norm_L(const double x[], int n); 
double computeMaxNorm_L(double *b, int n);
double computeMaxError_L(double *x, double *y, int n) ;
double computeNormalizedL2Norm_L(double *b, int n);
double computeNormalizedL2Error_L(double *x, double *y, int n) ;
double computeMaxError_L(double **x, double **y, int n, int m) ;
double computeNormalizedL2Error_L(double **x, double **y, int n, int m) ;

void printVector(int *x, int n);
void printVector(double *x, int n);
void printMatrix(double **A, int n, int m);
void saveVectorToFile(const double* vector, int size, const std::string& filename);
void saveMatrixToFile(const double** matrix, int rows, int cols, const std::string& filename);

void randomlyGeneratingData(double* b, int n, double min, double max) ;
void randomlyGeneratingData(double** A, int m, int n, double min, double max) ;

void MatrixAddingMatrix_L(double** A, double** B, double** result, int m, int n) ;
void MatrixMinusMatrix_L(double** A, double** B, double** result, int m, int n) ;
void MatrixMulVector_L(double** A, double* b, double* result, int m, int n);  // A is m*n while B is n*1
void MatrixMulMatrix_L(double** A, double** B, double** result, int m, int n, int l);  // A is m*n while B is n*l
void MatrixMulMatrix_openmp_L(double** A, double** B, double** result, int m, int n, int l);  // A is m*n while B is n*l
double InnerProduct_L(double* a, double* b, int n);
double VectorMatrixVector_L(double* b, double** A, double* c, int m, int n);  // b^T * A * c, where A is m*n

int inverse(double** A, double** result, int n); 

int solveLinearSystemByGaussElimination(double** A, double* b, double* x, int n);

int SolveLinearSystemByRichardson(double** A, double* b, double* x, int n, double omega, bool printOrNot);
int SolveLinearSystemByJacobi(double** A, double* b, double* x, int n, bool printOrNot); 
int SolveLinearSystemByGS(double** A, double* b, double* x, int n, bool printOrNot);  
int SolveLinearSystemBySOR(double** A, double* b, double* x, int n, double w, bool printOrNot);  

int SolveLinearSystemBySD(double **A, double *b, double* x, int n, bool printOrNot); 
int SolveLinearSystemByCG(double **A, double *b, double* x, int n, bool printOrNot); 

double makeQRdecompositionInGMRES(double **H, double **R, double *cs, double *sn, double *b, int m, double atol);
bool SolveLinearSystemByGMRES(double** A, double *b, double *u, int n, int max_m, int max_itr_num, double tol, int &itr_num);

bool SolveLinearSystemByBiCGStab(double** A, double* b, double* x, int n);

double Atan(double x, double y);  

#endif