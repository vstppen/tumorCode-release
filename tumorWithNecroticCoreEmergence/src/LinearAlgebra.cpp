#include <iostream> 
#include <fstream>
#include <cstdlib>
#include <cmath> 
// #include <mgl2/mgl.h>
#include <limits.h>
#include <ctime>    // For time()

#include <iomanip>    // 设置格式（例如固定小数点格式和精度）
#include <string>     // 字符串支持
#include <cstring>     // 字符串支持

#include "Const.h"
#include "LinearAlgebra.h"

double randDouble(void) 
{
  double v = (rand() % INT_MAX) / (INT_MAX - 1.0); 
  return v;
}

double randDouble(double a, double b)
{
  double v = a + (b - a) * randDouble(); 
  return v;
}

double Random(double a, double b) 
{
  static unsigned int seed = 0; 
  const int NUMBER_OF_ROTATION_BITS = 20;
  const unsigned int MAX_UNSIGNED_RANDOM = 1048575; 
  if (!seed) {
    unsigned int t = time(0); 
    int n = 0; 
    unsigned int temp = t; 
    while (temp) {
      temp /= 10; 
      n++; 
    }
    int m = 1; 
    for (int i = 0; i < (n * 3) / 4; i++) {
      m *= 10;
    }
    unsigned int s = t / m;
    unsigned int r = t - s * m; 
    seed = r * m + s + r;
  }

  srand(seed); 
  unsigned int k = rand()%MAX_UNSIGNED_RANDOM;
  double s = static_cast<double>(k) / (MAX_UNSIGNED_RANDOM);
  seed = rand();

  double r = a + s * (b - a);
  return r; 
}

int* createIntVector(int size) {
  int* vector = (int*)malloc(size * sizeof(int));
  if (vector == NULL) {
    perror("Memory allocation failed");
    exit(1);
  }
  for (int i = 0; i < size; i ++) {
    vector[i] = 0;
  }
  return vector;
}

double* createVector(int size) {
  double* vector = (double*)malloc(size * sizeof(double));
  if (vector == NULL) {
    perror("Memory allocation failed");
    exit(1);
  }
  for (int i = 0; i < size; i ++) {
    vector[i] = 0.0;
  }
  return vector;
}

void freeIntVector(int*& vector) {
  free(vector);
  vector = nullptr;
}

void freeVector(double*& vector) {
  free(vector);
  vector = nullptr;
}


double** createMatrix(int rows, int cols) {
  double** matrix = (double**)malloc(rows * sizeof(double*));
  if (matrix == NULL) {
    perror("Memory allocation failed");
    exit(1);
  }
  for (int i = 0; i < rows; i++) {
    matrix[i] = (double*)malloc(cols * sizeof(double));
    if (matrix[i] == NULL) {
      perror("Memory allocation failed");
      exit(1);
    }
  }
  for (int i = 0; i < rows; i ++) {
    for (int j =0; j < cols; j ++) {
      matrix[i][j] = 0.0;
    }
  }
  return matrix;
}

void freeMatrix(double**& matrix, int rows) {
    if (matrix == nullptr) return;  

    for (int i = 0; i < rows; i++) {
        if (matrix[i] != nullptr) {
            free(matrix[i]);  
            matrix[i] = nullptr;  
        }
    }

    free(matrix);  
    matrix = nullptr;  
}

bool** createBoolMatrix(int rows, int cols) {
  bool** matrix = (bool**)malloc(rows * sizeof(bool*));
  if (matrix == NULL) {
    perror("Memory allocation failed");
    exit(1);
  }
  for (int i = 0; i < rows; i++) {
    matrix[i] = (bool*)malloc(cols * sizeof(bool));
    if (matrix[i] == NULL) {
      perror("Memory allocation failed");
      exit(1);
    }
  }
  for (int i = 0; i < rows; i ++) {
    for (int j =0; j < cols; j ++) {
      matrix[i][j] = false;
    }
  }
  return matrix;
}

void freeBoolMatrix(bool**& matrix, int rows) {
    if (matrix == nullptr) return;  

    for (int i = 0; i < rows; i++) {
        if (matrix[i] != nullptr) {
            free(matrix[i]);  
            matrix[i] = nullptr;  
        }
    }

    free(matrix);  
    matrix = nullptr;  
}

int** createIntMatrix(int rows, int cols) {
  int** matrix = (int**)malloc(rows * sizeof(int*));
  if (matrix == NULL) {
    perror("Memory allocation failed");
    exit(1);
  }
  for (int i = 0; i < rows; i++) {
    matrix[i] = (int*)malloc(cols * sizeof(int));
    if (matrix[i] == NULL) {
      perror("Memory allocation failed");
      exit(1);
    }
  }
  for (int i = 0; i < rows; i ++) {
    for (int j =0; j < cols; j ++) {
      matrix[i][j] = 0;
    }
  }
  return matrix;
}

void freeIntMatrix(int**& matrix, int rows) {
    if (matrix == nullptr) return;  

    for (int i = 0; i < rows; i++) {
        if (matrix[i] != nullptr) {
            free(matrix[i]);  
            matrix[i] = nullptr;  
        }
    }

    free(matrix);  
    matrix = nullptr;  
}

void copy(double* dest, double* src, int n) {
    for (int i = 0; i < n; ++i) {
        dest[i] = src[i];
    }
}

double computeL2Norm_L(const double x[], int n) 
{
  double sum = 0.0; 
  for (int i = 0; i < n; i++) {
    sum += x[i] * x[i];
  }
  return sqrt(sum); 
}

double computeMaxNorm_L(double *b, int n)
{
  double max_b = 0.0;
  for (unsigned int i = 0; i < n; i++) {
    double d = fabs(b[i]); 
    if (d > max_b) {
      max_b = d; 
    }
  }
  return max_b;
}

double computeMaxError_L(double *x, double *y, int n) 
{
  double e = 0.0; 
  for (unsigned int i = 0; i < n; i++) {
    double d = fabs(x[i] - y[i]);
    if (d > e) {
      e = d;
    }
  }
  return e; 
}

double computeNormalizedL2Norm_L(double *b, int n)
{
  double e = 0.0; 
  for (unsigned int i = 0; i < n; i++) {
    e += (b[i]) * (b[i]);
  }
  e = sqrt(e) / sqrt(n);
  return e; 
}

double computeNormalizedL2Error_L(double *x, double *y, int n) 
{
  double e = 0.0; 
  for (unsigned int i = 0; i < n; i++) {
    e += (x[i] - y[i]) * (x[i] - y[i]);
  }
  e = sqrt(e) / sqrt(n);
  return e; 
}

double computeMaxError_L(double **x, double **y, int n, int m) {
    double e = 0.0; 
    for (unsigned int i = 0; i < n; i++) {
      for (unsigned int j = 0; j < m; j++) {
        double d = fabs(x[i][j] - y[i][j]);
        if (d > e) {
          e = d;
        }
      }
    }
    return e; 
}

double computeNormalizedL2Error_L(double **x, double **y, int n, int m) {
  double e = 0.0; 
  for (unsigned int i = 0; i < n; i++) {
    for (unsigned int j = 0; j < m; j++) {
      e += (x[i][j] - y[i][j]) * (x[i][j] - y[i][j]);
    }
  }
  e = sqrt(e) / sqrt(n * m);
  return e; 
}


void printVector(double *x, int n)
{
  for (unsigned int i = 0; i < n; i++) {
    std::cout << x[i] << " ";
  }
  std::cout << std::endl;
}

void printVector(int *x, int n)
{
  for (unsigned int i = 0; i < n; i++) {
    std::cout << x[i] << " ";
  }
  std::cout << std::endl;
}

void printMatrix(double **A, int n, int m) 
{ 
  for (unsigned int i = 0; i < n; i++) {
    for (unsigned int j = 0; j < m; j++) {
      if (fabs(A[i][j]) < EPSILON) {
        // printf("0.000  "); 
        printf("        "); 
      } 
      else {
        printf("%7.3f ", A[i][j]); 
      }
    }
  std::cout << std::endl; 
 }
 std::cout << std::endl; 
}


void saveVectorToFile(const double* vector, int size, const std::string& filename) {
    std::ofstream outFile(filename);

    if (!outFile.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
        return;
    }

    // 设置输出格式，增加数值精度
    outFile << std::fixed << std::setprecision(10);

    for (int i = 0; i < size; ++i) {
        outFile << vector[i];
        if (i < size - 1) {
            outFile << " ";  // Add a space between elements
        }
    }

    outFile << "\n";  // New line at the end
    outFile.close();

    std::cout << "Vector saved to " << filename << std::endl;
}


void saveMatrixToFile(const double** matrix, int rows, int cols, const std::string& filename) {
    std::ofstream outFile(filename);

    if (!outFile.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
        return;
    }

    // 设置输出格式，增加数值精度
    outFile << std::fixed << std::setprecision(10);

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            outFile << matrix[i][j];
            if (j < cols - 1) {
                outFile << " ";  // Add a space between columns
            }
        }
        outFile << "\n";  // New line at the end of each row
    }

    outFile.close();
    std::cout << "Matrix saved to " << filename << std::endl;
}


void randomlyGeneratingData(double* b, int n, double min, double max) 
{
  srand(time(0));

  for (int i = 0; i < n; i ++) {
    b[i] = randDouble(min, max);
  }
}


void randomlyGeneratingData(double** A, int m, int n, double min, double max) 
{
  srand(time(0));

  for (int i = 0; i < m; i ++) {
    for (int j = 0; j < n; j ++) {
      A[i][j] = randDouble(min, max);
    }
  }
}


void MatrixAddingMatrix_L(double** A, double** B, double** result, int m, int n) 
{
  for (int i = 0; i < m; i ++) {
    for (int j = 0; j < n; j ++) {
      result[i][j] = A[i][j] + B[i][j];
    }
  }
}


void MatrixMinusMatrix_L(double** A, double** B, double** result, int m, int n) 
{
  for (int i = 0; i < m; i ++) {
    for (int j = 0; j < n; j ++) {
      result[i][j] = A[i][j] - B[i][j];
    }
  }
}

void MatrixMulVector_L(double** A, double* b, double* result, int m, int n)  // A is m*n while B is n*1
{
  for (int i = 0; i < m; i ++) {
    result[i] = 0.0;
    for (int k = 0; k < n; k ++) {
      result[i] += A[i][k] * b[k];
    }
  }
}

void MatrixMulMatrix_L(double** A, double** B, double** result, int m, int n, int l)  // A is m*n while B is n*l
{
  // #pragma omp parallel for
  for (int i = 0; i < m; i ++) {
    for (int j = 0; j < l; j ++) {
      result[i][j] = 0.0;
      for (int k = 0; k < n; k ++) {
        result[i][j] += A[i][k] * B[k][j];
      }
    }
  }
}

void MatrixMulMatrix_openmp_L(double** A, double** B, double** result, int m, int n, int l)  // A is m*n while B is n*l
{
  #pragma omp parallel for 
  for (int i = 0; i < m; i ++) {
    for (int j = 0; j < l; j ++) {
      result[i][j] = 0.0;
      for (int k = 0; k < n; k ++) {
        result[i][j] += A[i][k] * B[k][j];
      }
    }
  }
}

double InnerProduct_L(double* a, double* b, int n)
{
  double result = 0.0;
  for (int i = 0; i < n; i ++) {
    result += a[i] * b[i];
  }
  return result;
}


double VectorMatrixVector_L(double* b, double** A, double* c, int m, int n)  // b^T * A * c, where A is m*n
{
  double* Ac = createVector(m);
  MatrixMulVector_L(A, c, Ac, m, n);
  double result = InnerProduct_L(b, Ac, m);

  freeVector(Ac);
  return result;
}


int inverse(double** A, double** result, int n) 
{
  double** a = createMatrix(n, n);
  for (int i = 0; i < n; i ++) {
    for (int j = 0; j < n; j ++) {
      a[i][j] = A[i][j];
    }
  }

  for ( int i = 0; i < n; i ++) {
    result[i][i] = 1.0;
  }

  for ( int j = 0; j < n; j ++)  {  // let a[i, j] = 0 for all i != j
    int p = j;  
    double maxV = fabs(a[j][j]);  
    for ( int i = j + 1; i < n; i ++) { // 找到第j列中元素绝对值最大行  
      if ( maxV < fabs(a[i][j]) ) {  
        p = i;  
        maxV = fabs(a[i][j]);  
      }  
    }  
   
  if ( maxV < 1E-12) {  // matrix is not invertible  
    for ( int i = 0; i < n; i ++) {
      for (int k = 0; k < n; k ++) {
        result[i][k] = 0.0;
      }
    }
    return 0;
  }  
 
  if ( j!= p ) {   
    for (int k = 0; k < n; k ++) {  // swap row j and row p
      double temp = a[j][k];
      a[j][k] = a[p][k];
      a[p][k] = temp;
      double temp2 = result[j][k];
      result[j][k] = result[p][k];
      result[p][k] = temp2;
    }
  }  
          
  double d = a[j][j];  
  for (int i = 0; i < n; i ++) { 
    a[j][i] /= d;  // row(j) = row(j) / a[j, j]
    result[j][i] /= d;
  }  
    
  for (int i = 0; i < n; i ++) {  
    if ( i != j ) {  
      double q = a[i][j];  
        for (int k = 0; k < n; k ++) {  
          a[i][k] -= q * a[j][k];  // row(i) = row(i) - a[i, j] * row(j)
          result[i][k] -= q * result[j][k];
        }  
      }  
    }  
  }   

  freeMatrix(a, n);
  return 1;
}


int solveLinearSystemByGaussElimination(double** A, double* b, double* x, int n) 
{
  double** A_copy = createMatrix(n, n);
  double* b_copy = createVector(n);
  for (int i = 0; i < n; i ++) {
    b_copy[i] = b[i];
    for (int j = 0; j < n; j ++) {
      A_copy[i][j] = A[i][j];
    }
  }

  for (unsigned int k = 0; k < (n - 1); k++) {
    for (unsigned int i = k + 1; i < n; i++) {
      if (fabs(A_copy[k][k]) < 1E-10) {
        return 0;
      }
      float lik = A_copy[i][k] / A_copy[k][k];
      for (unsigned int j = k; j < n; j++) {
        A_copy[i][j] = A_copy[i][j] - lik * A_copy[k][j];
      }
      b_copy[i] = b_copy[i] - lik * b_copy[k];
    }
  }

  for (int k = (n - 1); k >= 0; k--) {
    float sum = 0.0;
    for (unsigned int j = k + 1; j < n; j++) {
      sum += A_copy[k][j] * x[j];
    }
    x[k] = (b[k] - sum) / A_copy[k][k]; 
  }

  freeMatrix(A_copy, n);

  return 1;
}

int SolveLinearSystemByRichardson(double** A, double* b, double* x, int n, double omega, bool printOrNot)
{
  double rtol = 1.0E-8; // relative tolerance
  double res_norm0 = computeMaxNorm_L(b, n);

  double* r = createVector(n); // residual
  double atol = res_norm0 * rtol;
  int max_itr_num = 1000000;

  int count = 0;
  double res_norm = 0.0;

  // r = b - Ax
  for (unsigned int i = 0; i < n; i++) {
    double sum = 0.0;
    for (unsigned int j = 0; j < n; j++) {
      sum += A[i][j] * x[j];
    }
    r[i] = b[i] - sum;
  }

  res_norm = computeMaxNorm_L(r, n);

  do {
    count++;

    // x = x + omega * r
    for (unsigned int i = 0; i < n; i++) {
      x[i] += omega * r[i];
    }

    // r = b - A * x
    for (unsigned int i = 0; i < n; i++) {
      double sum = 0.0;
      for (unsigned int j = 0; j < n; j++) {
        sum += A[i][j] * x[j];
      }
      r[i] = b[i] - sum;
    }

    res_norm = computeMaxNorm_L(r, n);

    if (printOrNot == 1 && count % 1000 == 0) {
      std::cout << "After iterating " << count << " times, the residual is " << res_norm 
      << " and res_norm_" << count << " / res_norm_" << count - 1 << "=" << res_norm / res_norm0 << std::endl;
    }

    res_norm0 = res_norm;

  } while ((res_norm > atol) && (count < max_itr_num)) ;

  if (printOrNot == 1) {
    std::cout << "After iterating " << count << " times, the residual is " << res_norm  
    << " and the iteration stops." << std::endl;
  }

  freeVector(r);
  return count;
}

int SolveLinearSystemByJacobi(double** A, double* b, double* x, int n, bool printOrNot) 
{
  double rtol = 1.0E-8; // relative tolerance

  double res_norm0 = computeMaxNorm_L(b, n);

  double* x_old = createVector(n);
  double *r = createVector(n); // residual

  for (unsigned int i = 0; i < n; i++) {
    r[i] = b[i];
  }

  double atol = res_norm0 * rtol; // absolute tolerance
  int max_itr_num = 1000000;
  
  int count = 0;
  double res_norm = 0.0;

  do {
    count++;

    for (unsigned int i = 0; i < n; i++) {
      double sum = 0.0;
      for (unsigned int j = 0; j < i; j++) {
        sum += A[i][j] * x_old[j]; 
      }
      for (unsigned int j = i + 1; j < n; j++) {
        sum += A[i][j] * x_old[j]; 
      }  
      x[i] = (b[i] - sum) / A[i][i];
    }

    for (unsigned int i = 0; i < n; i++) {
      double sum = 0.0; 
      for (unsigned int j = 0; j < n; j++) {
        sum += A[i][j] * x[j];
      }
      r[i] = b[i] - sum;
    }
    
    for (unsigned int i = 0; i < n; i++) {
      x_old[i] = x[i];
    }

    res_norm = computeMaxNorm_L(r, n);

    if (printOrNot == 1 and count % 1000 == 0) {
      std::cout << "After iterating " << count << " times, the residual is " << res_norm 
      << " and res_norm_" << count << " / res_norm_" << count - 1 << "=" << res_norm / res_norm0 << std::endl;
    }

    res_norm0 = res_norm;

  } while ((res_norm > atol) && (count < max_itr_num)) ;

  if (printOrNot == 1) {
    std::cout << "After iterating " << count << " times, the residual is " << res_norm  
    << " and the iteration stops." << std::endl;
  }

  // if (count > max_itr_num -10) {
  //   std::cout << "f" << std::endl;
  // }

  freeVector(r);
  freeVector(x_old);

  return count;
}

int SolveLinearSystemByGS(double** A, double* b, double* x, int n, bool printOrNot)  
{
  double rtol = 1.0E-8; // relative tolerance

  double res_norm0 = computeMaxNorm_L(b, n);

  double *r = new double[n]; // residual

  for (unsigned int i = 0; i < n; i++) {
    r[i] = b[i];
  }

  double atol = res_norm0 * rtol; // absolute tolerance
  int max_itr_num = 1000000;
  
  int count = 0;
  double res_norm = 0.0;

  do {
    count++;

    for (unsigned int i = 0; i < n; i++) {
      double sum = 0.0;
      for (unsigned int j = 0; j < i; j++) {
        sum += A[i][j] * x[j]; // using the current (new) value
      }
      for (unsigned int j = i + 1; j < n; j++) {
        sum += A[i][j] * x[j]; // using the previous (old) value
      }  
      x[i] = (b[i] - sum) / A[i][i];
    }

    for (unsigned int i = 0; i < n; i++) {
      double sum = 0.0; 
      for (unsigned int j = 0; j < n; j++) {
        sum += A[i][j] * x[j];
      }
      r[i] = b[i] - sum;
    }
    
    res_norm = computeMaxNorm_L(r, n);

    if (printOrNot == 1 and count % 1000 == 0) {
      std::cout << "After iterating " << count << " times, the residual is " << res_norm 
      << " and res_norm_" << count << " / res_norm_" << count - 1 << "=" << res_norm / res_norm0 << std::endl;
    }

    res_norm0 = res_norm;

  } while ((res_norm > atol) && (count < max_itr_num)) ;

  if (printOrNot == 1) {
    std::cout << "After iterating " << count << " times, the residual is " << res_norm  
    << " and the iteration stops." << std::endl;
  }

  delete[] r;

  return count;
}


int SolveLinearSystemBySOR(double** A, double* b, double* x, int n, double w, bool printOrNot)  
{
  double rtol = 1.0E-8; // relative tolerance

  double res_norm0 = computeMaxNorm_L(b, n);

  double* x_old = createVector(n);
  double *r = createVector(n); // residual

  for (unsigned int i = 0; i < n; i++) {
    x[i] = 0.0;
    r[i] = b[i];
  }

  double atol = res_norm0 * rtol; // absolute tolerance
  int max_itr_num = 1000000;
  
  int count = 0;
  double res_norm = 0.0;

  do {
    count++;

    for (unsigned int i = 0; i < n; i++) {
      double sum = 0.0;
      for (unsigned int j = 0; j < i; j++) {
        sum += A[i][j] * x[j]; 
      }
      for (unsigned int j = i + 1; j < n; j++) {
        sum += A[i][j] * x[j]; 
      }  
      x[i] = (b[i] - sum) / A[i][i];
      x[i] = (1 - w) * x_old[i] + w * x[i]; 
    }

    for (unsigned int i = 0; i < n; i++) {
      double sum = 0.0; 
      for (unsigned int j = 0; j < n; j++) {
        sum += A[i][j] * x[j];
      }
      r[i] = b[i] - sum;
    }
    
    for (unsigned int i = 0; i < n; i++) {
      x_old[i] = x[i];
    }

    res_norm = computeMaxNorm_L(r, n);

    if (printOrNot == 1 and count % 1000 == 0) {
      std::cout << "After iterating " << count << " times, the residual is " << res_norm 
      << " and res_norm_" << count << " / res_norm_" << count - 1 << "=" << res_norm / res_norm0 << std::endl;
    }

    res_norm0 = res_norm;

  } while ((res_norm > atol) && (count < max_itr_num)) ;

  if (printOrNot == 1) {
    std::cout << "After iterating " << count << " times, the residual is " << res_norm  
    << " and the iteration stops." << std::endl;
  }

  freeVector(r);
  freeVector(x_old);

  return count;
}


int SolveLinearSystemBySD(double **A, double *b, double* x, int n, bool printOrNot) 
{
  double rtol = 1.0E-8; // relative tolerance

  double res_norm0 = computeMaxNorm_L(b, n);

  double atol = res_norm0 * rtol; // absolute tolerance
  int max_itr_num = 100000;

  double* r = createVector(n);
  double* d = createVector(n);

  // initial the residual
  for (int i = 0; i < n; i ++) {
    x[i] = 0.0;  // initial guess of x is zero-vector
    r[i] = b[i];
    // for (int j = 0; j < n; j ++) {
    //   r[i] -= A[i][j] * x[j];
    // }
  }

  double res_norm = InnerProduct_L(r, r, n);
  double alpha = 0.0;
  double* Ad = createVector(n);

  int count = 0;

  do {
    count++;

    for (int i = 0; i < n; i ++) {
      d[i] = r[i];
    }

    alpha = InnerProduct_L(r, r, n) / VectorMatrixVector_L(r, A, r, n, n);

    MatrixMulVector_L(A, d, Ad, n, n);

    for (int i = 0; i < n; i ++) {
      x[i] += alpha * d[i];
      r[i] -= alpha * Ad[i];
    }

    if (printOrNot == 1 and count % 10 == 0) {
      std::cout << "After iterating " << count << " times, the residual is " << res_norm 
      << " and res_norm_" << count << " / res_norm_" << count - 1 << "=" << res_norm / res_norm0 << std::endl;
    }

    res_norm = computeMaxNorm_L(r, n);

  } while ((res_norm > atol) && (count < max_itr_num)) ;
    
  freeVector(r);
  freeVector(d);  
  freeVector(Ad);

  return count;
}


int SolveLinearSystemByCG(double **A, double *b, double* x, int n, bool printOrNot) 
{
  double rtol = 1.0E-10; // relative tolerance

  double res_norm0 = computeMaxNorm_L(b, n);

  double atol = res_norm0 * rtol; // absolute tolerance
  int max_itr_num = 10000;

  double* r_old = createVector(n);
  double* r = createVector(n);
  double* d = createVector(n);

  // initial the residual
  for (int i = 0; i < n; i ++) {
    x[i] = 0.0;  // initial guess of x is zero-vector
    r_old[i] = b[i];
    r[i] = b[i];
  }

  double res_norm = InnerProduct_L(r, r, n);
  double alpha = 0.0;
  double* Ad = createVector(n);

  int count = 0;

  {
    count++;

    for (int i = 0; i < n; i ++) {
      d[i] = r[i];
    }

    alpha = InnerProduct_L(r, r, n) / VectorMatrixVector_L(r, A, r, n, n);

    MatrixMulVector_L(A, d, Ad, n, n);

    for (int i = 0; i < n; i ++) {
      x[i] += alpha * d[i];
      r[i] -= alpha * Ad[i];
    }

    if (printOrNot == 1 and count % 10 == 0) {
      std::cout << "After iterating " << count << " times, the residual is " << res_norm 
      << " and res_norm_" << count << " / res_norm_" << count - 1 << "=" << res_norm / res_norm0 << std::endl;
    }

  }  // run SD once

  double beta = 0.0;

  do {
    count++;

    beta = InnerProduct_L(r, r, n) / InnerProduct_L(r_old, r_old, n);

    for (int i = 0; i < n; i ++) {
      r_old[i] = r[i];
    }

    for (int i = 0; i < n; i ++) {
      d[i] = r[i] + beta * d[i];
    }

    alpha = InnerProduct_L(r, r, n) / VectorMatrixVector_L(d, A, d, n, n);

    MatrixMulVector_L(A, d, Ad, n, n);

    for (int i = 0; i < n; i ++) {
      x[i] += alpha * d[i];
      r[i] -= alpha * Ad[i];
    }

    if (printOrNot == 1 and count % 1 == 0) {
      std::cout << "After iterating " << count << " times, the residual is " << res_norm 
      << " and res_norm_" << count << " / res_norm_" << count - 1 << "=" << res_norm / res_norm0 << std::endl;
    }

    res_norm = computeMaxNorm_L(r, n);

  } while ((res_norm > atol) && (count < max_itr_num)) ;

  if (count > max_itr_num - 10) {
    std::cout << count << std::endl;
  }
    
  freeVector(r_old);
  freeVector(r);
  freeVector(d);  
  freeVector(Ad);

  return count;
}


double makeQRdecompositionInGMRES(double **H, double **R, double *cs, double *sn, double *b, int m, double atol)
{
    int m1 = m + 1;

    int ell = m - 1;
    for (int j = 0; j < m1; j++) {
        R[j][ell] = H[j][ell]; 
    }

    for (int k = 0; k < ell; k++) {
        double t = cs[k] * R[k][ell] + sn[k] * R[k + 1][ell];
        R[k + 1][ell] = - sn[k] * R[k][ell] + cs[k] * R[k + 1][ell]; 
        R[k][ell] = t;
    }

    if (fabs(R[ell + 1][ell]) > atol) {
        double x = R[ell][ell];
        double y = R[ell + 1][ell];
        double r = sqrt(x * x + y * y); 
        cs[ell] = x / r;
        sn[ell] = y / r; 

        int k = ell; 

        double t = cs[k] * R[k][ell] + sn[k] * R[k + 1][ell]; 
        R[k + 1][ell] = - sn[k] * R[k][ell] + cs[k] * R[k + 1][ell]; 
        R[k][ell] = t;

        t = cs[k] * b[k] + sn[k] * b[k + 1];
        b[k + 1] = - sn[k] * b[k] + cs[k] * b[k + 1]; 
        b[k] = t;
    }

    double res = b[m];

    return res; 
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
bool SolveLinearSystemByGMRES(double** A, double *b, double *u, int n, int max_m, int max_itr_num, double tol, int &itr_num)
{
    double atol = 1.0E-15;

    itr_num = 0; 

    double *r = new double[n];
    double *w = new double[n];
    for (int i = 0; i < n; i++) {
        r[i] = w[i] = 0.0;
    }

    MatrixMulVector_L(A, u, w, n, n);
    for (int i = 0; i < n; i++) {
        r[i] = b[i] - w[i]; 
    }   // r = b - Au
    // printVector(w, n);
    double denom = sqrt(n);
    double norm_r0 = computeL2Norm_L(r, n);
    if ((norm_r0 / denom) < atol) {
        delete[] r; r = 0;
        delete[] w; w = 0;
        return true;
    }   // if r is small enough, u is the solution to Au = b

    int max_m1 = max_m + 1;
    int max_m2 = max_m + 2;

    double **V = new double*[max_m2]; 
    double **H = new double*[max_m2];
    double **R = new double*[max_m2]; 

    for (int i = 0; i < max_m2; i++) {
        V[i] = new double[n + 1]; 
        H[i] = new double[max_m1];
        R[i] = new double[max_m1];
    }

    for (int i = 0; i < max_m2; i++) {
        for (int j = 0; j <= n; j++) {
            V[i][j] = 0.0;
        }
        for (int j = 0; j < max_m1; j++) {
            H[i][j] = 0.0; 
            R[i][j] = 0.0; 
        }
    }

    double *cs = new double[max_m1]; 
    double *sn = new double[max_m1]; 
    for (int i = 0; i < max_m1; i++) {
        cs[i] = 1.0;
        sn[i] = 0.0; 
    }

    double *z = new double[max_m1]; 
    for (int i = 0; i < max_m1; i++) {
        z[i] = 0.0; 
    }

    int m = 0; 
    double beta = norm_r0; 
    double r_beta = 1.0 / beta;
    for (int j = 0; j < n; j++) {
        V[m][j] = r[j] * r_beta; 
    }

    double *c = new double [max_m2]; 
    for (int i = 1; i < max_m2; i++) {
        c[i] = 0.0; 
    }
    c[0] = beta; 

    int done = 0; 
    tol = tol * denom;

    // double* temp = createVector(n);
    // double* rr = createVector(n);
    while ((itr_num < max_itr_num) && (!done)) {
        // std::cout << "GMRES iteration with itr_num = " << itr_num ;

        MatrixMulVector_L(A, V[m], w, n, n);
        for (int j = 0; j <= m; j++) {
            H[j][m] = InnerProduct_L(V[j], w, n);
            for (int i = 0; i < n; i++) {
                w[i] -= V[j][i] * H[j][m];
            }
        }
        H[m + 1][m] = computeL2Norm_L(w, n);

        if (fabs(H[m + 1][m]) > atol) {
            for (int i = 0; i < n; i++) {
                V[m + 1][i] = w[i] / H[m + 1][m];
            }
        } else {
            done = 1; 
        }

        double rz = makeQRdecompositionInGMRES(H, R, cs, sn, c,  m + 1, atol); 

        itr_num++; 

        if (fabs(rz) < tol) {

            for (int i = m; i >= 0; i--) {
                double s = 0.0;
                for (int j = m; j > i; j--) {
                    s += R[i][j] * z[j]; 
                }
                z[i] = (c[i] - s) / R[i][i]; 
            }

            for (int i = 0; i < n; i++) {
                double sum = 0.0; 
                for (int j = 0; j <= m; j++) {
                    sum += V[j][i] * z[j];
                }
                u[i] += sum;
            }

            done = 1; 

        } else {

            if (m == max_m) {

                for (int i = m; i >= 0; i--) {
                    double s = 0.0;
                    for (int j = m; j > i; j--) {
                        s += R[i][j] * z[j]; 
                    }
                    z[i] = (c[i] - s) / R[i][i];
                }

                for (int i = 0; i < n; i++) {
                    double sum = 0.0; 
                    for (int j = 0; j <= m; j++) {
                        sum += V[j][i] * z[j];
                    }
                    u[i] += sum; 
                }
 
                MatrixMulVector_L(A, u, w, n, n);
                for (int i = 0; i < n; i++) {
                    r[i] = b[i] - w[i]; 
                }
                norm_r0 = computeL2Norm_L(r, n);
                if ((norm_r0 / sqrt(n)) < atol) {
                    done = 1; 
                } else {
                    m = 0; 
                    beta = norm_r0;
                    double r_beta = 1.0 / beta;
                    for (int j = 0; j < n; j++) {
                        V[m][j] = r[j] * r_beta; 
                    }
                    for (int i = 1; i < max_m2; i++) {
                        c[i] = 0.0; 
                    }
                    c[0] = beta; 
                }
            } else {
                m++; 
            }
        }

        // computeMatrixVectorProduct(G, u, temp);   // w = Au
        // for (int i = 0; i < n; i++) {
        //     rr[i] = b[i] - temp[i]; 
        // }   // r = b - Au
        // std::cout << ", max-norm residual = " << computeMaxNorm(rr, n) << std::endl;
    }
    // freeVector(temp);
    // freeVector(rr);

    // computeMatrixVectorProduct(G, u, w);   // w = Au
    // for (int i = 0; i < n; i++) {
    //     r[i] = b[i] - w[i]; 
    // }   // r = b - Au
    // std::cout << "max-norm residual = " << computeMaxNorm(r, n) << std::endl;


    if (itr_num >= max_itr_num) {
        std::cout << "GMRES iteration failed." << std::endl;
        exit(EXIT_FAILURE);
        return false; 
    }

    delete[] cs; cs = 0; 
    delete[] sn; sn = 0; 

    for (int i = 0; i < max_m2; i++) {
        delete[] V[i]; V[i] = 0;
        delete[] H[i]; H[i] = 0; 
        delete[] R[i]; R[i] = 0;
    }
    delete[] H; H = 0;
    delete[] R; R = 0;
    delete[] V; V = 0;

    delete[] r; r = 0;
    delete[] w; w = 0;

    delete[] z; z = 0; 
    delete[] c; c = 0; 

    return true;
}


bool SolveLinearSystemByBiCGStab(double** A, double* b, double* x, int n) {
    const double TOL = 1e-8;
    const int MAX_ITER = 1000;

    double *r = new double[n];
    double *r_hat = new double[n];
    double *p = new double[n];
    double *v = new double[n];
    double *s = new double[n];
    double *t = new double[n];

    // r = b - A*x
    MatrixMulVector_L(A, x, r, n, n);
    for (int i = 0; i < n; ++i)
        r[i] = b[i] - r[i];

    copy(r_hat, r, n);
    copy(p, r, n);
    std::fill(v, v + n, 0.0);

    double rho_old = 1, alpha = 1, omega = 1;
    double rho_new, beta;

    for (int iter = 0; iter < MAX_ITER; ++iter) {
        rho_new = InnerProduct_L(r_hat, r, n);
        if (std::fabs(rho_new) < 1e-14) {
            std::cerr << "Breakdown: rho too small\n";
            break;
        }

        if (iter > 0) {
            beta = (rho_new / rho_old) * (alpha / omega);
            for (int i = 0; i < n; ++i)
                p[i] = r[i] + beta * (p[i] - omega * v[i]);
        }

        MatrixMulVector_L(A, p, v, n, n);
        double v_dot_rhat = InnerProduct_L(v, r_hat, n);
        if (std::fabs(v_dot_rhat) < 1e-14) {
            std::cerr << "Breakdown: v·r_hat too small\n";
            break;
        }
        alpha = rho_new / v_dot_rhat;

        for (int i = 0; i < n; ++i)
            s[i] = r[i] - alpha * v[i];

        double s_norm = std::sqrt(InnerProduct_L(s, s, n));
        if (s_norm < TOL) {
            for (int i = 0; i < n; ++i)
                x[i] += alpha * p[i];
            delete[] r; delete[] r_hat; delete[] p;
            delete[] v; delete[] s; delete[] t;
            return true;
        }

        MatrixMulVector_L(A, s, t, n, n);
        double t_dot_s = InnerProduct_L(t, s, n);
        double t_dot_t = InnerProduct_L(t, t, n);
        if (std::fabs(t_dot_t) < 1e-14) {
            std::cerr << "Breakdown: t·t too small\n";
            break;
        }
        omega = t_dot_s / t_dot_t;

        for (int i = 0; i < n; ++i)
            x[i] += alpha * p[i] + omega * s[i];

        for (int i = 0; i < n; ++i)
            r[i] = s[i] - omega * t[i];

        double r_norm = std::sqrt(InnerProduct_L(r, r, n));
        if (r_norm < TOL) {
            delete[] r; delete[] r_hat; delete[] p;
            delete[] v; delete[] s; delete[] t;
            return true;
        }

        if (std::fabs(omega) < 1e-14) {
            std::cerr << "Breakdown: omega too small\n";
            break;
        }

        rho_old = rho_new;
    }

    std::cerr << "SolveLinearSystemByBiCGStab did not converge.\n";

    double* Ax = new double[n];
    MatrixMulVector_L(A, x, Ax, n, n);
    double res_norm = 0;
    for (int i = 0; i < n; ++i)
        res_norm += (Ax[i] - b[i]) * (Ax[i] - b[i]);
    res_norm = std::sqrt(res_norm);
    std::cout << "Final residual ||Ax - b|| = " << res_norm << "\n";
    delete[] Ax;

    delete[] r; delete[] r_hat; delete[] p;
    delete[] v; delete[] s; delete[] t;
    return false;
}


double Atan(double x, double y)  // modification for atan(y / x)
{
  double theta = 0.0;
  if (fabs(x) < EPSILON12) {
    if (y > EPSILON12) {
      theta = M_PI_2; 
    } else if (y + EPSILON12 < 0) {
      theta = M_PI + M_PI_2; 
    } else {
      theta = 0; 
    }
  } else {
    theta = atan(y / x);
    if (x < 0) {
      theta += M_PI;
    } else {
      if (y + EPSILON12 < 0) {
        theta += M_2PI;
      }
    }
  }
  return(theta); 
}