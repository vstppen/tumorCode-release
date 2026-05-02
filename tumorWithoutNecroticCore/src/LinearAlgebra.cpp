#include <iostream> 
#include <fstream>
#include <cstdlib>
#include <cmath> 
// #include <mgl2/mgl.h>
#include <limits.h>
#include <ctime>    // For time()

#include <iomanip>    // 设置格式（例如固定小数点格式和精度）
#include <string>     // 字符串支持

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

double computeMaxNorm(double *b, int n)
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

double computeMaxError(double *x, double *y, int n) 
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

double computeL2Norm(double *b, int n)
{
  double e = 0.0; 
  for (unsigned int i = 0; i < n; i++) {
    e += (b[i]) * (b[i]);
  }
  e = sqrt(e) / sqrt(n);
  return e; 
}

double computeL2Error(double *x, double *y, int n) 
{
  double e = 0.0; 
  for (unsigned int i = 0; i < n; i++) {
    e += (x[i] - y[i]) * (x[i] - y[i]);
  }
  e = sqrt(e) / sqrt(n);
  return e; 
}

double computeMaxError(double **x, double **y, int n, int m) {
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

double computeL2Error(double **x, double **y, int n, int m) {
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


void MatrixAddingMatrix(double** A, double** B, double** result, int m, int n) 
{
  for (int i = 0; i < m; i ++) {
    for (int j = 0; j < n; j ++) {
      result[i][j] = A[i][j] + B[i][j];
    }
  }
}


void MatrixMinusMatrix(double** A, double** B, double** result, int m, int n) 
{
  for (int i = 0; i < m; i ++) {
    for (int j = 0; j < n; j ++) {
      result[i][j] = A[i][j] - B[i][j];
    }
  }
}


void MatrixMulVector(double** A, double* b, double* result, int m, int n)  // A is m*n while B is n*1
{
  for (int i = 0; i < m; i ++) {
    result[i] = 0.0;
    for (int k = 0; k < n; k ++) {
      result[i] += A[i][k] * b[k];
    }
  }
}


void MatrixMulMatrix(double** A, double** B, double** result, int m, int n, int l)  // A is m*n while B is n*l
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

void MatrixMulMatrix_openmp(double** A, double** B, double** result, int m, int n, int l)  // A is m*n while B is n*l
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

double InnerProduct(double* a, double* b, int n)
{
  double result = 0.0;
  for (int i = 0; i < n; i ++) {
    result += a[i] * b[i];
  }
  return result;
}


double VectorMatrixVector(double* b, double** A, double* c, int m, int n)  // b^T * A * c, where A is m*n
{
  double* Ac = createVector(m);
  MatrixMulVector(A, c, Ac, m, n);
  double result = InnerProduct(b, Ac, m);

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


int SolveLinearSystemByJacobi(double** A, double* b, double* x, int n, bool printOrNot) 
{
  double rtol = 1.0E-8; // relative tolerance

  double res_norm0 = computeMaxNorm(b, n);

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

    res_norm = computeMaxNorm(r, n);

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

  double res_norm0 = computeMaxNorm(b, n);

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
    
    res_norm = computeMaxNorm(r, n);

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

  double res_norm0 = computeMaxNorm(b, n);

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

    res_norm = computeMaxNorm(r, n);

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

  double res_norm0 = computeMaxNorm(b, n);

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

  double res_norm = InnerProduct(r, r, n);
  double alpha = 0.0;
  double* Ad = createVector(n);

  int count = 0;

  do {
    count++;

    for (int i = 0; i < n; i ++) {
      d[i] = r[i];
    }

    alpha = InnerProduct(r, r, n) / VectorMatrixVector(r, A, r, n, n);

    MatrixMulVector(A, d, Ad, n, n);

    for (int i = 0; i < n; i ++) {
      x[i] += alpha * d[i];
      r[i] -= alpha * Ad[i];
    }

    if (printOrNot == 1 and count % 10 == 0) {
      std::cout << "After iterating " << count << " times, the residual is " << res_norm 
      << " and res_norm_" << count << " / res_norm_" << count - 1 << "=" << res_norm / res_norm0 << std::endl;
    }

    res_norm = computeMaxNorm(r, n);

  } while ((res_norm > atol) && (count < max_itr_num)) ;
    
  freeVector(r);
  freeVector(d);  
  freeVector(Ad);

  return count;
}


int SolveLinearSystemByCG(double **A, double *b, double* x, int n, bool printOrNot) 
{
  double rtol = 1.0E-10; // relative tolerance

  double res_norm0 = computeMaxNorm(b, n);

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

  double res_norm = InnerProduct(r, r, n);
  double alpha = 0.0;
  double* Ad = createVector(n);

  int count = 0;

  {
    count++;

    for (int i = 0; i < n; i ++) {
      d[i] = r[i];
    }

    alpha = InnerProduct(r, r, n) / VectorMatrixVector(r, A, r, n, n);

    MatrixMulVector(A, d, Ad, n, n);

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

    beta = InnerProduct(r, r, n) / InnerProduct(r_old, r_old, n);

    for (int i = 0; i < n; i ++) {
      r_old[i] = r[i];
    }

    for (int i = 0; i < n; i ++) {
      d[i] = r[i] + beta * d[i];
    }

    alpha = InnerProduct(r, r, n) / VectorMatrixVector(d, A, d, n, n);

    MatrixMulVector(A, d, Ad, n, n);

    for (int i = 0; i < n; i ++) {
      x[i] += alpha * d[i];
      r[i] -= alpha * Ad[i];
    }

    if (printOrNot == 1 and count % 1 == 0) {
      std::cout << "After iterating " << count << " times, the residual is " << res_norm 
      << " and res_norm_" << count << " / res_norm_" << count - 1 << "=" << res_norm / res_norm0 << std::endl;
    }

    res_norm = computeMaxNorm(r, n);

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


// // y = (|x|, 0, ..., 0) or (-|x|, 0, ..., 0) and v = (x - y) / |x - y|
// bool computeHouseholderVector(double *x, int m, double *y, double *v)  
// {
//   for (int i = 0; i < m; i++) {
//     y[i] = 0.0;
//   }

//   double norm_x = computeL2Norm(x, m); 
//   if (norm_x / sqrt(m) < EPSILON12) {
//     for (int i = 0; i < m; i++) {
//       y[i] = x[i]; 
//       v[i] = 0.0; 
//     }
//     return false; 
//   }

//   if (x[0] > EPSILON12) {
//     y[0] = - norm_x;
//   } else {
//     y[0] = norm_x;
//   }

//   for (int i = 0; i < m; i++) {
//     v[i] = x[i];
//   }
//   v[0] -= y[0];

//   double norm_v = sqrt(norm_x * norm_x - x[0] * x[0] + v[0] * v[0]); 
//   if (norm_v < EPSILON12) {
//     for (int i = 0; i < m; i++) {
//       v[i] = 0.0;
//     }
//   } else {
//     double r_norm_v = 1.0 / norm_v;
//     for (int i = 0; i < m; i++) {
//       v[i] *= r_norm_v;
//     }
//   }

//   return true; 
// }


// // b = (|a|, 0, ..., 0) or (-|a|, 0, ..., 0) and w = (a - b) / |a - b| 
// bool computeHouseholderVector(double* a, double* w, int n) 
// {
//   double len = computeL2Norm(a, n); 
//   if (len / sqrt(n) < EPSILON12) {
//     for (int i = 0; i < n; i++) {
//       w[i] = 0.0; 
//     }
//     return false;
//   }

//   if (a[0] > 0.0) {
//     w[0] = a[0] + len;
//   } else {
//     w[0] = a[0] - len; 
//   }
//   for (int i = 1; i < n; i++) {
//     w[i] = a[i]; 
//   }

//   double len2 = computeL2Norm(w, n);
//   if (len2 > EPSILON) {
//     for (int i = 0; i < n; i++) {
//       w[i] /= len2; 
//     }
//   }

//   return true; 
// }


// /**
//  *****************************************************************************
//  * w: the Householder vector                                                 *
//  * n: the dimension of the vectors                                           *
//  * u: the vector before reflection                                           *
//  * v: the vector after reflection                                            *
//  * v = (I - 2 w w^T) * u                                                     *
//  *****************************************************************************
//  **/
// void makeHouseholderTransform(double* w, double* u, int n, double* v)
// {
//   double p = InnerProduct(w, u, n); 
//   for (int i = 0; i < n; i++) {
//     double c = p * w[i]; 
//     v[i] = u[i] - (c + c);
//   }
// }

// /**
//  *****************************************************************************
//  * w: the Householder vector                                                 *
//  * n: the dimension of the vectors                                           *
//  * u: the vector before reflection                                           *
//  * v: the vector after reflection                                            *
//  * v = (I - 2 w w^T) * u                                                     *
//  *****************************************************************************
//  **/
// void makeHouseholderTransform(double* w, int n, double* u, double* v) 
// {
//   double p = InnerProduct(w, u, n);
//   for (int i = 0; i < n; i++) {
//     double c = p * w[i]; 
//     v[i] = u[i] - (c + c); 
//   }
// }

// bool SolveLinearSystemByQRdecomposition(double** A, double* rhs, double* u, int N) 
// {
//   int n = N;
//   double* b = createVector(N); 
//   for (int i = 0; i < n; i++) {
//     b[i] = rhs[i]; 
//   }

//   double* v = createVector(N); 
//   double* x = createVector(N); 
//   double* y = createVector(N); 

//   double** R = createMatrix(N, N);  

//   for (int i = 0; i < n; i++) {
//     for (int j = 0; j < n; j++) {
//       R[i][j] = A[i][j];
//     }
//   }

//   for (int k = 0, m = n; k < n - 1; k++, m--) {
//     for (int i = 0; i < m; i++) {
//       x[i] = R[i + k][k]; 
//     }
//     bool status = computeHouseholderVector(x, m, y, v); 
//     if (status) {
//       for (int i = 0; i < m; i++) {
//         R[i + k][k] = y[i]; 
//       }

//       for (int j = k + 1; j < n; j++) {
//         for (int i = 0; i < m; i++) {
//           x[i] = R[i + k][j];
//         }
//         makeHouseholderTransform(v, m, x, y); 
//         for (int i = 0; i < m; i++) {
//           R[i + k][j] = y[i];
//         }
//       }

//       makeHouseholderTransform(v, m, b + k, y); 
//       for (int i = 0; i < m; i++) {
//         b[k + i] = y[i];
//       }

//     } else {
//       return false; 
//     }
//   }

//   bool status = true;
//   if (fabs(R[n - 1][n - 1]) > EPSILON12) {
//     for (int i = n - 1; i >= 0; i--) {
//       double sum = 0.0; 
//       for (int j = i + 1; j < n; j++) {
//         sum += R[i][j] * u[j]; 
//       }
//       u[i] = (b[i] - sum) / R[i][i]; 
//     }
//   } else {
//     status = false;
//   }

//   freeMatrix(R, N);
//   freeVector(y);
//   freeVector(x);
//   freeVector(v);
//   freeVector(b);
//   return status;
// }