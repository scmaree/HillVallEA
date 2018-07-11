#pragma once

/*

HillVallEA

By S.C. Maree
s.c.maree[at]amc.uva.nl
github.com/SCMaree/HillVallEA

*/

#include "hillvallea_internal.hpp"
#include "param.hpp"

namespace hillvallea
{
  
  class population_t;

  // sample parameters using normal distribution
  // returns the number of trials before an in-range sample is found.
  // on too many fails, sample uniform.
  //----------------------------------------------
  int sample_normal(vec_t & sample, const size_t problem_size, const vec_t & mean, const matrix_t & chol, const vec_t & lower_param_range, const vec_t & upper_param_range, rng_pt rng);
  int sample_normal_univariate(vec_t & sample, const size_t problem_size, const vec_t & mean, const matrix_t & chol, const vec_t & lower_param_range, const vec_t & upper_param_range, rng_pt rng);

  void sample_uniform(vec_t & sample, const size_t problem_size, const vec_t & lower_user_range, const vec_t & upper_user_range,  rng_pt rng);
  
  // check if the parameter is in range
  //-----------------------------------------
  bool in_range(const vec_t & sample, const vec_t & lower_param_range, const vec_t & upper_param_range);
  bool boundary_repair(vec_t & sample, const vec_t & lower_param_range, const vec_t & upper_param_range);

  // Cholesky decomposition
  // based on BLAS / LINPACK library functions
  //-----------------------------------------
  void choleskyDecomposition(const matrix_t & cov, matrix_t & chol);
  void choleskyDecomposition_univariate(const matrix_t & cov, matrix_t & chol);
  void *Malloc(long size);
  double **matrixNew(int n, int m);
  double vectorDotProduct(const double *vector0, const double *vector1, int n0);
  double *matrixVectorMultiplication(const double **matrix, const double *vector, int n0, int n1);
  double **matrixMatrixMultiplication(const double **matrix0, const double **matrix1, int n0, int n1, int n2);
  int blasDSWAP(int n, double *dx, int incx, double *dy, int incy);
  int blasDAXPY(int n, double da, double *dx, int incx, double *dy, int incy);
  void blasDSCAL(int n, double sa, double x[], int incx);
  int linpackDCHDC(double a[], int lda, int p, double work[], int ipvt[]);
  double **choleskyDecomposition(double **matrix, int n);
  int linpackDTRDI(double t[], int ldt, int n);
  double **matrixLowerTriangularInverse(double **matrix, int n);
  
}
