#ifndef _LINEAR_ALGEBRA_H_
#define _LINEAR_ALGEBRA_H_

#include <lapacke.h>

/*
 *  Invert dense matrices in row major order of dimension n*n
 */
inline
void matrix_inverse(double* out, const double* in, lapack_int n) {
  lapack_int info(0);
  lapack_int* ipiv = new lapack_int[n];
  
  std::copy(in, in + n * n, out);
  info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, out, n, ipiv);
  if (info != 0)
    throw std::string("matrix is singular.");
  
  info = LAPACKE_dgetri(LAPACK_ROW_MAJOR, n, out, n, ipiv);
  if (info != 0)
    throw std::string("matrix inversion failed.");

  delete [] ipiv;
}

inline
void matrix_inverse(double* m, lapack_int n) {
  lapack_int info(0);
  lapack_int* ipiv = new lapack_int[n];
  
  info = LAPACKE_dgetrf(LAPACK_COL_MAJOR, n, n, m, n, ipiv);
  if (info != 0)
    throw std::string("matrix is singular.");
  
  info = LAPACKE_dgetri(LAPACK_COL_MAJOR, n, m, n, ipiv);
  if (info != 0)
    throw std::string("matrix inversion failed.");

  delete [] ipiv;
}

#endif /* _LINEAR_ALGEBRA_H_ */
