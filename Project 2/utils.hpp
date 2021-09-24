#ifndef __utils_hpp__
#define __utils_hpp__
 
#include <armadillo>
#include <sstream>
#include <string>
#include <iostream>
#include <cstring>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <assert.h>
#include <iomanip>


double max_offdiag_symmetric(const arma::mat& A, int& k, int& l);
arma::mat create_tridiagonal(int N, double a, double d);
void jacobi_rotate(arma::mat& A, arma::mat& R, int k, int l);
void jacobi_eigensolver (arma::mat& A, double eps, arma::vec& eigenvalues, arma::mat& eigenvectors, const int maxiter, int& iterations, bool& converged);


#endif  // end of include guard 
