#include "utils_5.hpp"

using namespace std;
using namespace arma;


  int N;
  arma::mat A = arma::mat(N, N, arma::fill::randu);

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// A function that finds the max off-diag element of a symmetric matrix A

double max_offdiag_symmetric(const arma::mat& A, int& k, int& l)
{

  // Get size of the matrix A. 
  int size_A =  A.n_rows;
  
  // Possible consistency checks:
  // Check that A is square and larger than 1x1. Here you can for instance use A.is_square(), 
  assert (A.is_square()== true && size_A>1);

  // Initialize references k and l to the first off-diagonal element of A
  k=0;
  l=1;

  // Initialize a double variable 'maxval' to A(k,l). We'll use this variable 
  // to keep track of the largest off-diag element.
  
  double maxval = abs (A(k,l));

  // Loop through all elements in the upper triangle of A (not including the diagonal)
  // When encountering a matrix element with larger absolute value than the current value of maxval,
  // update k, l and max accordingly.
  
  for (int i=1; i< size_A; i++){
  	for (int j=0; j<i; j++){
  		if (abs (A(j,i)) > maxval ) {
  			k=j;
  			l=i;
  			maxval= abs(A(j,i));
  		}
  	}
  }

  return maxval;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Performs a single Jacobi rotation, to "rotate away"
// the off-diagonal element at A(k,l).
// - Assumes symmetric matrix, so we only consider k < l
// - Modifies the input matrices A and R

void jacobi_rotate(arma::mat& A, arma::mat& R, int k, int l) {
	double tau;
	double t,c,s;
	double A_m;
	double R_m;
	double maxvalue;

		if (A(k,l)==0) {
			t=0;
			s=0;
			c=1;
		}
		else { 
			tau = (A(l,l) - A(k,k))/(2*A(k,l));
			if (tau>=0) {
				t = - tau + sqrt(1+pow(tau,2));
			}
			else {
				t = - tau - sqrt(1+pow(tau,2));
			}
			c=1/sqrt(1+pow(t,2));
			s=c*t;
		}

		A_m = A(k,k);	
		A(k,k) = A_m*pow(c,2) -2*A(k,l)*c*s +A(l,l)*pow(s,2);
		A(l,l) = A(l,l)*pow(c,2) + 2*A(k,l)*c*s +A_m*pow(s,2);
		A(k,l) = A(l,k)=0;


		for (int i=0; i<A.n_rows; i++) {
			if (i != k && i != l) {
				A_m = A(i,k); 
				A(i,k) = A_m*c - A(i,l)*s;
				A(k,i) = A(i,k);
				A(i,l) = A(i,l)*c + A_m*s;
				A(l,i) = A(i,l);
			}

			R_m = R(i,k);
			R(i,k) = R_m*c - R(i,l)*s;
			R(i,l) = R(i,l)*c + R_m*s;
		}
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Jacobi method eigensolver:
void jacobi_eigensolver (arma::mat& A, double eps, arma::vec& eigenvalues, arma::mat& eigenvectors, const int maxiter, int& iterations, bool& converged) {
    int k;
    int l;
    mat R = mat(A.n_rows, A.n_rows, fill::eye);
    double maxvalue = max_offdiag_symmetric(A, k, l);
    
    // Runs jacobi_rotate until max off-diagonal element < eps
    while (abs(A(k,l))> eps) {
        jacobi_rotate(A, R, k, l) ;
        maxvalue = max_offdiag_symmetric(A, k, l);	
        iterations = iterations + 1;   // Writes the number of iterations to the integer "iterations"
        assert (iterations<maxiter);   // Stops if the number of iterations reaches "maxiter"
    }
    converged = true;
        
    for (int i=0;i<A.n_rows; i++) {
        eigenvalues(i) = A(i,i);
    }
    
    eigenvectors = R;    // Each row of R is an eigenvector
    
    // Writes the eigenvalues as entries in the vector "eigenvalues"
    // Writes the eigenvectors as columns in the matrix "eigenvectors"
    eigenvalues=sort(eigenvalues);
    eigenvectors=sort(eigenvectors);

}
