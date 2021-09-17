#include"utils_5.hpp"

using namespace std;
using namespace arma;


// Define a function to create a tridiagonal matrix
arma::mat create_tridiagonal(int N, double a, double d) {
  // Start from identity matrix
  arma::mat A = arma::mat(N, N, arma::fill::eye);

  // Fill the first row
  A(0,0) = d;
  A(0,1) = a;

  // Loop that fills rows 2 to n-1 
  for (int i=1; i<N-1; i++){
    for (int j=0; j<=N-1; j++){ 

      if(i==j){
         A(i, j) = d;
      }
      if(j == i+1){
         A(i, j) = a;
      }
      if(j == i-1){
         A(i, j)= a;

      }
    } 
  }

  // Fill last row 
  A(N-1, N-2) = a;
  A(N-1, N-1) = d;

  return A;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// A function that finds the max off-diag element of a symmetric matrix A.

double max_offdiag_symmetric(const arma::mat& A, int& k, int& l)
{
  // Get size of the matrix A. 
  
  int size_A =  A.n_rows;
  
  // Possible consistency checks:
  // Check that A is square and larger than 1x1. Here you can for instance use A.is_square(), 
  assert (A.is_square()== true && size_A>1);

  // Initialize references k and l to the first off-diagonal element of A
  k= 0;
  l=1;

  // Initialize a double variable 'maxval' to A(k,l). We'll use this variable 
  // to keep track of the largest off-diag element.
  
  double maxval = abs (A(k,l));

  // Loop through all elements in the upper triangle of A (not including the diagonal)
  // When encountering a matrix element with larger absolute value than the current value of maxval,
  // update k, l and max accordingly.
  
  for (int i=1; i< size_A; i++) {
  	for (int j=0; j<i; j++) {
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
	double epsilon= 1.0e-8;
	double tau;
	double t,c,s;
	double A_m;
	double R_m;
	double maxvalue;
	
	while (abs(A(k,l))> epsilon) {
		if (A(k,l))==0) {
			t=0;
			s=0;
			c=1;
		}
		else { 
			tau= (A(l,l) - A(k,k))/(2*A(k,l));
			if (tau>=0) {
				t = - tau + sqrt(1+pow(tau,2));
			}
			else {
				t =  - tau - sqrt(1+pow(tau,2));
			}
			c=1/sqrt(1+pow(t,2));
			s=c*t;
			}
		
		A_m=A(k,k);	
		A(k,k)= A_m*pow(c,2) -2*A(k,l)*c*s +A(l,l)*pow(s,2);
		A(l,l)= A(l,l)*pow(c,2) + 2*A(k,l)*c*s +A_m*pow(s,2);
		A(k,l)=A(l,k)=0;
		
		
		for (i=0; i<N; i++) {
			if (i != k && i != l) {
				A_m=A(i,k); 
				A(i,k) = A_m*c - A(i,l)*s;
				A(k,i) = A(i,k);
				A(i,l)=A(i,l)*c + A_m*s;
				A(l,i)=A(i,l);
			}
			
			R_m= R(i,k);
			R(i,k)= R_m*c - R(i,l)*s;
			R(i,l)=R(i,l)*c + R_m*s;
		}
		
		maxvalue = max_offdiag_symmetric(A, k, l);
		
	}
		
				
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Jacobi method eigensolver:
// - Runs jacobo_rotate until max off-diagonal element < eps
// - Writes the eigenvalues as entries in the vector "eigenvalues"
// - Writes the eigenvectors as columns in the matrix "eigenvectors"
//   (The returned eigenvalues and eigenvectors are sorted using arma::sort_index)
// - Stops if it the number of iterations reaches "maxiter"
// - Writes the number of iterations to the integer "iterations"
// - Sets the bool reference "converged" to true if convergence was reached before hitting maxiter

void jacobi_eigensolver (const arma::mat& A, double eps, arma::vec& eigenvalues, arma::mat& eigenvectors, const int maxiter, int& iterations, bool& converged);		
	
	
	
	
	
	
	
	
	
	
	
	
	