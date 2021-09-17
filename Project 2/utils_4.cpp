#include"utils_4.hpp"

using namespace std;
using namespace arma;

// A function that finds the max off-diag element of a symmetric matrix A.
// - The matrix indices of the max element are returned by writing to the  
//   int references k and l (row and column, respectively)
// - The value of the max element A(k,l) is returned as the function
//   return value
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