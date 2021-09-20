#include "utils_3.hpp"

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
