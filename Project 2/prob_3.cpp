#include <iostream>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <armadillo>
#include <cmath>
#include <cstdlib>

using namespace std;
using namespace arma;

arma::mat create_tridiagonal(int N, double a, double d);

int main(){
  
int n = 7;
int N = n-1;
double h = 1./n; 
double a = -1./(pow(h, 2));
double d = 2./(pow(h, 2));


// Calculate analytical eigenvalues and eigenvectors
vec eigvals_analytic(N);
mat eigvecs_analytic(N,N);    


for(int i=1; i<=N; i++){
 for(int j=1; j<=N; j++){
	eigvecs_analytic(j-1,i-1) = sin((M_PI*i*j)/(N+1));
 }
	eigvals_analytic(i-1) = d + 2*a*cos((M_PI*i)/(N+1));
}

eigvals_analytic = normalise(eigvals_analytic);
eigvecs_analytic = normalise(eigvecs_analytic);

eigvecs_analytic.print("eigvecs = ");
eigvals_analytic.t().print("eigvals = ");


// Create tridiagonal matrix
mat A = create_tridiagonal(N, a, d);
cout << "A = " << endl << A << endl;
cout << endl << endl;


// Solve matrix equation with armadillo
vec eigvals_arma(N);
mat eigvecs_arma(N,N);

eig_sym(eigvals_arma, eigvecs_arma, A);

eigvals_arma = normalise(eigvals_arma);
eigvecs_arma = normalise(eigvecs_arma);

eigvecs_arma.print("eigvecs = ");
eigvals_arma.t().print("eigvals = ");


  return 0;
}


// Define a function to create a tridiagonal matrix
arma::mat create_tridiagonal(int N, double a, double d)
{
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
