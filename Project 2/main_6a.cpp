#include "utils.hpp"

using namespace std;
using namespace arma;

int main (){ 

	int N;

	ofstream myfile;
	myfile.open ("iterations_6a.txt");

   if (!myfile ) { // file couldn't be opened
      cerr << "Error: file could not be opened" << endl;
      exit(1);
   }

	double eps = 1.0e-8;
	int maxiter = 1.0e6;

for(N=10; N<=200; N+=10){

	double h = 1./(N+1); 
	double a = -1./(pow(h, 2));
	double d = 2./(pow(h, 2));
	int iterations = 0;
	bool converged = false;

// Solve matrix equation with Jacobi method
	vec eigenvalues(N);
	mat eigenvectors(N,N);

	mat A = create_tridiagonal(N, a, d);

	jacobi_eigensolver (A, eps, eigenvalues, eigenvectors, maxiter, iterations, converged);

// Print the number of iterations
	cout << "number of iterations for N = " << N << ": " << iterations << endl;

// Write resluts on a file 
	myfile << N << " " << iterations <<"\n"; 

}
	myfile.close();
  
	return 0;

}



