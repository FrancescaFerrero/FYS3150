#include "utils_5.hpp"

using namespace std;
using namespace arma;

int main (){ 

	int N =6;
	double h = 1./N+1; 
	double a = -1./(pow(h, 2));
	double d = 2./(pow(h, 2));
	double eps= 1.0e-8;
	int maxiter = 1.0e6;
	int iterations=0;
	bool converged= false;
	vec eigenvalues (N);
	mat eigenvectors (N,N);

	mat A = create_tridiagonal(N, a, d);
	
	jacobi_eigensolver (A, eps, eigenvalues, eigenvectors, maxiter, iterations, converged);
	
// Calculate analytical eigenvalues and eigenvectors
	vec eigvals_analytic(N);
	mat eigvecs_analytic(N,N);    


	for(int i=1; i<=N; i++){
	for(int j=1; j<=N; j++){
		eigvecs_analytic(i-1,j-1) = sin((M_PI*i*j)/(N+1));
 	}
		eigvals_analytic(i-1) = d + 2*a*cos((M_PI*i)/(N+1));
	}

	eigvals_analytic = normalise(eigvals_analytic);
	eigvecs_analytic = normalise(eigvecs_analytic);
	eigenvalues = normalise(eigenvalues);
	eigenvectors = normalise(eigenvectors);
	
	cout<< "eigenvalues_jacobi - eigenvalues_analytic = " << endl << eigenvalues.t() - eigvals_analytic.t()<<endl;
	cout<< "eigenvectors_jacobi - eigenvectors_analytic = " << endl << eigenvectors - eigvecs_analytic<<endl;
	cout<<"matrix A: "<<endl<<A<<endl;
	
	
	
	
	return 0;
	
}
	
	
	
	