#include "utils.hpp"

using namespace std;
using namespace arma;

int main (){
  
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

eigvecs_analytic.print("eigvecs_analytic = ");
eigvals_analytic.t().print("eigvals_analytic = ");


// Create tridiagonal matrix 
	mat A = create_tridiagonal(N, a, d);
	cout << "A = " << endl << A << endl;
	cout << endl << endl;


// Solve matrix equation with Armadillo
	vec eigvals_arma(N);
	mat eigvecs_arma(N,N);

	eig_sym(eigvals_arma, eigvecs_arma, A);

	eigvals_arma = normalise(eigvals_arma);
	eigvecs_arma = normalise(eigvecs_arma);

//eigvecs_arma.print("eigvecs_arma = ");
//eigvals_arma.t().print("eigvals_arma = ");


// Check Armadillo vs analytic
	cout << "eigenvalues_arma - eigenvalues_analytic = " << endl << eigvals_arma.t() - eigvals_analytic.t() << endl;
	cout << "eigenvectors_arma - eigenvectors_analytic = " << endl << abs(eigvecs_arma) - abs(eigvecs_analytic) << endl;

	return 0;

}

