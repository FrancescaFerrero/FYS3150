#include "utils.hpp"

using namespace std;
using namespace arma;

int main (int argc, char* argv[]){ 
	int n = atoi(argv[1]);   // to choose  n=10 or n=100
	int N = n-1; 
	double h = 1./(N+1); 
	double a = -1./(pow(h, 2));
	double d = 2./(pow(h, 2));
	double eps = 1.0e-8;
	int maxiter = 1.0e6;
	int iterations = 0;
	bool converged = false;

// Solve matrix equation with Jacobi method
	vec eigenvalues(N);
	mat eigenvectors(N,N);
	vec x(N+1);

	mat A = create_tridiagonal(N, a, d);

	jacobi_eigensolver (A, eps, eigenvalues, eigenvectors, maxiter, iterations, converged);

	eigenvalues = normalise(eigenvalues);
	eigenvectors = normalise(eigenvectors);


for (int i=0; i<=N; i++){
	x(i)= h*i;
}


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



/////////////////////////////////////////////////////////////////////

// Write eigenvectors on a file 
	ofstream myfile1a;
	string filename1 = "eigenvectors_" + to_string(n) + ".txt";
	myfile1a.open (filename1);

   if (!myfile1a) { // file couldn't be opened
      cerr << "Error: file could not be opened" << endl;
      exit(1);
   }

	myfile1a << eigenvectors; 

	myfile1a.close();


//////////////////////////////////////////////////////////////////

// Write eigenvectors on a file 
	ofstream myfile1;
	string filename2 = "eigenvectors_analytical_" + to_string(n) + ".txt";
	myfile1.open (filename2);

   if (!myfile1) { // file couldn't be opened
      cerr << "Error: file could not be opened" << endl;
      exit(1);
   }

	myfile1 << eigvecs_analytic; 

	myfile1.close();


////////////////////////////////////////////////////////////////////

// Write x values on a file 
	ofstream myfile2;
	string filename3 = "x_" + to_string(n) + ".txt";
	myfile2.open (filename3);

   if (!myfile2) { // file couldn't be opened
      cerr << "Error: file could not be opened" << endl;
      exit(1);
   }
   
	for (int i=0; i<=N; i++) {
		myfile2 << x(i) << "\n"; 
	}

	myfile2.close();


	return 0;

}

