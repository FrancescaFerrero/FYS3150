#include <iostream>
#include<cstring>
#include<fstream>
#include<iomanip>
#include<armadillo>
#include<cmath>
#include <cstdlib>

using namespace std;

int main() { 
	int n=1.0e3;  // n number of points in the grid, excluding the boundary points (0,1)
	arma::vec a = arma::vec(n-1).fill(-1.);
	arma::vec c = arma::vec(n-1).fill(-1.);
	arma::vec b = arma::vec(n).fill(2.);
	arma::vec g = arma::vec(n);
	double h= 1./(n+1);   //stepsize
	
	arma::vec b_tilde = arma::vec(n);
	arma::vec g_tilde = arma::vec(n);
	arma::vec v = arma::vec(n);
	arma::vec x = arma::vec(n);
	
	// fill in the x vector
	for ( int i = 0; i<n ; i++) {
   		x(i)= h*(i+1);
   }
   
	//fill in the g vector 
	for ( int i = 0; i<n ; i++) {
   		g(i)= pow(h,2)*100*exp(-10*x(i));
   }
   
   //forward substitution 
   b_tilde(0)=b(0);
   g_tilde(0)=g(0);
   
	for ( int i = 1; i<n ; i++) {
   		b_tilde(i)= b(i) - a(i-1)*c(i-1)/b_tilde(i-1) ;					// a(i-1) and not a(i) as we wrote in class! Because in class the first 'a' was
   		g_tilde(i)= g(i) - a(i-1)*g_tilde(i-1)/b_tilde(i-1) ;			//denoted as a2
   }
   
//    //backward substitution
   v(n-1)=g_tilde(n-1)/b_tilde(n-1);
  
   for ( int i = n-2; i>=0 ; --i) {
   		v(i)= (g_tilde(i) - c(i)*v(i+1))/b_tilde(i);
    }

//////////////////writing the solution to a file 

	ofstream myfile;
	myfile.open ("solutions_7_1000.txt");
	
   if (!myfile ) { // file couldn't be opened
      cerr << "Error: file could not be opened" << endl;
      exit(1);
   }
	
	for (int i=0; i<n; i++) {
		myfile << x(i)<<" "<<v(i)<<"\n"; 
	}
	
	myfile.close();
  
  
  
  
  return 0;
}