#include <iostream>
#include<cstring>
#include<fstream>
#include<iomanip>
#include<armadillo>
#include<cmath>
#include <time.h>
#include <cstdlib>

using namespace std;

int main() { 

// opening file

	ofstream myfile;
	myfile.open ("times_general_alg.txt");
	
   if (!myfile ) { // file couldn't be opened
      cerr << "Error: file could not be opened" << endl;
      exit(1);
   }
   
   ofstream myfile2;
	myfile2.open ("times_special_alg.txt");
	
   if (!myfile2 ) { // file couldn't be opened
      cerr << "Error: file could not be opened" << endl;
      exit(1);
   }
	
	
////GENERAL ALGORITHM ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	int n=10;  // n number of points in the grid, excluding the boundary points (0,1)
	
	
	while (n<=1.0e6){
	
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
	
		double sum_t=0; 
	

		for (int k=0;k<10;k++){

// Start measuring time
			clock_t t1 = clock();
   
 //forward substitution 
 			b_tilde(0)=b(0);
 			g_tilde(0)=g(0);
   
			for ( int i = 1; i<n ; i++) {
   				b_tilde(i)= b(i) - a(i-1)*c(i-1)/b_tilde(i-1) ;					
   				g_tilde(i)= g(i) - a(i-1)*g_tilde(i-1)/b_tilde(i-1) ;			//denoted as a2
   			}
   
//backward substitution
   			v(n-1)=g_tilde(n-1)/b_tilde(n-1);
  
   			for ( int i = n-2; i>=0 ; --i) {
   				v(i)= (g_tilde(i) - c(i)*v(i+1))/b_tilde(i);
   			}
    
    
// Stop measuring time
    		clock_t t2 = clock();
    
    // Calculate the elapsed time.
  			double duration_seconds = ((double) (t2 - t1)) / CLOCKS_PER_SEC;
    
   			sum_t = sum_t + duration_seconds;
    }
    myfile << n<<" "<<sum_t/10<<"\n"; 
    n=n*10;
  }
    
////SPECIAL ALGORITHM ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
	n=10;  // n number of points in the grid, excluding the boundary points (0,1)
	
	
	while (n<=1.0e6){
	
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
	
		double sum_t=0; 
	

		for (int k=0;k<10;k++){

// Start measuring time
			clock_t t1 = clock();
			
			g_tilde(0)=g(0);
   
			for ( int i = 1; i<n ; i++) {
   				g_tilde(i)= g(i) + 1.*i/(i+1)*g_tilde(i-1) ;			// a(i-1) and not a(i) as we wrote in class! Because in class the first 'a' was denoted as a2
  			 }
  
  //backward substitution
   v(n-1)=g_tilde(n-1)*1.*n/(n+1.);
  
  			 for ( int i = n-2; i>=0 ; --i) {
   				v(i)= (g_tilde(i) + v(i+1)) *1.*(i+1)/(i+2);
    			}  	  
    			
   // Stop measuring time
    		clock_t t2 = clock();
    
    // Calculate the elapsed time.
  			double duration_seconds = ((double) (t2 - t1)) / CLOCKS_PER_SEC;
    
   			sum_t = sum_t + duration_seconds;
    }
    myfile2 << n<<" "<<sum_t/10<<"\n"; 
    n=n*10;
  }
    
    
    myfile.close();
    myfile2.close();
    return 0;
}