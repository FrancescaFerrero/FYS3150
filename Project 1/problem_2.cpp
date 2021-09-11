#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <fstream>
using namespace std;

double function_u(double x);

int main()
{

  // Parameters for output formatting
  int width = 18;
  int prec  = 10;
  
  // Output a header 
  std::cout << "#" << std::setw(width-1) << "points x"
            << std::setw(width) << "function u"
            << std::endl;


int n = 1000;  // n number of points in the grid
double h = 1./(n+1);   // stepsize
double x[n];
double u[n];

for(int i = 0; i<= n; i++){
	x[0] = 0;
	x[i+1] = x[i] + h; 
	u[i] = function_u(x[i]);

	// Output to screen
	std::cout << std::setw(width) << std::setprecision(prec) << std::scientific << x[i] // print the elements of x
              << std::setw(width) << std::setprecision(prec) << std::scientific << u[i] // print the elements of u(x)
              << std::endl;

// Writing the solution on a file 

	ofstream myfile;
	myfile.open ("exact_solutions.txt");
	
   if (!myfile ) { // file couldn't be opened
      cerr << "Error: file could not be opened" << endl;
      exit(1);
   }
	
	for (int i=0; i<=n; i++) {
		myfile << x[i]<<" "<<u[i]<<"\n"; 
	}

	myfile.close();


}

return 0;
}

double function_u(double x){
	return 1 - (1 - exp(-10))*x - exp(-10*x);
}
