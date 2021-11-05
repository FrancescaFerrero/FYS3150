#ifndef __utils_hpp__
#define __utils_hpp__
 
#include <armadillo>
#include <ctime>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <vector>      // For vector
#include <string>      // For string


int idx (int i, int limit, int add);
double energy (arma::imat spin_matrix);
double mag (arma::imat spin_matrix);
void Metropolis(arma::imat& spin_matrix, double& E, double& M, arma::vec Bf);
	
	
	
#endif  // end of include guard 
