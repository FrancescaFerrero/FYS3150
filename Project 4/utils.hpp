#ifndef __utils_hpp__
#define __utils_hpp__
 
#include<armadillo>
#include <ctime>
#include<cmath>
#include <iostream>
#include<cstring>
#include<fstream>
#include<iomanip>
#include <cstdlib>


int idx (int i, int limit, int add);
double energy_spin (arma::imat spin_matrix);
double mag_spin (arma::imat spin_matrix);
void Metropolis(arma::imat& spin_matrix, double& E, double& M, arma::vec Bf);


	
	
	
#endif  // end of include guard 