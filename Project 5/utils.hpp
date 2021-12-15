#ifndef __utils_hpp__
#define __utils_hpp__
 
#include <armadillo>
#include <ctime>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <complex>
#include <cstdlib>
#include <vector>      // For vector
#include <string>      // For string


int idx (int i, int j, int M);
void fill_AB (int M_2, double h, double dt, arma::mat V, arma::sp_mat &A_im, arma::sp_mat &B_im );
arma::mat set_pot (int M_2, double h, int slit, double v_0);
void Gauss_pack (int M_2, double h, double x_c, double y_c, double sigma_x, double sigma_y, double p_x, double p_y, arma::cx_mat &U_0);

	
#endif  // end of include guard 
