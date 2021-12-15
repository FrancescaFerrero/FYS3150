#include "utils.hpp"

using namespace std;
using namespace arma;

/*
translates a pair of indices (i,j) of the matrix into a single index (k), that gives the corresponding position in a vector 
*/
int idx (int i, int j, int M_2) { 
	int k = 0;
	k = j*M_2 + i%(M_2);
	return k;
}

/*
Fill the imaginary part of the two matrices A and B passed by reference
*/
void fill_AB (int M_2, double h, double dt, mat V, sp_mat &A_im, sp_mat &B_im ) {

	int M_squared = M_2*M_2; 
	double r = 0.5*dt/(h*h);
	vec r_2 (M_2 - 1, fill::value(r));		//for the submatrix on the diagonal
	vec r_6 (M_squared - M_2, fill::value(r));	
//defining imaginary part of vectors a and b
	vec v = V.as_col();
	vec a_im (M_squared);
	for (int j = 0; j < M_squared; j++){
		a_im(j) = 4*r + 0.5*dt*v(j);
	}
	vec b_im = - a_im;
// defining the submatrix 3x3
	mat R_square (M_2, M_2, fill::zeros);
	R_square.diag(1) = r_2;
	R_square.diag(-1) = r_2;
// building the imaginary part of matrices A and B
	for (int i = 0; i <= M_squared - M_2; i = i + M_2){
		A_im(span(i,i+M_2-1), span(i,i+M_2-1)) = - R_square;
		B_im(span(i,i+M_2-1), span(i,i+M_2-1)) = R_square;
	}
	A_im.diag(-M_2) = - r_6;
	A_im.diag(M_2) = - r_6;
	A_im.diag() = a_im;
	B_im.diag(-M_2) = r_6;
	B_im.diag(M_2) = r_6;
	B_im.diag() = b_im;
}	

/*
create a matrix of the initial conditions of the Gaussian wavepacket
*/
void Gauss_pack (int M_2, double h, double x_c, double y_c, double sigma_x, double sigma_y, double p_x, double p_y, cx_mat &U_0){
	cx_double im_unit (0, 1.0);
	
	for (int i = 0; i < M_2; i++){
		for (int j = 0; j < M_2; j++){ 
			U_0(i,j) = exp(-((j*h-x_c)*(j*h-x_c)/(2*sigma_x*sigma_x)) - ((i*h-y_c)*(i*h-y_c)/(2*sigma_y*sigma_y)) + im_unit*p_x*(j*h-x_c) + im_unit*p_y*(i*h-y_c));
		}
	}
}


/*
Initialize potential V
*/
mat set_pot (int M_2, double h, int slit, double v_0) {
	mat V (M_2, M_2, fill::zeros);	

	int n_thickness = 0.02/h;	//length in number of cells of the wall thickness
	//column index where the wall starts
	int idx_wall_x = (M_2/2 - n_thickness/2) - 1; // -1 because the first index in matrix is 0 

	int n_aperture = 0.05/h;	//length in number of cells of the slit aperture
	int n_central_wall = 0.05/h;	//length in number of cells of the wall separating two slits
	int n_wall = (M_2 - slit*n_aperture - (slit - 1)*n_central_wall)/2; //length in number of cells of the upper and lower wall
	//row index where the wall ends
	int idx_wall_y_up  = n_wall -1;
	int idx_wall_y_down = M_2 - n_wall - 1;

	if (slit == 1){
		//add upper and lower walls 
		V(span(0, idx_wall_y_up), span(idx_wall_x, idx_wall_x + n_thickness)).fill(v_0); 
		V(span(idx_wall_y_down, M_2 - 1), span(idx_wall_x, idx_wall_x + n_thickness)).fill(v_0);
	}
	if (slit == 2){
		//add upper and lower walls
		V(span(0, idx_wall_y_up), span(idx_wall_x, idx_wall_x + n_thickness)).fill(v_0); 
		V(span(idx_wall_y_down, M_2 - 1), span(idx_wall_x, idx_wall_x + n_thickness)).fill(v_0);
		//add one central wall
		V(span(idx_wall_y_up + n_aperture, idx_wall_y_up + n_aperture + n_central_wall), span(idx_wall_x, idx_wall_x + n_thickness)).fill(pow(10,10));
	}
		if (slit == 3){
		//add upper and lower walls
		V(span(0, idx_wall_y_up), span(idx_wall_x, idx_wall_x + n_thickness)).fill(v_0); 
		V(span(idx_wall_y_down, M_2 - 1), span(idx_wall_x, idx_wall_x + n_thickness)).fill(v_0);
		//add 2 central walls
		V(span(idx_wall_y_up + n_aperture, idx_wall_y_up + n_aperture + n_central_wall), span(idx_wall_x, idx_wall_x + n_thickness)).fill(v_0);
		V(span(idx_wall_y_up + 2*n_aperture + n_central_wall, idx_wall_y_up + 2*n_aperture + 2*n_central_wall), span(idx_wall_x, idx_wall_x + n_thickness)).fill(v_0);
	}
	return V;
}
