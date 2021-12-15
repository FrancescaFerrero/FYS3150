#include "utils.hpp"

using namespace std;
using namespace arma;

int main (int argc, char* argv[]){ 

	// read simulation parameters as command-line input
	//spatial and temporal step sizes
	double h = atof(argv[1]);
	double dt = atof(argv[2]);
	//total time
	double T = atof(argv[3]);

	double x_c = atof(argv[4]);	//coordinate of the centre of the initial wave packet
	double sigma_x = atof(argv[5]);	//initial widths of the wave packet in the x direction
	double p_x = atof(argv[6]);	//wave packet momenta
	double y_c = atof(argv[7]);
	double sigma_y = atof(argv[8]);
	double p_y = atof(argv[9]);

	double v_0 = atof(argv[10]);	//constant potential inside the barrier

	int M = 1/h;
	int M_2 = M - 2;
	int slit = 2; //number of slits
	mat V = set_pot (M_2, h, slit, v_0);	//set potential matrix
	cx_mat U_0 (M_2, M_2, fill::zeros);		// initial state
	Gauss_pack (M_2, h, x_c, y_c, sigma_x, sigma_y, p_x, p_y, U_0);	//Gaussian wave packet

	sp_mat A_im (M_2*M_2, M_2*M_2);
	sp_mat B_im (M_2*M_2, M_2*M_2);

	//imaginary part of A and B
	fill_AB (M_2, h, dt, V, A_im, B_im);

	//put imaginary and real part together
	mat identity = eye(M_2*M_2, M_2*M_2);
	sp_mat identity_sparse = sp_mat(identity);
	sp_cx_mat A (identity_sparse, A_im);
	sp_cx_mat B (identity_sparse, B_im);
	cx_vec U_n (M_2*M_2);	//vector of state
	cx_vec b (M_2*M_2);

	cx_vec U_0_flat = U_0.as_col();
	U_0_flat = U_0_flat/norm(U_0_flat);		//normalize initial state
	b = B*U_0_flat;		//get the first b
	U_0 = reshape(U_0_flat, M_2, M_2);
	mat probability (T/dt + 1, 2);
	probability(0,0) = real(cdot(U_0_flat, U_0_flat));
	probability(0,1) = 0;		//column for time

	for (int j = 1; j <= T/dt; j ++) {
		//solve the matrix equation
		U_n = spsolve(A, b);
		b = B*U_n;
		probability(j,0) = real(cdot(U_n, U_n));
		probability(j,1) = j*dt;
	}

	if (v_0 == 0){ 
		probability.save("probability.txt", arma_ascii);		//without potential
	}
	else{
		probability.save("probability_ds.txt", arma_ascii);
	}
	

}