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

	int M = 1/h + 1;
	int M_2 = M - 2;
	int slit = atoi(argv[11]); //number of slits
	mat V = set_pot (M_2, h, slit, v_0);	//set potential matrix
	cx_mat U_0 (M_2, M_2, fill::zeros);	// initial state
	cx_mat U_n_mat (M_2, M_2, fill::zeros);
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
	U_0_flat = U_0_flat/norm(U_0_flat);
	b = B*U_0_flat;
	U_0 = reshape(U_0_flat, M_2, M_2);
	cube probability (M_2, M_2, T/dt + 1);
	
	cx_mat U_0_c = conj(U_0);		//compute complex conjugate
	probability.slice(0) = real(U_0_c % U_0);		//Schur product
	probability.slice(0).save("probability0_"+std::to_string(slit)+".bin");
	U_0.save("U0_state_"+std::to_string(slit)+".bin");

	for (int j = 1; j <= T/dt; j ++) {
		//solve the matrix equation
		U_n = spsolve(A, b);
		b = B*U_n;
		U_n_mat = reshape(U_n, M_2, M_2);
		probability.slice(j) = real(conj(U_n_mat)%U_n_mat);
		if (j*dt == 0.001){
			probability.slice(j).save("probability1_"+std::to_string(slit)+".bin");
			U_n_mat.save("U1_state_"+std::to_string(slit)+".bin");
		}
		if (j*dt == 0.002){
			probability.slice(j).save("probability2_"+std::to_string(slit)+".bin");
			U_n_mat.save("U2_state_"+std::to_string(slit)+".bin");
		}
	}
	
	probability.save("probability_"+std::to_string(slit)+".bin");

}
