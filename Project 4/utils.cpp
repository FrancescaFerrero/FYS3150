#include "utils.hpp"

using namespace std;
using namespace arma;

// for periodic boundary conditions
int idx (int i, int limit, int add) { 
	return (i + limit + add) % limit;
}


// compute energy per spin
double energy_spin (imat spin_matrix) {
	double energy_spin = 0;
	int L = arma::size(spin_matrix)(0);
	for(int y = 0; y < L; y++) {			//loop over the whole matrix
		for (int x = 0; x < L; x++){
		energy_spin = - spin_matrix(y, x)*(spin_matrix(idx(y, L, -1), x) + spin_matrix(y, idx(x, L, 1))); //taking into account right and bottom neighbours	
		}
	}
	return energy_spin/(pow(L,2))	;			//return the energy per spin
}


// compute magnetization per spin
double mag_spin (imat spin_matrix) {
	double mag_spin = 0;
	int L = arma::size(spin_matrix)(0);
	for(int y = 0; y < L; y++) {			//loop over the whole matrix
		for (int x = 0; x < L; x++){
		mag_spin += spin_matrix(y, x);
		}
	}
	return abs(mag_spin)/(pow(L,2))	;			//return the absolute value of the magnetization per spin
}


// do a Monte Carlo cycle
 void Metropolis(imat& spin_matrix, double& E, double& M, vec Bf) {
 	int L = arma::size(spin_matrix)(0);
	for(int y =0; y < L; y++) { 				
		for (int x= 0; x < L; x++){
			// we loop over the whole matrix but we choose the lattice positions x and y randomly
			ivec rand_idx =  randi<ivec>(2, distr_param(0, L-1));
			int ix = rand_idx(0); 	
			int iy = rand_idx(1); 
			//compute the energy difference caused by flipping the spin
			int deltaE = 2*spin_matrix(iy,ix)*(spin_matrix(iy, idx(ix,L,-1))+ spin_matrix(idx(iy,L,-1),ix) + spin_matrix(iy,idx(ix,L,1)) + spin_matrix(idx(iy,L,1), ix));
			// the Metropolis test
			if ( randu<double>() <= Bf (deltaE+8) ) {
				spin_matrix(iy, ix) *= -1; // flip the spin
				// update energy and magnetization
				M = abs(M + 2*spin_matrix(iy,ix));
      			E += deltaE;
    		}
		}
 	}
}	 
		
			
	