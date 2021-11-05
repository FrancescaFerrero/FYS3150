#include "utils.hpp"

using namespace std;
using namespace arma;

// for periodic boundary conditions
int idx (int i, int limit, int add) { 
	return (i + limit + add) % limit;
}


// compute energy 
double energy (imat spin_matrix) {
	double energy = 0;
	int L = arma::size(spin_matrix)(0);
	for(int y = 0; y < L; y++) {	//loop over the whole matrix
		for (int x = 0; x < L; x++){
			energy -= spin_matrix(y,x)*(spin_matrix(idx(y,L,-1),x) + spin_matrix(y,idx(x,L,1))); //taking into account right and bottom neighbours	
		}
	}
	return energy;	//return the energy 
}


// compute magnetization
double mag (imat spin_matrix) {
	double mag = 0;
	int L = arma::size(spin_matrix)(0);
	for(int y = 0; y < L; y++) {	//loop over the whole matrix
		for (int x = 0; x < L; x++){
			mag += spin_matrix(x,y);
		}
	}
	return mag;	//return the magnetization
}


// do a Monte Carlo cycle
 void Metropolis(imat& spin_matrix, double& E, double& M, vec Bf) {
 	int L = arma::size(spin_matrix)(0);
	for (int n = 0; n < L*L; n++){
		// we loop over the whole matrix but we choose the lattice positions x and y randomly
		int ix = randi(distr_param(0, L-1)); 	
		int iy = randi(distr_param(0, L-1)); 
		// compute the energy difference caused by flipping the spins
		int deltaE = 2*spin_matrix(iy,ix)*(spin_matrix(iy, idx(ix,L,-1))+ spin_matrix(idx(iy,L,-1),ix) + spin_matrix(iy,idx(ix,L,1)) + spin_matrix(idx(iy,L,1), ix));
		if(deltaE > 0){
			// the Metropolis test
			if ( randu<double>() <= Bf (deltaE+8) ) {
				spin_matrix(iy, ix) *= -1; // flip the spins
				// update energy and magnetization
				M += (2*spin_matrix(iy,ix));
				E += deltaE;
      		}
      	}
      	else {
      		spin_matrix(iy, ix) *= -1; // flip the spins
      		M += (2*spin_matrix(iy,ix));
      		E += deltaE;      			
      	}
	}
}	 

