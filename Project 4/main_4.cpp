#include "utils.hpp"

using namespace std;
using namespace arma;

int main (){ 
	
	int L = 2;		// dimension of the LxL matrix
	arma_rng::set_seed_random();
	imat spin_matrix =  randi<imat>(L, L, distr_param(0,1))*2 - 1;	// create a random configuration of spin
	double T = 1;
	int n_cycles = 10;		// number of MC cycles
	vec average (4, fill::zeros); 	// initialize a vector to store averages
	
	// initial energy and magnetization per spin
	double E = energy_spin(spin_matrix);
	double M = mag_spin(spin_matrix);
	
	// set up array for the Boltmann factor
	vec Bf (17, fill::zeros);
	for( int i = -8; i <= 8; i += 4){ 
		Bf(i+8) = exp(-i/T);
	}
	
	// Monte Carlo 
	for (int i = 1; i <= n_cycles; i++){
		Metropolis(spin_matrix, E, M, Bf);
		// update averages
		average(0) += E; 
		average(1) += pow(E,2); 
		average(2) += M; 
		average(3) += pow(M,2);
	}
	
	average = average*(1/n_cycles)
	
	cout<< average(0)<< endl;


	return 0;

}