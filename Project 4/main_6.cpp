#include "utils.hpp"

using namespace std;
using namespace arma;

int main (int argc, char* argv[]){ 
	
	int L = atoi(argv[1]);		// set dimension of the LxL matrix
	double T = atof(argv[2]);	// set temperature
	arma_rng::set_seed_random();
	
	// to store data in a file
	ofstream myfile;
	// set random initial state
	imat spin_matrix = randi<imat>(L, L, distr_param(0,1))*2 - 1;
	myfile.open ("6_rand_"+std::to_string((int)T)+"_"+std::to_string(L)+".txt");
	
	int max_cycles = 2000000;		// number of MC cycles
	vec average (4, fill::zeros); 	// initialize a vector to store averages
	int burnin_time = 1000;	
	
	// initial energy and magnetization
	double E = energy(spin_matrix);
	double M = mag(spin_matrix);
	
	// set up array for the Boltzmann factor
	vec Bf (17, fill::zeros);
	for( int i = -8; i <= 8; i += 4){ 
		Bf(i+8) = exp(-1.*i/T);
	}
	
	for (int i = 0; i < burnin_time; i++){
		Metropolis(spin_matrix, E, M, Bf); 
	}

	for (int n_cycles = 1; n_cycles <= max_cycles; n_cycles++){	// do it for <max_cycles> cycles
		// Monte Carlo 
		Metropolis(spin_matrix, E, M, Bf);
		myfile<<E*1./(L*L)<<endl;

	}
	
	myfile.close();
	
	return 0;
}
