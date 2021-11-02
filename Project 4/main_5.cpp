#include "utils.hpp"

using namespace std;
using namespace arma;

int main (){ 
	
	int L = 20;		// dimension of the LxL matrix
	arma_rng::set_seed_random();
	imat spin_matrix =  randi<imat>(L, L, distr_param(0,1))*2 - 1;	// create a random configuration of spin
	//imat spin_matrix = ones<imat>(L,L);		// ordered initial state 
	double T = 1;
	int max_cycles = 100;		// number of MC cycles
	vec average (4, fill::zeros); 	// initialize a vector to store averages
	
	// initial energy and magnetization per spin
	double E = energy_spin(spin_matrix);
	double M = mag_spin(spin_matrix);
	
	// set up array for the Boltmann factor
	vec Bf (17, fill::zeros);
	for( int i = -8; i <= 8; i += 4){ 
		Bf(i+8) = exp(-i/T);
	}
	
	//to store data in a file T=1
	ofstream myfile;
	myfile.open ("5a_rand_1.txt");
//	myfile.open ("5a_ord_1.txt");
	
	for (int n_cycles = 1; n_cycles <= max_cycles; n_cycles++){	// loop over the different numbers of cycles
		// Monte Carlo 
		for (int i = 1; i <= n_cycles; i++){
			Metropolis(spin_matrix, E, M, Bf);
			// update averages
			average(0) += E; 
			average(1) += pow(E,2); 
			average(2) += M; 
			average(3) += pow(M,2);
		}	
	average = average*(1./n_cycles);
	myfile<<n_cycles<<" "<<average(0)<<" "<<average(2)<<endl;
	}
	
	myfile.close();
	
	//////////////////////////////////////////////////////////////////////////////////////
	//to store data in a file T=2.4
	T = 2.4;
	ofstream myfile1;
	spin_matrix =  randi<imat>(L, L, distr_param(0,1))*2 - 1;
	myfile1.open ("5a_rand_2.txt");
//	spin_matrix.ones();
//	myfile1.open ("5a_ord_2.txt");
	
	// Boltzmann factor
	for( int i = -8; i <= 8; i += 4){ 
		Bf(i+8) = exp(-i/T);
	}
	
	E = energy_spin(spin_matrix);
	M = mag_spin(spin_matrix);
	average.zeros();
	
	for (int n_cycles = 1; n_cycles <= max_cycles; n_cycles++){	// loop over the different numbers of cycles
		// Monte Carlo 
		for (int i = 1; i <= n_cycles; i++){
			Metropolis(spin_matrix, E, M, Bf);
			// update averages
			average(0) += E; 
			average(1) += pow(E,2); 
			average(2) += M; 
			average(3) += pow(M,2);
		}	
	average = average*(1./n_cycles);
	myfile1<<n_cycles<<" "<<average(0)<<" "<<average(2)<<endl;
	}
	
	myfile1.close();

	return 0;

}