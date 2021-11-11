#include "utils.hpp"

using namespace std;
using namespace arma;

int main (int argc, char* argv[]){ 
	
	int L = atoi(argv[1]);		// set dimension of the LxL matrix
	double T = atof(argv[2]);	// set temperature
	//arma_rng::set_seed_random();
	
	// to store data in a file
	ofstream myfile;
	// set ordered or random initial state
	string state = argv[3];
	imat spin_matrix;
	if (state == "ordered"){
		spin_matrix = ones<imat>(L,L);
		myfile.open ("6_ord_"+std::to_string((int)T)+"_"+std::to_string(L)+".txt");
	}
	else{
		spin_matrix = randi<imat>(L, L, distr_param(0,1))*2 - 1;
		myfile.open ("6_rand_"+std::to_string((int)T)+"_"+std::to_string(L)+".txt");
	}
	
	int max_cycles = 100000;		// number of MC cycles
	vec average (4, fill::zeros); 	// initialize a vector to store averages
	int burning_time = 200000;	
	
	// initial energy and magnetization
	double E = energy(spin_matrix);
	double M = mag(spin_matrix);
	
	double Cv; 	 // specific heat capacity 
	
	
	// set up array for the Boltzmann factor
	vec Bf (17, fill::zeros);
	for( int i = -8; i <= 8; i += 4){ 
		Bf(i+8) = exp(-1.*i/T);
	}
	

	for (int i = 0; i < burning_time; i++){
		Metropolis(spin_matrix, E, M, Bf); 
	}

	for (int n_cycles = 1; n_cycles <= max_cycles; n_cycles++){	// loop over the different numbers of cycles
		// Monte Carlo 
		Metropolis(spin_matrix, E, M, Bf);
			// update averages
		average(0) += E; 
		average(1) += E*E; 
		average(2) += fabs(M); 
		average(3) += M*M;

		myfile<<average(0)*(1./(n_cycles*L*L))<<endl;

	}
	
	myfile.close();
	
	average(0) = average(0)*(1./(max_cycles*L*L));
	average(1) = average(1)*(1./(max_cycles*L*L));

	// Compute specific heat capacity and susceptibility
	Cv = (average(1) - pow(average(0),2)*(L*L))/(T*T);

	cout<<"Cv = "<<Cv<<endl;

	
	return 0;
}
