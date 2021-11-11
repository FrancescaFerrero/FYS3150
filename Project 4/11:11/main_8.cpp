#include "utils.hpp"

using namespace std;
using namespace arma;

int main (int argc, char* argv[]){ 
	
	//arma_rng::set_seed_random();
	int L = atoi(argv[1]);	// set dimension of the LxL matrix
	double T;		// set temperature
	int n_cycles = 1000;	//number of MC cycles
	vec average (4, fill::zeros); 	// initialize a vector to store averages	
	imat spin_matrix;	// lattice of spins

	// initial energy and magnetization
	double E = energy(spin_matrix);
	double M = mag(spin_matrix);	
	double X;	 // susceptibility
	double Cv; 	 // specific heat capacity 

	// to store data in a file
	ofstream myfile;
	// set ordered or random initial state
	string state = argv[2];
	if (state == "ordered"){
		spin_matrix = ones<imat>(L,L);
		myfile.open ("8_ord_"+std::to_string(L)+".txt");
	}
	else{
		spin_matrix = randi<imat>(L, L, distr_param(0,1))*2 - 1;
		myfile.open ("8_rand_"+std::to_string(L)+".txt");
	}

	
for(T = 2.1; T <= 2.4;  T+=0.01){	// loop over the temperature values	
	// set up array for the Boltzmann factor
	vec Bf (17, fill::zeros);
	for( int i = -8; i <= 8; i += 4){ 
		Bf(i+8) = exp(-1.*i/T);
	}
	if (state == "ordered"){
		spin_matrix = ones<imat>(L,L);
	}
	else{
		spin_matrix = randi<imat>(L, L, distr_param(0,1))*2 - 1;
	}
	average.zeros();
	E = energy(spin_matrix);
	M = mag(spin_matrix);
		
	// Monte Carlo 
	for (int i = 1; i <= n_cycles; i++){
		Metropolis(spin_matrix, E, M, Bf);
		// update averages
		average(0) += E; 
		average(1) += E*E; 
		average(2) += fabs(M); 
		average(3) += M*M;
	}		
	// compute final average normalize to the number of spins 
	average = average*(1./(n_cycles*L*L));
	average(1) = average(1)/(L*L);
	average(3) = average(3)/(L*L);		

	// Compute specific heat capacity and susceptibility
	Cv = (average(1) - pow(average(0),2))/(T*T);
	X = (average(3) - pow(average(2),2))/T;

	myfile << T << " " << average(0) << " " << average(2) << " " << Cv << " " << X << endl;
	
}

	myfile.close();
		
	return 0;
}
