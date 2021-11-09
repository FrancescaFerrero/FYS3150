#include "utils.hpp"
#include "omp.h"  // OpenMP header

using namespace std;
using namespace arma;

int main (int argc, char* argv[]){ 
	
	//arma_rng::set_seed_random();
	int L = atoi(argv[1]);	// set dimension of the LxL matrix
	int n_cycles = 1000;	//number of MC cycles
	vec average (4, fill::zeros); 	// initialize a vector to store averages	
	imat spin_matrix;	// lattice of spins
	double T;
	double min_T = 2.1;
	double max_T = 2.4;
	double step_size = 0.01;
	int n_step = (max_T - min_T)/step_size + 1;	// n_step points correspond to (n_step - 1) intervals

	// initial energy and magnetization
	double E = energy(spin_matrix);
	double M = mag(spin_matrix);	
	double X;	 // susceptibility
	double Cv; 	 // specific heat capacity 


	#pragma omp parallel // Start parallel region
  	{
  		// Each thread will get its own output file name
	    const int my_thread = omp_get_thread_num();
	    ofstream myfile;

		// set ordered or random initial state
		string state = argv[2];
		if (state == "ordered"){
			spin_matrix = ones<imat>(L,L);
			myfile.open ("7_ord_"+std::to_string(L)+ ".thread_" + to_string(my_thread) + ".txt", ofstream::trunc); // ofstream::trunc makes sure old content is deleted
		}
		else{
			spin_matrix = randi<imat>(L, L, distr_param(0,1))*2 - 1;
			myfile.open ("7_rand_" + std::to_string(L) + ".thread_" + to_string(my_thread) + ".txt", ofstream::trunc);
		}

	  

		#pragma omp for // Start parallel region
		for(int int_T = 0; int_T <= n_step; int_T += 1){	// loop over the temperature values	
			T = (int_T*step_size) + min_T;
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

		} // End parallelized loop over T

		myfile.close();

	}	// End entire parallel region
		
	return 0;
}
