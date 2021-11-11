#include "utils.hpp"
#include "omp.h"  // OpenMP header

using namespace std;
using namespace arma;

int main (int argc, char* argv[]){ 
	
	//arma_rng::set_seed_random();
	int L = atoi(argv[1]);	// set dimension of the LxL matrix
	double min_T = 2.1;
	double max_T = 2.4;
	double step_size = 0.01;
	int burnin_time = 200000;	
	int max_cycles = 1000;		// number of MC cycles
	int n_step = (max_T - min_T)/step_size + 1;	// n_step points correspond to (n_step - 1) intervals

	// initial energy and magnetization
	double E;
	double M;
	double X;	 // susceptibility
	double Cv; 	 // specific heat capacity 

	mat results = mat(n_step+1, 5, fill::zeros);

	ofstream myfile;
	myfile.open ("7_rand_" + std::to_string(L) + ".txt");

	#ifdef _OPENMP
	{
		#pragma omp parallel  // Start parallel region
	  	{

			#pragma omp for // Start parallelize loop
			for(int int_T = 0; int_T <= n_step; int_T += 1){	// loop over the temperature values	
				
				double T = (int_T*step_size) + min_T;

				vec average (4, fill::zeros);	// initialize a vector to store averages
				imat spin_matrix = randi<imat>(L, L, distr_param(0,1))*2 - 1;	// set random initial state
				E = energy(spin_matrix);
				M = mag(spin_matrix);
				
				// set up array for the Boltzmann factor
				vec Bf (17, fill::zeros);
				for( int i = -8; i <= 8; i += 4){ 
					Bf(i+8) = exp(-1.*i/T);
				}


				for (int i = 0; i < burnin_time; i++){
					Metropolis(spin_matrix, E, M, Bf); 
				}
					
				// Monte Carlo 
				for (int n_cycles = 1; n_cycles <= max_cycles; n_cycles++){
					Metropolis(spin_matrix, E, M, Bf);
					// update averages
					average(0) += E; 
					average(1) += E*E; 
					average(2) += fabs(M); 
					average(3) += M*M;
				}		

				
				// compute final average and normalize to the number of spins 
				average(0) = average(0)*(1./(max_cycles*L*L));
				average(2) = average(2)*(1./(max_cycles*L*L));
				average(1) = average(1)*(1./(max_cycles*L*L));
				average(3) = average(3)*(1./(max_cycles*L*L));
			

				// Compute specific heat capacity and susceptibility per spin
				Cv = (average(1) - pow(average(0),2)*(L*L))/(T*T);
				X = (average(3) - pow(average(2),2)*(L*L))/T;

				results(int_T,0) = T;
	      		results(int_T,1) = average(0);	//avg energy per spin
	      		results(int_T,2) = average(2);	//avg magnetization per spin
	      		results(int_T,3) = Cv;
	      		results(int_T,4) = X;


			} // End parallelized loop over T
		}	// End entire parallel region
	}


	#else
	{
		imat spin_matrix;	// lattice of spins
	  	vec average (4, fill::zeros);

		for(int int_T = 0; int_T <= n_step; int_T += 1){	// loop over the temperature values	
			
			double T = (int_T*step_size) + min_T;
			
			average.zeros(); 	// initialize a vector to store averages
			spin_matrix = randi<imat>(L, L, distr_param(0,1))*2 - 1;	// set random initial state
			E = energy(spin_matrix);
			M = mag(spin_matrix);

			// set up array for the Boltzmann factor
			vec Bf (17, fill::zeros);
			for( int i = -8; i <= 8; i += 4){ 
				Bf(i+8) = exp(-1.*i/T);
			}


			for (int i = 0; i < burnin_time; i++){
				Metropolis(spin_matrix, E, M, Bf); 
			}
				
			// Monte Carlo 
			for (int n_cycles = 1; n_cycles <= max_cycles; n_cycles++){
				Metropolis(spin_matrix, E, M, Bf);
				// update averages
				average(0) += E; 
				average(1) += E*E; 
				average(2) += fabs(M); 
				average(3) += M*M;
			}		

			
			// compute final average and normalize to the number of spins 
			average(0) = average(0)*(1./(max_cycles*L*L));
			average(2) = average(2)*(1./(max_cycles*L*L));
			average(1) = average(1)*(1./(max_cycles*L*L));
			average(3) = average(3)*(1./(max_cycles*L*L));
		

			// Compute specific heat capacity and susceptibility per spin
			Cv = (average(1) - pow(average(0),2)*(L*L))/(T*T);
			X = (average(3) - pow(average(2),2)*(L*L))/T;

			results(int_T,0) = T;
      		results(int_T,1) = average(0);	//avg energy per spin
      		results(int_T,2) = average(2);	//avg magnetization per spin
      		results(int_T,3) = Cv;
      		results(int_T,4) = X;
		}
	}
	#endif

	myfile << results;
	myfile.close();

	return 0;
	}
