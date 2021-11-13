#include "utils.hpp"
#include "omp.h"  // OpenMP header

using namespace std;
using namespace arma;

int main (int argc, char* argv[]){ 
	
	int L = atoi(argv[1]);	// set dimension of the LxL matrix
	double min_T = 2.1;
	double max_T = 2.4;
	double step_size = 0.003; //(max_T - min_T)/100
	int burnin_time = 1000;	
	int max_cycles = 2097152;		// number of MC cycles 2^21
	int n_step = (max_T - min_T)/step_size + 1;	// n_step points correspond to (n_step - 1) intervals

	// initial energy and magnetization
	double E;
	double M;
	double X;	 // susceptibility
	double Cv; 	 // specific heat capacity 

	mat results = mat(n_step+1, 5, fill::zeros);

	ofstream myfile;
	myfile.open ("7_rand_" + std::to_string(L) + ".txt");

	#pragma omp parallel private(E, M, X, Cv)  // Start parallel region
  	{	
  		const int my_thread = omp_get_thread_num();
  		arma_rng::set_seed_random();

  		double start = omp_get_wtime();
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
			
			if (my_thread == 0) cout << T<<endl;  // to keep track of where we are during execution

			// compute final average and normalize to the number of spins 
			double avg_e = average(0)/double(L*L);
			avg_e = avg_e/(double)(max_cycles);
			double avg_e2 = average(1)/double(L*L);
			avg_e2 = avg_e2/(double)(max_cycles);
			double avg_m = average(2)/double(L*L);
			avg_m = avg_m/(double)(max_cycles);
			double avg_m2 = average(3)/double(L*L);
			avg_m2 = avg_m2/(double)(max_cycles);

			// Compute specific heat capacity and susceptibility per spin
			Cv = (avg_e2 - pow(avg_e,2)*(L*L))/(T*T);
			X = (avg_m2 - pow(avg_m,2)*(L*L))/T;

			results(int_T,0) = T;
      		results(int_T,1) = avg_e;	//avg energy per spin
      		results(int_T,2) = avg_m;	//avg magnetization per spin
      		results(int_T,3) = Cv;
      		results(int_T,4) = X;

		} // End parallelized loop over T
		double end = omp_get_wtime();
		// to know how long it takes __ every thread is doing the same type of operations, so they will basically take the same time.
		double timeused = end - start;
		if (my_thread == 0) cout << "timeused = " << timeused << endl;
	}	// End entire parallel region

	myfile << results;
	myfile.close();

	return 0;
}
