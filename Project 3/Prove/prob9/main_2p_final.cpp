#include "PenningTrap.hpp"
#include "Particle.hpp"
#include <algorithm>
#include <ctime>


using namespace std;
using namespace arma;

int main (){ 

	ofstream myfile_t;
	myfile_t.open ("time.txt");

   if (!myfile_t ) { // file couldn't be opened
      cerr << "Error: file could not be opened" << endl;
      exit(1);
   }
   
//arma_rng::set_seed_random();
//arma_rng::set_seed(1)
double n_step = 1.0e5;
double t_tot = 100.;
double dt = t_tot/n_step;
int n_particles=2;

vec position_0 (3);
vec velocity_0 (3);
mat  r_step (3,n_particles);
mat v_step (3,n_particles);
cube R (3, n_step, n_particles, fill::zeros);
cube V (3, n_step, n_particles, fill::zeros);

PenningTrap my_trap_inter = PenningTrap(9.65e1, 9.65e8, 1.0e4,1);
PenningTrap my_trap_nointer = PenningTrap(9.65e1, 9.65e8, 1.0e4,0);
Particle my_particle = Particle(position_0, velocity_0, 1, 40.08);

for (int k=0; k<n_particles;k++) { 
 	//position_0= vec(3).randn() * 0.1 * my_trap_inter.d_;
 	//velocity_0= vec(3).randn() * 0.1 * my_trap_inter.d_;
	if (k==0){
		position_0={0,0,0};
	}
	else {
		position_0={10,10,10};
	}
	
	velocity_0.fill(1);

my_particle.position_=position_0;
my_particle.velocity_=velocity_0;
my_trap_inter.add_particle(my_particle);
my_trap_nointer.add_particle(my_particle);

cout<<"initial position particle "<< k+1<<": "<< my_particle.position()<<endl;

}


if (n_particles>1){

// with interactions 

	for (int i=0; i<n_step; i++) {
		for (int j=0; j<my_trap_inter.particle_collection.size(); j++){									//loop over particles
			if (i==0){
				R.slice(j).col(i)=my_trap_inter.particle_collection[j].position();
				V.slice(j).col(i)=my_trap_inter.particle_collection[j].velocity();
				r_step.col(j)=my_trap_inter.particle_collection[j].position();
				v_step.col(j)=my_trap_inter.particle_collection[j].velocity();
			}
			else {
				my_trap_inter.evolve_RK4(dt,j, r_step, v_step);
				//my_trap_inter.evolve_forward_Euler(dt,j, r_step, v_step);
				R.slice(j).col(i)=r_step.col(j);
				V.slice(j).col(i)=v_step.col(j);
			}
		}
		}
	
	ofstream myfile1_inter;
	myfile1_inter.open ("r1_inter.txt");
	myfile1_inter<<R.slice(0).t();
	myfile1_inter.close();

	ofstream myfile2_inter;
	myfile2_inter.open ("r2_inter.txt");
	myfile2_inter<<R.slice(1).t();
	myfile2_inter.close();
	
	ofstream myfile3_inter;
	myfile3_inter.open ("v1_inter.txt");
	myfile3_inter<<V.slice(0).t();
	myfile3_inter.close();
	
	ofstream myfile4_inter;
	myfile4_inter.open ("v2_inter.txt");
	myfile4_inter<<V.slice(1).t();
	myfile4_inter.close();
	
	cout<< "interaction "<< my_trap_inter.interaction_<<endl;	
	cout<< "force inter total "<< my_trap_inter.total_force(0,0)<<endl;
//	cout<< "force inter particle "<< my_trap_inter.force_particle(0,1)<<endl;
}

// no interactions
	for (int i=0; i<n_step; i++) {
		for (int j=0; j<my_trap_nointer.particle_collection.size(); j++){									//loop over particles
			if (i==0){
				R.slice(j).col(i)=my_trap_nointer.particle_collection[j].position();
				V.slice(j).col(i)=my_trap_nointer.particle_collection[j].velocity();
				r_step.col(j)=my_trap_nointer.particle_collection[j].position();
				v_step.col(j)=my_trap_nointer.particle_collection[j].velocity();
			}
			else {
				my_trap_nointer.evolve_RK4(dt,j, r_step, v_step);
				//my_trap_nointer.evolve_forward_Euler(dt,j, r_step, v_step);
				R.slice(j).col(i)=r_step.col(j);
				V.slice(j).col(i)=v_step.col(j);
			}
		}
		myfile_t<<i*dt<<endl;
		}
	
	cout<< "no interaction "<< my_trap_nointer.interaction_<<endl;
	cout<< "force no inter "<< my_trap_nointer.total_force(0,0)<<endl;
	
	ofstream myfile1_nointer;
	myfile1_nointer.open ("r1_nointer.txt");
	myfile1_nointer<<R.slice(0).t();
	myfile1_nointer.close();
	
	ofstream myfile3_nointer;
	myfile3_nointer.open ("v1_nointer.txt");
	myfile3_nointer<<V.slice(0).t();
	myfile3_nointer.close();
	
	if (n_particles>1){
	ofstream myfile2_nointer;
	myfile2_nointer.open ("r2_nointer.txt");
	myfile2_nointer<<R.slice(1).t();
	myfile2_nointer.close();
	
	ofstream myfile4_nointer;
	myfile4_nointer.open ("v2_nointer.txt");
	myfile4_nointer<<V.slice(1).t();
	myfile4_nointer.close();
	}
	
	myfile_t.close();
	



return 0;

}