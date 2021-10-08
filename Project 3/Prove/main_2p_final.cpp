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
   

double n_step = 1.0e4;
double t_tot = 100.;
double dt = t_tot/n_step;
int n_particles=2;

mat position1_evol_inter (n_step,3);
mat position2_evol_inter (n_step,3);
mat position1_evol_nointer (n_step,3);
mat position2_evol_nointer (n_step,3);
vec position_0 (3);
vec velocity_0 (3);
mat  r_step (3,n_particles);
mat v_step (3,n_particles);

PenningTrap my_trap_inter = PenningTrap(9.65e1, 9.65e8, 1.0e4,1);
PenningTrap my_trap_nointer = PenningTrap(9.65e1, 9.65e8, 1.0e4,0);
Particle my_particle = Particle(position_0, velocity_0, 1, 40.08);

for (int k=0; k<n_particles;k++) { 
 	position_0.randu()*1.0e4;
 	velocity_0.randu();

my_particle.position_=position_0;
my_particle.velocity_=velocity_0;
my_trap_inter.add_particle(my_particle);
my_trap_nointer.add_particle(my_particle);

cout<<"initial positions: "<< my_particle.position()<<endl;

}



// with interactions 
cube R (3, n_step, my_trap_inter.particle_collection.size(), fill::zeros);

	for (int i=0; i<n_step; i++) {
		for (int j=0; j<my_trap_inter.particle_collection.size(); j++){									//loop over particles
			if (i==0){
				R.slice(j).col(i)=my_trap_inter.particle_collection[j].position();
				r_step.col(j)=my_trap_inter.particle_collection[j].position();
				v_step.col(j)=my_trap_inter.particle_collection[j].velocity();
			}
			else {
				my_trap_inter.evolve_RK4(dt,j, r_step, v_step);
				R.slice(j).col(i)=r_step.col(j);
			}
		}
		myfile_t<<i*dt<<endl;
		}
	
	ofstream myfile1_inter;
	myfile1_inter.open ("xy1_inter.txt");
	myfile1_inter<<R.slice(0).t();
	myfile1_inter.close();
	
	ofstream myfile2_inter;
	myfile2_inter.open ("xy2_inter.txt");
	myfile2_inter<<R.slice(1).t();
	myfile2_inter.close();

// no interactions
	for (int i=0; i<n_step; i++) {
		for (int j=0; j<my_trap_nointer.particle_collection.size(); j++){									//loop over particles
			if (i==0){
				R.slice(j).col(i)=my_trap_nointer.particle_collection[j].position();
				r_step.col(j)=my_trap_nointer.particle_collection[j].position();
				v_step.col(j)=my_trap_nointer.particle_collection[j].velocity();
			}
			else {
				my_trap_nointer.evolve_RK4(dt,j, r_step, v_step);
				R.slice(j).col(i)=r_step.col(j);
			}
		}
		}
	
	
	ofstream myfile1_nointer;
	myfile1_nointer.open ("xy1_nointer.txt");
	myfile1_nointer<<R.slice(0).t();
	myfile1_nointer.close();
	
	ofstream myfile2_nointer;
	myfile2_nointer.open ("xy2_nointer.txt");
	myfile2_nointer<<R.slice(1).t();
	myfile2_nointer.close();
	
	
	myfile_t.close();
	



return 0;

}