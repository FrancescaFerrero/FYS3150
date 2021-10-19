#include "PenningTrap.hpp"
#include "Particle.hpp"
#include <algorithm>
#include <ctime>


using namespace std;
using namespace arma;

int main (){ 
   
   	
	ofstream myfile_omega;
	myfile_omega.open ("omega.txt");
	
	
//arma_rng::set_seed_random();

double n_step = 1.0e4;
double t_tot = 500.;
double dt = t_tot/n_step;
double n_particles=100.;
double t=0;
vector<double> p_inside;

vec position_0 (3);
vec velocity_0 (3);
mat  r_step (3,n_particles);
mat v_step (3,n_particles);
mat pos_0 (3, n_particles);
mat vel_0 (3, n_particles);
double V_in=0.0025*9.64852558*1.0e7;


//PenningTrap my_trap_inter = PenningTrap(9.65e1, V_in, 500.,1);
PenningTrap my_trap_nointer = PenningTrap(9.65e1, V_in, 500.,0,0);
Particle my_particle = Particle(position_0, velocity_0, 1, 40.08);

for (int k=0; k<n_particles;k++) { 
 	position_0= vec(3).randn() * 0.1 * my_trap_nointer.d_;
 	velocity_0= vec(3).randn() * 0.1 * my_trap_nointer.d_;

	my_particle.position_=position_0;
	my_particle.velocity_=velocity_0;
	//my_trap_inter.add_particle(my_particle);
	my_trap_nointer.add_particle(my_particle);
	
	pos_0.col(k)=position_0;
	vel_0.col(k)=velocity_0;

	//cout<<"initial position particle "<< k+1<<": "<< my_particle.position()<<endl;
}



// no interactions

for (double om =0.2; om<=2.5; om += 0.002){				//loop over omega_v
	my_trap_nointer.omega_v_= om;
	for (int i=0; i<n_step; i++) {
		for (int j=0; j<my_trap_nointer.particle_collection.size(); j++){									//loop over particles
			if (i==0){
				r_step.col(j)=pos_0.col(j);
				v_step.col(j)=vel_0.col(j);
			}
			else {
				t=dt*(i);
				my_trap_nointer.evolve_RK4(dt,j, r_step, v_step,t);
			}
		}
	}
	p_inside.push_back(my_trap_nointer.count_particles(r_step)/n_particles);
	myfile_omega<<om<<endl;
}
		

	myfile_omega.close();

	
	ofstream myfile1;
	myfile1.open ("p_inside_0.1.txt");
	for (int i = 0; i < p_inside.size(); i++) {
	myfile1<<p_inside.at(i)<<endl;
	}
	myfile1.close();
	

return 0;

}

