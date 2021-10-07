#include "PenningTrap.hpp"
#include "Particle.hpp"
#include <algorithm>
#include <ctime>




using namespace std;
using namespace arma;

int main (){ 

double n = 1.0e2;
double t = 1;
double dt = t/n;

mat position_evol_FE (n,3);
mat velocity_evol_FE (n,3);
mat position_evol_RK (n,3);
mat velocity_evol_RK (n,3);
vec position_0 (3);
vec velocity_0 (3);

PenningTrap my_trap = PenningTrap(9.65e1, 9.65e8, 1.0e4);
Particle my_particle = Particle(position_0, velocity_0, 1, 40.08);

for (int k=0; k<3;k++) { 
 	position_0.randu()*1.0e4;
 	velocity_0.randu();

my_particle.position_=position_0;
my_particle.velocity_=velocity_0;
my_trap.add_particle(my_particle);

cout<<"initial positions: "<< my_particle.position()<<endl;

}


	for (int i=0; i<n; i++) {
		for (int j=0; j<my_trap.particle_collection.size(); j++){
		if (i==0){
			position_evol_FE.row(i)=my_trap.particle_collection[j].position().t();
			velocity_evol_FE.row(i)=my_trap.particle_collection[j].velocity().t();
			position_evol_RK.row(i)=my_trap.particle_collection[j].position().t();
			velocity_evol_RK.row(i)=my_trap.particle_collection[j].velocity().t();
		}
		else {
		my_trap.evolve_forward_Euler(dt, j);
		position_evol_FE.row(i)=my_trap.particle_collection[j].position().t();
		velocity_evol_FE.row(i)=my_trap.particle_collection[j].velocity().t();
		
		my_trap.evolve_RK4(dt,j);
		position_evol_RK.row(i)=my_trap.particle_collection[j].position().t();
		velocity_evol_RK.row(i)=my_trap.particle_collection[j].velocity().t();
		}
		
		}
	}


vec F_tot=my_trap.total_force(0,0);

F_tot.print("total force: ");
position_evol_FE.print("positions FE: ");
velocity_evol_FE.print("velocities FE: ");

position_evol_RK.print("positions RK: ");
velocity_evol_RK.print("velocities RK: ");

cout<< "Differences: "<<position_evol_FE-position_evol_RK <<endl;


return 0;

}


/*
srand(unsigned(std::time(nullptr)));
generate(position_0.begin(), position_0.end(), rand);
srand(unsigned(std::time(nullptr)));
generate(velocity_0.begin(), velocity_0.end(), rand);
*/