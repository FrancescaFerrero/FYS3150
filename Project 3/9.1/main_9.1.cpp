#include "PenningTrap.hpp"
#include "Particle.hpp"
#include <algorithm>
#include <ctime>

using namespace std;
using namespace arma;

int main (){ 

double n = 1.0e6;
double t = 100;
double dt = t/n;
//std::vector<double> dt = {1.0e-3, 1.0e-4, 1.0e-5, 1.0e-6};
int i;

mat position_evol_FE (n,3);
mat velocity_evol_FE (n,3);
mat position_evol_RK (n,3);
mat velocity_evol_RK (n,3);
vec position_0 (3);
vec velocity_0 (3);

PenningTrap my_trap = PenningTrap(9.65e1, 9.65e8, 1.0e4);
Particle my_particle = Particle(position_0, velocity_0, 1, 40.08);

//for (int k=0; k<3;k++) { 
	
position_0.randu()*1.0e4;
velocity_0.randu();

my_particle.position_=position_0;
my_particle.velocity_=velocity_0;
my_trap.add_particle(my_particle);

// Open file to save position values
ofstream myfile1;
myfile1.open ("r_values.txt");

if (!myfile1 ) { // file couldn't be opened
 cerr << "Error: file could not be opened" << endl;
 exit(1);
}


cout<<"initial positions: "<< my_particle.position()<<endl;

//}


for (i=0; i<n; i++) {
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
		
			my_trap.evolve_RK4(dt,j,0);
			position_evol_RK.row(i)=my_trap.particle_collection[j].position().t();
			velocity_evol_RK.row(i)=my_trap.particle_collection[j].velocity().t();
		}
	
	}
	
	
//cout<< "Time: "<< i*dt <<endl;
}


// Write time and positions on a file 
myfile1 << position_evol_RK;  
//myfile1 << position_evol_FE;  
myfile1.close();


// Open file to save time values
ofstream myfile2;
myfile2.open ("t_values.txt");

if (!myfile2 ) { // file couldn't be opened
 cerr << "Error: file could not be opened" << endl;
 exit(1);
}
for (i=0; i<n; i++){
 myfile2 << i*dt << "\n"; 
}

myfile2.close();


// Relative error
//rel_err = position_evol_RK/


return 0;
}


/*
srand(unsigned(std::time(nullptr)));
generate(position_0.begin(), position_0.end(), rand);
srand(unsigned(std::time(nullptr)));
generate(velocity_0.begin(), velocity_0.end(), rand);
*/
