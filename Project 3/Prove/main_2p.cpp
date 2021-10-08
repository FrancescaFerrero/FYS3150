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
   
   ofstream myfile_inter;
   ofstream myfile_nointer;
   

double n = 1.0e6;
double t = 100;
double dt = t/n;

mat position_evol_inter (n,3);
mat position_evol_nointer (n,3);
vec position_0 (3);
vec velocity_0 (3);

PenningTrap my_trap_inter = PenningTrap(9.65e1, 9.65e8, 1.0e4,1);
PenningTrap my_trap_nointer = PenningTrap(9.65e1, 9.65e8, 1.0e4,0);
Particle my_particle = Particle(position_0, velocity_0, 1, 40.08);

for (int k=0; k<2;k++) { 
 	position_0.randu()*1.0e4;
 	velocity_0.randu();

my_particle.position_=position_0;
my_particle.velocity_=velocity_0;
my_trap_inter.add_particle(my_particle);
my_trap_nointer.add_particle(my_particle);

cout<<"initial positions interactions: "<< my_particle.position()<<endl;

}

// with interactions 
	for (int i=0; i<n; i++) {
		for (int j=0; j<my_trap_inter.particle_collection.size(); j++){									//loop over particles
			if (i==0){
				position_evol_inter.row(i)=my_trap_inter.particle_collection[j].position().t();
			}
			else {
				my_trap_inter.evolve_RK4(dt,j);
				position_evol_inter.row(i)=my_trap_inter.particle_collection[j].position().t();
			}
			
			string filename = "xy" + to_string(j) + "_inter.txt";
			myfile_inter.open (filename,ios_base::app);					// append instead of overwrite
			myfile_inter<<position_evol_inter;
		}
	}

myfile_inter.close();

// no interactions
	for (int i=0; i<n; i++) {
		for (int j=0; j<my_trap_nointer.particle_collection.size(); j++){									//loop over particles
			if (i==0){
				position_evol_nointer.row(i)=my_trap_nointer.particle_collection[j].position().t();
			}
			else {
				my_trap_nointer.evolve_RK4(dt,j);
				position_evol_nointer.row(i)=my_trap_nointer.particle_collection[j].position().t();
			}
			
			string filename = "xy" + to_string(j) + "_nointer.txt";
			myfile_nointer.open (filename,ios_base::app);					// append instead of overwrite
			myfile_nointer<<position_evol_inter;
		}
		
		myfile_t<<i*dt<<endl;
	}
	
	myfile_nointer.close();
	myfile_t.close();
	



return 0;

}