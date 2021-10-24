#include "PenningTrap.hpp"
#include "Particle.hpp"
#include <algorithm>
#include <ctime>
#include <cmath>
#include <iomanip>

using namespace std;
using namespace arma;

int main (int argc, char* argv[]){ 
   
double n_step = atoi(argv[1]); 
double t_tot = 100.;  // total time
double dt = t_tot/n_step; // time step
int n_particles = 1;  // number of particles
int j = 0;

vec position_0 (3);
vec velocity_0 (3);
mat r_step (3, n_particles);
mat v_step (3, n_particles);
cube R (3, n_step, n_particles, fill::zeros);
cube V (3, n_step, n_particles, fill::zeros);

mat r_step_eu (3, n_particles);
mat v_step_eu (3, n_particles);
cube R_eu (3, n_step, n_particles, fill::zeros);
cube V_eu (3, n_step, n_particles, fill::zeros);

mat R_an (3, n_step);

PenningTrap my_trap_nointer = PenningTrap(9.65e1, 9.65e8, 1.0e4, 0); // Set a penning trap without interactions between particles
Particle my_particle = Particle(position_0, velocity_0, 1, 40.08); // Set a particle of Ca++


// Open file to save time values
ofstream myfile_t;
myfile_t.open ("t_"+std::to_string((int)n_step)+".txt");

if (!myfile_t) { // file couldn't be opened
   cerr << "Error: file could not be opened" << endl;
   exit(1);
}

// Set initial position and velocity
position_0 = {10, 0, 5};	
velocity_0 = {0, 8, 0};

// Assign initial position and velocity
my_particle.position_ = position_0;
my_particle.velocity_ = velocity_0;

// Add the particle to the trap
my_trap_nointer.add_particle(my_particle);


// no interactions
	for (int i=0; i<n_step; i++) {
			if (i==0){
				// Runge-Kutta
				R.slice(j).col(i) = my_trap_nointer.particle_collection[j].position();
				V.slice(j).col(i) = my_trap_nointer.particle_collection[j].velocity();
				r_step.col(j) = my_trap_nointer.particle_collection[j].position();
				v_step.col(j) = my_trap_nointer.particle_collection[j].velocity();
				
				// Forward Euler
				R_eu.slice(j).col(i) = my_trap_nointer.particle_collection[j].position();
				V_eu.slice(j).col(i) = my_trap_nointer.particle_collection[j].velocity();
				r_step_eu.col(j) = my_trap_nointer.particle_collection[j].position();
				v_step_eu.col(j) = my_trap_nointer.particle_collection[j].velocity();
			}
			else {
				// Runge-Kutta
				my_trap_nointer.evolve_RK4(dt, j, r_step, v_step);
				R.slice(j).col(i) = r_step.col(j);
				V.slice(j).col(i) = v_step.col(j);
				
				// Forward Euler
				my_trap_nointer.evolve_forward_Euler(dt, j, r_step_eu, v_step_eu);
				R_eu.slice(j).col(i) = r_step_eu.col(j);
				V_eu.slice(j).col(i) = v_step_eu.col(j);
			}
		myfile_t<<i*dt<<endl;
	}
	
	myfile_t.close();

double omega_0 = 9.65e1/40.08;
double omega_z = sqrt(2*9.65/40.08); 

double omega_p = (omega_0+sqrt(pow(omega_0,2)-2*pow(omega_z,2)))/2;
double omega_m = (omega_0-sqrt(pow(omega_0,2)-2*pow(omega_z,2)))/2;  

double A_p =  (velocity_0(1)+omega_m*position_0(0))/(omega_m-omega_p);
double A_m = -(velocity_0(1)+omega_p*position_0(0))/(omega_m-omega_p);	
	
// Calculate analytic r values
for (int i=0; i<n_step; i++){
	R_an(0,i) = A_p*cos(omega_p*dt*i) + A_m*cos(omega_m*dt*i);
	R_an(1,i) = -A_p*sin(omega_p*dt*i) - A_m*sin(omega_m*dt*i);
	R_an(2,i) = position_0(2)*cos(omega_z*dt*i);
}
	
	
// Calculate absoloute and relative error - RK and FE
vec abs_err_RK(n_step);
vec abs_err_FE(n_step);
vec rel_err_RK(n_step);
vec rel_err_FE(n_step);
  
for (int i=0; i<n_step; i++){
abs_err_RK(i) = norm(R_an.col(i)-R.slice(0).col(i),2);
abs_err_FE(i) = norm(R_an.col(i)-R_eu.slice(0).col(i),2);

rel_err_RK(i) = abs_err_RK(i)/norm(R_an.col(i)); 
rel_err_FE(i) = abs_err_FE(i)/norm(R_an.col(i));

}
	ofstream myfile6;
	myfile6.open ("abs_err_RK."+std::to_string((int)n_step)+".txt");
	myfile6<< fixed << setprecision(16) << abs_err_RK;
	myfile6.close();

	ofstream myfile7;
	myfile7.open ("abs_err_FE."+std::to_string((int)n_step)+".txt");
	myfile7<< fixed << setprecision(16) << abs_err_FE;
	myfile7.close();

	ofstream myfile8;
	myfile8.open ("rel_err_RK."+std::to_string((int)n_step)+".txt");
	myfile8<< fixed << setprecision(16) << rel_err_RK;
	myfile8.close();

	ofstream myfile9;
	myfile9.open ("rel_err_FE."+std::to_string((int)n_step)+".txt");
	myfile9<< fixed << setprecision(16) << rel_err_FE;
	myfile9.close();

return 0;
}

