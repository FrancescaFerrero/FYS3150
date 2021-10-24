// Definitions for the functions in the Penning trap class

#include <vector>      // For vector
#include <string>      // For string
#include <stdlib.h>    // For rand (from C). 
#include <stdexcept>   // For runtime_error
#include "PenningTrap.hpp"
#include "Particle.hpp"

using namespace arma;
using namespace std;


// Constructor 
PenningTrap::PenningTrap(double B0_in, double V0_in, double d_in, int interaction, double omega_v_in)
{
  B0_ = B0_in;
  V0_ = V0_in;
  d_ = d_in;
  interaction_ = interaction;
  omega_v_ = omega_v_in;
}

// Method that adds a particle to the Penning Trap by copying an input Particle
void PenningTrap::add_particle(Particle p_in){
  particle_collection.push_back(p_in);
}

// Method that counts how many particles are inside the trap
double PenningTrap::count_particles(mat r_step) {

    vec r(3);
    double counter = 0;
    
    for (int i=0; i<particle_collection.size(); i++){
        r=r_step.col(i);
       if (sqrt(pow(r(0),2)+pow(r(1),2)+pow(r(2),2))< d_){
           counter += 1;
        }
    }
    return counter;
}

// Calculate the external electric field at r=(x,y,z)
arma::vec PenningTrap::external_E_field(arma::vec r, double t){

double V_d = V0_/pow(d_,2);
vec E_field (3);
double f = 0.1;   //amplitude

if (sqrt(pow(r(0),2)+pow(r(1),2)+pow(r(2),2))> d_){
    E_field.zeros();
}
else {
    E_field(0) = V_d*r(0)*(1+f*cos(omega_v_*t));  
    E_field(1) = V_d*r(1)*(1+f*cos(omega_v_*t)); 
    E_field(2) = -(2*V_d)*r(2)*(1+f*cos(omega_v_*t));   
}

return E_field;  
}


// Calculate the external magnetic field at r=(x,y,z)
arma::vec PenningTrap::external_B_field(arma::vec r){

vec B_field = vec(3);

if (sqrt(pow(r(0),2)+pow(r(1),2)+pow(r(2),2))> d_){
    B_field.zeros();
}
else {

B_field(0) = 0;
B_field(1) = 0; 
B_field(2) = B0_;
}
return B_field; 
}


// Force on particle_i from particle_j
arma::vec PenningTrap::force_particle(int i, int j){

double k = 1.38935333e5;
double c = k*particle_collection[i].charge();
arma::vec dr = particle_collection[i].position()-particle_collection[j].position();
double dr_abs = pow(sqrt(pow(dr(0),2)+pow(dr(1),2)+pow(dr(2),2)),3);
vec F_particle (3);

for(int n=0; n<=2; n++){
F_particle(n) = particle_collection[j].charge()*(dr(n)/dr_abs);

}


return F_particle*c;  
}


// The total force on particle_i from the external fields
arma::vec PenningTrap::total_force_external(int i, double t){

vec E_field = external_E_field(particle_collection[i].position(),t);
vec B_field = external_B_field(particle_collection[i].position());
vec F_ext (3);

F_ext(0) = particle_collection[i].charge()*(E_field(0)+particle_collection[i].velocity()(1)*B_field(2)); 
F_ext(1) = particle_collection[i].charge()*(E_field(1)-particle_collection[i].velocity()(0)*B_field(2)); 
F_ext(2) = particle_collection[i].charge()*E_field(2);

return F_ext;
}


// The total force on particle_i from the other particles
arma::vec PenningTrap::total_force_particles(int i){


vec F_particles (3);

for(int n=0; n<=2; n++){
    F_particles(n) = 0;
 for(int m=0; m<particle_collection.size(); m++){
  if(m!=i){
 	F_particles(n) += force_particle(i, m)(n);
 	   }
 }
}

return F_particles;
}


// The total force on particle_i from both external fields and other particles
arma::vec PenningTrap::total_force(int i, double t){

vec F_tot (3);

for(int n=0; n<=2; n++){
    F_tot(n) = total_force_external(i,t)(n);
    if (interaction_ == 1){
        F_tot(n) += total_force_particles(i)(n);
    }
 
}

return F_tot;
}


// Evolve the system one time step (dt) using Runge-Kutta 4th order
void PenningTrap::evolve_RK4(double dt, int i , mat& r_step, mat& v_step, double t) {                 
    vec k1_v, k1_r, k2_v, k2_r, k3_v, k3_r, k4_v, k4_r;
    double m = 40.08;
    vec r_old (3);
    vec v_old (3);

    r_old = r_step.col(i);
    v_old = v_step.col(i);
    
    particle_collection[i].position_= r_old; 
    particle_collection[i].velocity_ = v_old; 
    k1_r = dt*particle_collection[i].velocity();
    k1_v = dt*total_force(i,t)/m;
    
    particle_collection[i].position_ = r_old +k1_r/2;
    particle_collection[i].velocity_ = v_old + k1_v/2;  
    
    k2_r = dt*particle_collection[i].velocity();
    k2_v = dt*total_force(i,t + 0.5*dt)/m;
    
    particle_collection[i].position_ = r_old +k2_r/2;
    particle_collection[i].velocity_ = v_old + k2_v/2;  
    
    k3_r = dt*particle_collection[i].velocity();
    k3_v = dt*total_force(i,t + 0.5*dt)/m;
    
    particle_collection[i].position_ = r_old +k3_r;
    particle_collection[i].velocity_ = v_old + k3_v;  
    
    k4_r = dt*particle_collection[i].velocity();
    k4_v = dt*total_force(i,t + dt)/m;
    
    r_step.col(i)= r_old + (1./6)*(k1_r + 2*k2_r + 2*k3_r + k4_r);
    v_step.col(i)= v_old + (1./6)*(k1_v + 2*k2_v + 2*k3_v + k4_v);
    
    particle_collection[i].position_= r_old;
    particle_collection[i].velocity_ = v_old;  
    
    

}


//Evolve the system one time step (dt) using Forward Euler
void PenningTrap::evolve_forward_Euler(double dt, int i, mat& r_step, mat& v_step){

vec r_old (3);
vec v_old (3);

r_old = r_step.col(i);
v_old = v_step.col(i);

r_step.col(i) = r_old + dt*v_old;
v_step.col(i) = v_old + dt*total_force(i,0)/40.08;
particle_collection[i].position_ = r_old;
particle_collection[i].velocity_ = v_old;

}



///////////////////////////////////////////////////
