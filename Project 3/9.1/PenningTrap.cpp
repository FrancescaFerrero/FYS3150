// Definitions for the functions in the LotteryMachine class

#include <vector>      // For vector
#include <string>      // For string
#include <stdlib.h>    // For rand (from C). For more powerful random number generation in C++ we should use <random>
#include <stdexcept>   // For runtime_error
#include "PenningTrap.hpp"
#include "Particle.hpp"

using namespace arma;
using namespace std;


// Constructor that takes a vector of Particle_collection as input
PenningTrap::PenningTrap(double B0_in, double V0_in, double d_in)
{
  B0_ = B0_in;
  V0_ = V0_in;
  d_ = d_in;
}


// Method that adds a particle to the Penning Trap by copying an input Particle
void PenningTrap::add_particle(Particle p_in)
{
  particle_collection.push_back(p_in);
}


// Calculate the external electric field at r=(x,y,z)
arma::vec PenningTrap::external_E_field(arma::vec r){

double V_d = 9.65;
vec E_field (3);


E_field(0) = V_d*r(0);  
E_field(1) = V_d*r(1); 
E_field(2) = -(2*V_d)*r(2);   

return E_field;  
}


// Calculate the external magnetic field at r=(x,y,z)
arma::vec PenningTrap::external_B_field(){

vec B_field = vec(3);

B_field(0) = 0;
B_field(1) = 0; 
B_field(2) = B0_;

return B_field; 
}


// Force on particle_i from particle_j
arma::vec PenningTrap::force_particle(int i, int j){

double k = 1.38935333e5;
double c = k*particle_collection[i].charge();
arma::vec dr = particle_collection[i].position()-particle_collection[j].position();
double dr_abs = pow(sqrt(pow(dr(0),2)+pow(dr(1),2)+pow(dr(2),2)),3);
vec F_particle (3);

for(int n=0; n<=2;n++){
F_particle(n) = particle_collection[j].charge()*(dr(n)/dr_abs);

}

return F_particle*c;  
}


// The total force on particle_i from the external fields
arma::vec PenningTrap::total_force_external(int i){

vec E_field = external_E_field(particle_collection[i].position());
vec B_field=external_B_field();
vec F_ext (3);

F_ext(0) = particle_collection[i].charge()*(E_field(0)+particle_collection[i].velocity()(1)*B_field(2)); 
F_ext(1) = particle_collection[i].charge()*(E_field(1)-particle_collection[i].velocity()(0)*B_field(2)); 
F_ext(2) = particle_collection[i].charge()*E_field(2);

return F_ext;
}


// The total force on particle_i from the other particles
arma::vec PenningTrap::total_force_particles(int i){

double k = 1.38935333e5;
double c = k*particle_collection[i].charge();
vec F_particles (3);

for(int n=0; n<=2;n++){
F_particles(n) = 0;
 for(int m=0;m<particle_collection.size();m++){
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

for(int n=0; n<=2;n++){
 F_tot(n) = total_force_particles(i)(n) + total_force_external(i)(n);
 
}

return F_tot;
}


// Evolve the system one time step (dt) using Runge-Kutta 4th order
void PenningTrap::evolve_RK4(double dt, int i, double t){                 
    vec k1_v,k1_r,k2_v,k2_r, k3_v,k3_r, k4_v,k4_r;
    double m =40.08;
    vec r_step (3);
    vec v_step(3);
     vec r_old (3);
      vec v_old (3);
      

    k1_r=dt*particle_collection[i].velocity();
    k1_v=dt*total_force(i,t)/m;
    
    r_old= particle_collection[i].position();
    v_old= particle_collection[i].velocity();
    
    particle_collection[i].position_=r_old +k1_r/2;
    particle_collection[i].velocity_ = v_old + k1_v/2;  
    
    k2_r= dt*particle_collection[i].velocity();
    k2_v= dt*total_force(i,t)/m;
    
    particle_collection[i].position_=r_old +k2_r/2;
    particle_collection[i].velocity_ = v_old + k2_v/2;  
    
    k3_r= dt*particle_collection[i].velocity();
    k3_v= dt*total_force(i,t)/m;
    
    particle_collection[i].position_=r_old +k3_r;
    particle_collection[i].velocity_ = v_old + k3_v;  
    
    k4_r= dt*particle_collection[i].velocity();
    k4_v= dt*total_force(i,t)/m;
    
    r_step= r_old + (1./6)*(k1_r + 2*k2_r +2*k3_r + k4_r);
    v_step= v_old + (1./6)*(k1_v + 2*k2_v +2*k3_v + k4_v);
    
    particle_collection[i].position_=r_step;
    particle_collection[i].velocity_ = v_step;

}


//Evolve the system one time step (dt) using Forward Euler
void PenningTrap::evolve_forward_Euler(double dt, int i){

vec position_step (3);
vec velocity_step (3);

position_step = particle_collection[i].position() + dt*particle_collection[i].velocity();
velocity_step= particle_collection[i].velocity() + dt*total_force(i,0)/40.08;
particle_collection[i].position_ =position_step;
particle_collection[i].velocity_ =velocity_step;

}



/*
void PenningTrap::evolve_forward_Euler(double dt, int i){

vec position_step (3);
vec velocity_step (3);
vec const_force (3, fill::ones);

position_step = particle_collection[i].position() + dt*particle_collection[i].velocity();
velocity_step= particle_collection[i].velocity() + dt*const_force/40.08;
particle_collection[i].position_ =position_step;
particle_collection[i].velocity_ =velocity_step;

}


// Evolve the system one time step (dt) using Runge-Kutta 4th order
void PenningTrap::evolve_RK4(double dt, int i){                 
    vec k1_v,k1_r,k2_v,k2_r, k3_v,k3_r, k4_v,k4_r;
    int t =0;
    double m =40.08;
    vec r_step (3);
    vec v_step(3);
    vec r_old (3);
    vec v_old (3);
    vec const_force (3, fill::ones);
      

    k1_r=dt*particle_collection[i].velocity();
    k1_v=dt*const_force/m;
    
    r_old= particle_collection[i].position();
    v_old= particle_collection[i].velocity();
    
    particle_collection[i].position_=r_old +k1_r/2;
    particle_collection[i].velocity_ = v_old + k1_v/2;  
    
    k2_r= dt*particle_collection[i].velocity();
    k2_v= dt*const_force/m;
    
    particle_collection[i].position_=r_old +k2_r/2;
    particle_collection[i].velocity_ = v_old + k2_v/2;  
    
    k3_r= dt*particle_collection[i].velocity();
    k3_v= dt*const_force/m;
    
    particle_collection[i].position_=r_old +k3_r;
    particle_collection[i].velocity_ = v_old + k3_v;  
    
    k4_r= dt*particle_collection[i].velocity();
    k4_v= dt*const_force/m;
    
    r_step= r_old + (1./6)*(k1_r + 2*k2_r +2*k3_r + k4_r);
    v_step= v_old + (1./6)*(k1_v + 2*k2_v +2*k3_v + k4_v);

}

*/
