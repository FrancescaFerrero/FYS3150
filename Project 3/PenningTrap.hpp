// The LotteryMachine class

#ifndef __PenningTrap_hpp__  
#define __PenningTrap_hpp__

#include <string>
#include <armadillo>
#include <vector>
#include "Particle.hpp"  // Some of the declarations below need the Particle type


class PenningTrap
{
  // Public stuff
  public:

    double B0_;
    double V0_;
    double d_;

    // The PenningTrap holds a collection of Particles
    std::vector<Particle> particle_collection;

    // Constructor that creates an empty PenningTrap
    PenningTrap();

    // Constructor
    PenningTrap(double B0_in, double V0_in, double d_in);

    // Method that returns the magnetic field strength B0
    //double B0();
    
    // Method that returns the applied potential V0
    //double V0(); 

    // Method that returns the characteristic dimension d
    //double d();    

  // Add a particle to the trap
  void add_particle(Particle p_in);

  // External electric field at point r=(x,y,z)
  arma::vec external_E_field(arma::vec r);  

  // External magnetic field at point r=(x,y,z)
  arma::vec external_B_field(); 
  
  // Force on particle_i from particle_j
  arma::vec force_particle(int i, int j);

  // The total force on particle_i from the external fields
  arma::vec total_force_external(int i);

  // The total force on particle_i from the other particles
  arma::vec total_force_particles(int i);

  // The total force on particle_i from both external fields and other particles
  arma::vec total_force(int i, double t);

  // Evolve the system one time step (dt) using Runge-Kutta 4th order
  void evolve_RK4(double dt, int i);

  // Evolve the system one time step (dt) using Forward Euler
  void evolve_forward_Euler(double dt, int i);

};


#endif



