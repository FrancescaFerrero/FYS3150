// The Particle class

#ifndef __Particle_hpp__  
#define __Particle_hpp__

#include <armadillo>


class Particle
{
  // Public stuff
  public:
    
    double charge_;
    double mass_;
    arma::vec position_;
    arma::vec velocity_;

    // Constructor
    Particle(arma::vec position_in, arma::vec velocity_in, double charge_in, double mass_in);

    // Method that returns the charge
    double charge();

    // Method that returns the mass
    double mass();
    
    // Method that returns the position
    arma::vec position();
    
    // Method that returns the velocity
    arma::vec velocity();

};

#endif
