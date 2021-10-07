// Definitions for the functions in the Particle class

#include "Particle.hpp"


// Constructor
Particle::Particle(arma::vec position_in, arma::vec velocity_in, double charge_in, double mass_in)
{
  charge_ = charge_in;
  mass_ = mass_in;
  position_ = position_in;
  velocity_ = velocity_in;
}


// Method that returns the charge
double Particle::charge()
{
  return charge_;
}

// Method that returns the mass
double Particle::mass()
{
  return mass_;
}

// Method that returns the vector position
arma::vec Particle::position()
{
  return position_;
}

// Method that returns the vector velocity
arma::vec Particle::velocity()
{
  return velocity_;
}



