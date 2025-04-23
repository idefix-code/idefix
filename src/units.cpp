// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include "units.hpp"
#include "idefix.hpp"
#include "input.hpp"


void idfx::Units::Init(Input &input) {
  if(input.CheckBlock("Units")) {
    this->_length = input.GetOrSet<real>("Units","length",0,1.0);
    this->_velocity = input.GetOrSet<real>("Units","velocity",0,1.0);
    this->_density = input.GetOrSet<real>("Units","density",0,1.0);
    this->ComputeUnits();
  }
}

void idfx::Units::SetDensity(const real in) {
  this->_density = in;
  this->ComputeUnits();
}

void idfx::Units::SetLength(const real in) {
  this->_length = in;
  this->ComputeUnits();
}

void idfx::Units::SetVelocity(const real in) {
  this->_velocity = in;
  this->ComputeUnits();
}

void idfx::Units::ShowConfig() {
  if(_isInitialized) {
    idfx::cout << "Units: [Length]      = " << this->_length << " cm" << std::endl;
    idfx::cout << "Units: [Velocity]    = " << this->_velocity << " cm/s" << std::endl;
    idfx::cout << "Units: [Density]     = " << this->_density << " g/cm3" << std::endl;
    idfx::cout << "Units: [Energy]     = " << this->_energy << " erg/cm3" << std::endl;
    idfx::cout << "Units: [Time]     = " << this->_time << " s" << std::endl;
    idfx::cout << "Units: [Temperature] = " << this->_Kelvin << " K" << std::endl;
    idfx::cout << "Units: [Mag. Field]  = " << this->_magField << " G" << std::endl;
  }
}

void idfx::Units::ComputeUnits() {
  this->_isInitialized = true;
  this->_magField = std::sqrt(4*M_PI*_density*_velocity*_velocity);
  this->_Kelvin = u * _velocity*_velocity/k_B;
  this->_energy = _density*_velocity*_velocity;
  this->_time = _length/_velocity;
}

