// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef UNITS_HPP_
#define UNITS_HPP_

#include "idefix.hpp"

// Units class implemented in CGS
// Constants taken from Astropy
class Input;
namespace idfx {
class Units {
 public:
  void Init(Input &input);
  

  const real u{1.6605390666e-24};                   // Atomic mass unit (g)
  const real m_p{1.67262192369e-24};                // Proton mass unit (g)
  const real m_n{1.67492749804e-24};                // neutron mass unit (g)
  const real m_e{9.1093837015e-28};                 // electron mass unit (g)
  const real k_B{1.380649e-16};                     // Boltzmann constant (erg/K)
  const real sigma_sb{5.6703744191844314e-05};      // Stephan Boltzmann constant (g/(K^4 s^3))
  const real ar{7.5646e-15};                        // Radiation constant = 4*sigma_sb/c (g/(K^4 s^2 cm))
  const real c{29979245800.0};                      // Speed of light (cm/s)
  const real M_sun{1.988409870698051e+33};          // Solar mass (g)
  const real R_sun{69570000000.0};                  // Solar radius (cm)
  const real M_earth{5.972167867791379e+27};        // Earth mass (g)
  const real R_earth{637810000.0};                  // Earth radius (cm)
  const real G{6.674299999999999e-8};               // Gravitatonal constant  (cm3 / (g s2))
  const real h{6.62607015e-27};                     // Planck constant (erg.s)
  const real pc{3.08568e+18};                       // Parsec (cm)
  const real au{1.49597892e13};                       // Astronomical unit (cm)


  // User-defined units, non user-modifiable
  KOKKOS_INLINE_FUNCTION real GetLength() const {return _length;}                      // L (cm)  = L (code) * Units::length
  KOKKOS_INLINE_FUNCTION real GetVelocity() const {return _velocity;}                  // V(cm/s) = V(code) * Units::velocity
  KOKKOS_INLINE_FUNCTION real GetDensity() const {return _density;}                    // density (g/cm^3) = density(code)
                                                    //                     * Units::density

  // Deduced units from user-defined units
  KOKKOS_INLINE_FUNCTION  real GetKelvin() const {return _Kelvin;}                      // T(K) = P(code)/rho(code) * mu * Units::Kelvin
  KOKKOS_INLINE_FUNCTION  real GetMagField() const {return _magField;}                  // B(G) = B(code) * Units::MagField
  KOKKOS_INLINE_FUNCTION  real GetEnergy() const {return _energy;}                  // energy(erg/cm3) = energy(code) * Units::energy
  KOKKOS_INLINE_FUNCTION  real GetTime() const {return _time;}                  // time(s) = time(code) * Units::time

  KOKKOS_INLINE_FUNCTION  bool GetIsInitialized() const {return _isInitialized;}

  // code-style change of the units
  void SetLength(const real);
  void SetVelocity(const real);
  void SetDensity(const real);

  // Show the configuration
  void ShowConfig();

 private:
  bool _isInitialized{false};
  // read-write variables
  real _length{1.0};
  real _velocity{1.0};
  real _density{1.0};

  // Deduced units from user-defined units
  real _Kelvin{1.0};
  real _magField{1.0};
  real _time{1.0};
  real _energy{1.0};

  // Recompute deduced user-units
  void ComputeUnits();
};
} // namespace idfx
#endif // UNITS_HPP_

