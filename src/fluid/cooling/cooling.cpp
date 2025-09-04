// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

// This source code is largely inspired from the viscous_flux of Pluto4.2
// ((c) P. Tzeferacos & A. Mignone)

// Implementation of monotonicity-preserving viscous flux following ZuHone et al.,
// ApJ

#include <string>

#include "idefix.hpp"
#include "cooling.hpp"
#include "dataBlock.hpp"
#include "fluid.hpp"

void RadCooling::InitUnits() {
  #if defined(UNIT_DENSITY) || defined(UNIT_MASS)
  #ifdef UNIT_DENSITY
  rho_unit = UNIT_DENSITY;
  #endif
  #ifdef UNIT_MASS
  mass_unit = UNIT_MASS;
  #endif
  #else
  #define UNIT_DENSITY    1.0
  rho_unit = 1.0;
  IDEFIX_WARNING("UNIT_DENSITY or UNIT_MASS "
                 "missing in definitions.hpp");
  #endif

  #if defined(UNIT_VELOCITY) || defined(UNIT_TIME)
  #ifdef UNIT_VELOCITY
  vel_unit = UNIT_VELOCITY;
  #endif
  #ifdef UNIT_TIME
  time_unit = UNIT_TIME;
  #endif
  #else
  #define UNIT_VELOCITY    1.0
  vel_unit = 1.0;
  IDEFIX_WARNING("UNIT_VELOCITY or UNIT_TIME "
                 "missing in definitions.hpp");
  #endif

  #ifdef UNIT_LENGTH
  len_unit = UNIT_LENGTH;
  #else
  #define UNIT_LENGTH    1.0
  len_unit = 1.0;
  IDEFIX_WARNING("UNIT_LENGTH "
                 "missing in definitions.hpp");
  #endif

  #ifndef UNIT_DENSITY
  rho_unit = mass_unit/pow(len_unit,3);
  #endif
  #ifndef UNIT_MASS
  mass_unit = rho_unit*pow(len_unit,3);
  #endif

  #ifndef UNIT_VELOCITY
  vel_unit = len_unit/time_unit;
  #endif
  #ifndef UNIT_TIME
  time_unit = len_unit/vel_unit;
  #endif
}

void RadCooling::InitArrays() {
   
  // Allocate and fill arrays when neede
  cooltable_data = LookupTable<1>(cooltable_filename,' ');
  temperature_data = cooltable_data.xinDev;
  Lambda_cool_data = cooltable_data.dataDev;

  delta_eng = IdefixArray3D<real>("delta_eng_cooling_source", data->np_tot[KDIR],
                                                              data->np_tot[JDIR],
                                                              data->np_tot[IDIR]);
}
void RadCooling::ShowConfig() {
  if(cooling_type.compare("tabulated")==0) {
    idfx::cout << "Optically thin radiative cooling: ENABLED"
               << std::endl;
    if(cooling_integration.compare("townsend") == 0) {
      idfx::cout << "Radiative cooling Integration: Townsend"
                 << std::endl;
    } else {
      IDEFIX_ERROR("Unknown radiative cooling integration mode");
    }
  } else {
    IDEFIX_ERROR("Unknown radiative cooling mode");
  }
  /*
    IDEFIX_WARNING("A sample warning "
                   "to be used with Idefix.");
  */
}

// This function computes the pressure drop due to optically thin radiative cooling
// using Townsend integration
void RadCooling::TownsendIntegration(real dt) {
  idfx::pushRegion("RadCooling::TownsendIntegration");
  IdefixArray1D<real> x1 = this->data->x[IDIR];
  IdefixArray1D<real> x2 = this->data->x[JDIR];
  IdefixArray1D<real> x3 = this->data->x[KDIR];
  IdefixArray3D<real> dV = this->data->dV;

  real kB = 1.380649e-16;
  real mu = 0.609;
  real mp = 1.6726e-24;
  real XH = 0.71;
  
  int ibeg, iend, jbeg, jend, kbeg, kend;
  ibeg = this->data->beg[IDIR];
  iend = this->data->end[IDIR];
  jbeg = this->data->beg[JDIR];
  jend = this->data->end[JDIR];
  kbeg = this->data->beg[KDIR];
  kend = this->data->end[KDIR];

  idefix_for("RadCoolingLoop",kbeg, kend, jbeg, jend, ibeg, iend,
    KOKKOS_LAMBDA (int k, int j, int i) {
      int T_indx_lo = 0, T_indx_hi=temperature_data.extent(0)-1;
      int T_indx_mid;
      // ideal gas eos is used
      real temperature = Vc(PRS,k,j,i)/Vc(RHO,k,j,i)*(mu*mp/kB)*pow(vel_unit,2);

      if ( (temperature<temperature_data(0)) || (temperature>temperature_data(temperature_data.extent(0)-1)) ) {
        char tmp[32];
          sprintf(tmp, "%e K", temperature);
          IDEFIX_ERROR("Temperature out of range: T="+std::string(tmp)); 
      }
      
      while (T_indx_lo<=T_indx_hi) {
        T_indx_mid = (T_indx_lo + T_indx_hi)/2;
        if (temperature < temperature_data(T_indx_mid)) {
          T_indx_hi = T_indx_mid-1;
        } else if (temperature > temperature_data(T_indx_mid)) {
          T_indx_lo = T_indx_mid+1;
        } else {
          T_indx_lo = T_indx_mid;
          T_indx_hi = T_indx_mid;
          break;
        }
      }
      if (T_indx_lo!=T_indx_hi) {
        // swap
        T_indx_mid = T_indx_lo;
        T_indx_lo = T_indx_hi;
        T_indx_hi = T_indx_mid;
      }

      real temperature_lo = temperature_data(T_indx_lo); 
      real temperature_hi = temperature_data(T_indx_hi);
      real Lambda_lo, Lambda_hi = Lambda_cool_data(T_indx_lo), Lambda_cool_data(T_indx_hi);
      // T_ref = temperature_hi
      real alpha = std::log(Lambda_hi/Lambda_lo)/std::log(temperature_hi/temperature_lo);
      real Lambda_T = Lambda_lo * pow((temperature/temperature_lo), alpha);
      // real gamma = eos.GetGamma(Vc(PRS,k,j,i),Vc(RHO,k,j,i));
      real t_cool = eos.GetInternalEnergy(Vc(PRS,k,j,i), Vc(RHO,k,j,i))*rho_unit*pow(vel_unit,2)/(pow(Vc(RHO,k,j,i)*rho_unit*XH/mp,2)*Lambda_T); // cgs

      real Y, del_prs;
      if (alpha!=1.0) {
        Y = (1.0-pow(temperature_hi/temperature, alpha-1))/(1-alpha);
      } else {
        Y = std::log(temperature_hi/temperature);
      }
      real inv_Y_arg = Y + (temperature/temperature_hi) * (Lambda_hi/Lambda_T) * dt/(t_cool/(len_unit/vel_unit));
      real T_fin;
      if (alpha!=1.0) {
        T_fin = pow(temperature_hi*(1.0 - (1-alpha)*inv_Y_arg), 1.0/(1-alpha));
      } else {
        T_fin = temperature_hi * std::exp(-inv_Y_arg);
      }
      // ideal gas eos is used
      del_prs = -Vc(RHO,k,j,i)/(mu*mp/kB)*std::min(temperature-T_fin, temperature-TcoolFloor)/pow(vel_unit,2);
      delta_eng(k,j,i) = eos.GetInternalEnergy(del_prs, Vc(RHO,k,j,i));
    });
  idfx::popRegion();
}

// This function computes the energy loss in the conservative variable due
// to optically thin radiative cooling
void RadCooling::CalculateCoolingSource(real dt) {
  idfx::pushRegion("RadCooling::CalculateCoolingSource");
  // calculate the change in internal thermal energy
  rad_cooling_int_func(dt);
  idfx::popRegion();
}
