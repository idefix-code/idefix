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
#include <cstdio>

#include "idefix.hpp"
#include "units.hpp"
#include "cooling.hpp"
#include "dataBlock.hpp"
#include "fluid.hpp"

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
  // Copy all required members to local device-capturable objects.
  // Capturing class members directly in a CUDA kernel would capture the host `this` pointer.
  auto Vc = this->Vc;
  auto delta_eng = this->delta_eng;
  auto temperature_data = this->temperature_data;
  auto Lambda_cool_data = this->Lambda_cool_data;
  auto eos = this->eos;
  real TcoolFloor = this->TcoolFloor;

  real kB = idfx::units.k_B;
  real mu = 0.609;
  real m_p = idfx::units.m_p;
  real XH = 0.71;

  // Code-to-physical (CGS) unit conversion factors
  real vel_unit = idfx::units.GetVelocity();
  real rho_unit = idfx::units.GetDensity();
  real len_unit = idfx::units.GetLength();

  int ibeg, iend, jbeg, jend, kbeg, kend;
  ibeg = this->data->beg[IDIR];
  iend = this->data->end[IDIR];
  jbeg = this->data->beg[JDIR];
  jend = this->data->end[JDIR];
  kbeg = this->data->beg[KDIR];
  kend = this->data->end[KDIR];

  idefix_for("RadCoolingLoop",kbeg, kend, jbeg, jend, ibeg, iend,
    KOKKOS_LAMBDA (int k, int j, int i) {
      int T_indx_lo = 0;
      int T_indx_hi = static_cast<int>(temperature_data.extent(0)) - 1;
      int T_indx_mid = 0;
      real temperature_max_data = temperature_data(T_indx_hi);
      real temperature_min_data = temperature_data(T_indx_lo);
      // ideal gas eos is used
      real temperature = Vc(PRS,k,j,i)/Vc(RHO,k,j,i)*(mu*m_p/kB)*pow(vel_unit,2);

      if (temperature<=TcoolFloor) {
        // not zero cooling but floor the temperature
        // delta_eng(k,j,i) = ZERO_F;
        real del_prs = -Vc(RHO,k,j,i)/(mu*m_p/kB)*(temperature-TcoolFloor)/pow(vel_unit,2);
        delta_eng(k,j,i) = eos.GetInternalEnergy(del_prs, Vc(RHO,k,j,i));
      } else if ((temperature < temperature_min_data) ||
                 (temperature > temperature_max_data)) {
        // tabulated data does not enclose the temperature value
        printf("RadCooling::TownsendIntegration Temperature out of range: T=%e, "
               "valid range=[%e, %e]\n",
               temperature, temperature_min_data, temperature_max_data);
        Kokkos::abort("RadCooling::TownsendIntegration Temperature out of range");
      } else {
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
        if (T_indx_lo != T_indx_hi) {
          // At loop exit, indices straddle the value as hi < lo.
          // Swap so lo <= hi before interpolation.
          T_indx_mid = T_indx_lo;
          T_indx_lo = T_indx_hi;
          T_indx_hi = T_indx_mid;
        }

        real temperature_lo = temperature_data(T_indx_lo);
        real temperature_hi = temperature_data(T_indx_hi);
        real Lambda_lo = Lambda_cool_data(T_indx_lo);
        real Lambda_hi = Lambda_cool_data(T_indx_hi);
        // T_ref = temperature_hi
        real Lambda_T;
        real alpha;
        if (temperature_hi == temperature_lo) {
          // Exact match in table: avoid 0/0 in slope estimate.
          alpha = ZERO_F;
          Lambda_T = Lambda_lo;
        } else {
          alpha = std::log(Lambda_hi/Lambda_lo)/std::log(temperature_hi/temperature_lo);
          Lambda_T = Lambda_lo * pow((temperature/temperature_lo), alpha);
        }
        // real gamma = eos.GetGamma(Vc(PRS,k,j,i),Vc(RHO,k,j,i));
        real eint = eos.GetInternalEnergy(Vc(PRS,k,j,i),
                                          Vc(RHO,k,j,i))
                                         *rho_unit*pow(vel_unit,2); // cgs
        real edot = pow(Vc(RHO,k,j,i)*rho_unit*XH/m_p,2)*Lambda_T; // cgs
        real t_cool = eint/edot; // cgs

        real Y, del_prs;
        if (alpha!=1.0) {
          Y = (1.0-pow(temperature_hi/temperature, alpha-1.0))/(1.0-alpha);
        } else {
          Y = std::log(temperature_hi/temperature);
        }
        real inv_Y_arg = Y + (temperature/temperature_hi)
                           * (Lambda_hi/Lambda_T)
                           * (dt/(t_cool/(len_unit/vel_unit)));
        real T_fin;
        if (alpha!=1.0) {
          T_fin = temperature_hi*pow((1.0 - (1.0-alpha)*inv_Y_arg), 1.0/(1.0-alpha));
        } else {
          T_fin = temperature_hi * std::exp(-inv_Y_arg);
        }
        if (T_fin>=TcoolFloor) {
          del_prs = -Vc(RHO,k,j,i)/(mu*m_p/kB)*(temperature-T_fin)/pow(vel_unit,2);
          delta_eng(k,j,i) = eos.GetInternalEnergy(del_prs, Vc(RHO,k,j,i));
        }
        } else {
          T_fin = TcoolFloor;
          del_prs = -Vc(RHO,k,j,i)/(mu*m_p/kB)*(temperature-T_fin)/pow(vel_unit,2);
          delta_eng(k,j,i) = eos.GetInternalEnergy(del_prs, Vc(RHO,k,j,i));
        }
      }
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
