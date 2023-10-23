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
#include "bragViscosity.hpp"
#include "dataBlock.hpp"
#include "fluid.hpp"



void BragViscosity::InitArrays() {
  // Allocate and fill arrays when needed
  #if GEOMETRY != CARTESIAN
    one_dmu = IdefixArray1D<real>("BragViscosity_1dmu", data->np_tot[JDIR]);
    IdefixArray1D<real> dmu = one_dmu;
    IdefixArray1D<real> th = data->x[JDIR];
    idefix_for("ViscousInitGeometry",1,data->np_tot[JDIR],
      KOKKOS_LAMBDA(int j) {
        real scrch =  FABS((1.0-cos(th(j)))*(sin(th(j)) >= 0.0 ? 1.0:-1.0)
                     -(1.0-cos(th(j-1))) * (sin(th(j-1)) > 0.0 ? 1.0:-1.0));
        dmu(j) = 1.0/scrch;
      });
  #endif
  bragViscSrc = IdefixArray4D<real>("BragViscosity_source", COMPONENTS, data->np_tot[KDIR],
                                                                data->np_tot[JDIR],
                                                                data->np_tot[IDIR]);
}
void BragViscosity::ShowConfig() {
  if(status.status==Constant) {
    idfx::cout << "Braginskii Viscosity: ENABLED with constant braginskii viscosity etaBrag="
                    << this->etaBrag << " ."<< std::endl;
  } else if (status.status==UserDefFunction) {
    idfx::cout << "Braginskii Viscosity: ENABLED with user-defined braginskii viscosity function."
                   << std::endl;
    if(!bragViscousDiffusivityFunc) {
      IDEFIX_ERROR("No braginskii viscosity function has been enrolled");
    }
  } else {
    IDEFIX_ERROR("Unknown braginskii viscosity mode");
  }
  if(this->status.isExplicit) {
    idfx::cout << "Braginskii Viscosity: uses an explicit time integration." << std::endl;
  } else if(this->status.isRKL) {
    idfx::cout << "Braginskii Viscosity: uses a Runge-Kutta-Legendre time integration."
                << std::endl;
  } else {
    IDEFIX_ERROR("Unknown time integrator for braginskii viscosity.");
  }
  if(haveSlopeLimiter) {
    idfx::cout << "Braginskii Viscosity: uses a slope limiter." << std::endl;
  }

  #if GEOMETRY == CYLINDRICAL || GEOMETRY == POLAR
    IDEFIX_WARNING("Braginskii viscosity in cylindrical and polar coordinates "
                   "has been implemented but never tested so far.");
  #endif
}

void BragViscosity::EnrollBragViscousDiffusivity(DiffusivityFunc myFunc) {
  if(this->status.status != UserDefFunction) {
    IDEFIX_WARNING("Braginskii viscous diffusivity enrollment requires Hydro/BragViscosity "
                 "to be set to userdef in .ini file");
  }
  this->bragViscousDiffusivityFunc = myFunc;
}

// This function computes the viscous flux and stores it in hydro->fluxRiemann
// (this avoids an extra array)
// Associated source terms, present in non-cartesian geometry are also computed
// and stored in this->viscSrc for later use (in calcRhs).
void BragViscosity::AddBragViscousFlux(int dir, const real t, const IdefixArray4D<real> &Flux) {
  idfx::pushRegion("BragViscosity::AddBragViscousFlux");
  switch(limiter) {
    case PLMLimiter::VanLeer:
      this->AddBragViscousFluxLim<PLMLimiter::VanLeer>(dir,t,Flux);
      break;
    case PLMLimiter::McLim:
      this->AddBragViscousFluxLim<PLMLimiter::McLim>(dir,t,Flux);
      break;
    case PLMLimiter::MinMod:
      this->AddBragViscousFluxLim<PLMLimiter::MinMod>(dir,t,Flux);
      break;
    default:
      IDEFIX_ERROR("The slope limiter for the Braginskii viscosity is not defined.");
      break;
  }
  idfx::popRegion();
}
