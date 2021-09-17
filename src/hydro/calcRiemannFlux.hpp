// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2021 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef HYDRO_CALCRIEMANNFLUX_HPP_
#define HYDRO_CALCRIEMANNFLUX_HPP_

#include "hydro.hpp"
#include "dataBlock.hpp"

#if MHD == YES
#include "hlldMHD.hpp"
#include "hllMHD.hpp"
#include "roeMHD.hpp"
#include "tvdlfMHD.hpp"
#else
#include "hllcHD.hpp"
#include "hllHD.hpp"
#include "tvdlfHD.hpp"
#include "roeHD.hpp"
#endif

// Compute Riemann fluxes from states
template <int dir>
void Hydro::CalcRiemannFlux(const real t) {
  idfx::pushRegion("Hydro::CalcRiemannFlux");

  if(hallStatus.status == UserDefFunction && dir == IDIR) {
    if(hallDiffusivityFunc)
      hallDiffusivityFunc(*data, t, xHall);
    else
      IDEFIX_ERROR("No user-defined Hall diffusivity function has been enrolled");
  }

  if(haveIsoSoundSpeed == UserDefFunction && dir == IDIR) {
    if(isoSoundSpeedFunc)
      isoSoundSpeedFunc(*data, t, isoSoundSpeedArray);
    else
      IDEFIX_ERROR("No user-defined isothermal sound speed function has been enrolled");
  }

  switch (mySolver) {
#if MHD == YES
    case TVDLF:
      switch(emf.averaging) {
        case ElectroMotiveForce::arithmetic:
          TvdlfMHD<dir,ElectroMotiveForce::arithmetic>();
          break;
        case ElectroMotiveForce::uct0:
          TvdlfMHD<dir,ElectroMotiveForce::uct0>();
          break;
        case ElectroMotiveForce::uct_contact:
          TvdlfMHD<dir,ElectroMotiveForce::uct_contact>();
          break;
        case ElectroMotiveForce::uct_hll:
          TvdlfMHD<dir,ElectroMotiveForce::uct_hll>();
          break;
        case ElectroMotiveForce::uct_hlld:
          TvdlfMHD<dir,ElectroMotiveForce::uct_hlld>();
          break;
        default:
          IDEFIX_ERROR("EMF averaging scheme not implemented");
      }
      break;
    case HLL:
      switch(emf.averaging) {
        case ElectroMotiveForce::arithmetic:
          HllMHD<dir,ElectroMotiveForce::arithmetic>();
          break;
        case ElectroMotiveForce::uct0:
          HllMHD<dir,ElectroMotiveForce::uct0>();
          break;
        case ElectroMotiveForce::uct_contact:
          HllMHD<dir,ElectroMotiveForce::uct_contact>();
          break;
        case ElectroMotiveForce::uct_hll:
          HllMHD<dir,ElectroMotiveForce::uct_hll>();
          break;
        case ElectroMotiveForce::uct_hlld:
          HllMHD<dir,ElectroMotiveForce::uct_hlld>();
          break;
        default:
          IDEFIX_ERROR("EMF averaging scheme not implemented");
      }
      break;
    case HLLD:
      switch(emf.averaging) {
        case ElectroMotiveForce::arithmetic:
          HlldMHD<dir,ElectroMotiveForce::arithmetic>();
          break;
        case ElectroMotiveForce::uct0:
          HlldMHD<dir,ElectroMotiveForce::uct0>();
          break;
        case ElectroMotiveForce::uct_contact:
          HlldMHD<dir,ElectroMotiveForce::uct_contact>();
          break;
        case ElectroMotiveForce::uct_hll:
          HlldMHD<dir,ElectroMotiveForce::uct_hll>();
          break;
        case ElectroMotiveForce::uct_hlld:
          HlldMHD<dir,ElectroMotiveForce::uct_hlld>();
          break;
        default:
          IDEFIX_ERROR("EMF averaging scheme not implemented");
      }
      break;
    case ROE:
      switch(emf.averaging) {
        case ElectroMotiveForce::arithmetic:
          RoeMHD<dir,ElectroMotiveForce::arithmetic>();
          break;
        case ElectroMotiveForce::uct0:
          RoeMHD<dir,ElectroMotiveForce::uct0>();
          break;
        case ElectroMotiveForce::uct_contact:
          RoeMHD<dir,ElectroMotiveForce::uct_contact>();
          break;
        case ElectroMotiveForce::uct_hll:
          RoeMHD<dir,ElectroMotiveForce::uct_hll>();
          break;
        case ElectroMotiveForce::uct_hlld:
          RoeMHD<dir,ElectroMotiveForce::uct_hlld>();
          break;
        default:
          IDEFIX_ERROR("EMF averaging scheme not implemented");
      }
      break;
#else
    case TVDLF:
      TvdlfHD<dir>();
      break;
    case HLL:
      HllHD<dir>();
      break;
    case HLLC:
      HllcHD<dir>();
      break;
    case ROE:
      RoeHD<dir>();
      break;
#endif
    default: // do nothing
      break;
  }

  idfx::popRegion();
}
#endif // HYDRO_CALCRIEMANNFLUX_HPP_
