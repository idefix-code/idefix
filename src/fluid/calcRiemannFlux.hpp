// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef FLUID_CALCRIEMANNFLUX_HPP_
#define FLUID_CALCRIEMANNFLUX_HPP_

#include "fluid.hpp"
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
template <typename Phys>
template <int dir>
void Fluid<Phys>::CalcRiemannFlux(const real t) {
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
      TvdlfMHD<dir>();
      break;
    case HLL:
      HllMHD<dir>();
      break;
    case HLLD:
      HlldMHD<dir>();
      break;
    case ROE:
      RoeMHD<dir>();
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
#endif // FLUID_CALCRIEMANNFLUX_HPP_
