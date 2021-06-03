// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2021 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include "hydro.hpp"
#include "dataBlock.hpp"

#if MHD == YES
#include "solversMHD.hpp"
#else
#include "solversHD.hpp"
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
      TvdlfMHD<dir>();
      break;
    case HLL:
      HllMHD<dir>();
      break;
    case HLLD:
      HlldMHD<dir>();
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
    case ROE:
      RoeHD<dir>();
      break;
#endif
    default: // do nothing
      break;
  }

  idfx::popRegion();
}
