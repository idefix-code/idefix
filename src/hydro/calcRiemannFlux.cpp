// ***********************************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************************

#include "../idefix.hpp"
#include "hydro.hpp"

#if MHD == YES
#include "solversMHD.hpp"
#else
#include "solversHD.hpp"
#endif

// Compute Riemann fluxes from states
void Hydro::CalcRiemannFlux(DataBlock & data, int dir, const real t) {
  idfx::pushRegion("Hydro::CalcRiemannFlux");

  if(haveHall == UserDefFunction && dir == IDIR) {
    if(hallDiffusivityFunc)
      hallDiffusivityFunc(data, t, data.xHall);
    else
      IDEFIX_ERROR("No user-defined Hall diffusivity function has been enrolled");
  }

  switch (mySolver) {
#if MHD == YES
    case TVDLF:
      switch(dir) {
        case IDIR:
          TvdlfMHD<IDIR,ARG_EXPAND(MX1,MX2,MX3),
                   ARG_EXPAND(BX1,BX2,BX3)>(data, this->gamma, this->C2Iso);
          break;
        case JDIR:
          TvdlfMHD<JDIR,ARG_EXPAND(MX2,MX1,MX3),
                   ARG_EXPAND(BX2,BX1,BX3)>(data, this->gamma, this->C2Iso);
          break;
        case KDIR:
          TvdlfMHD<KDIR,ARG_EXPAND(MX3,MX1,MX2),
                   ARG_EXPAND(BX3,BX1,BX2)>(data, this->gamma, this->C2Iso);
          break;
      }
      break;
    case HLL:
      switch(dir) {
        case IDIR:
          HllMHD<IDIR,ARG_EXPAND(MX1,MX2,MX3),
                 ARG_EXPAND(BX1,BX2,BX3)>(data, this->gamma, this->C2Iso,
                                               this->haveHall, this->xH);
          break;
        case JDIR:
          HllMHD<JDIR,ARG_EXPAND(MX2,MX1,MX3),
                 ARG_EXPAND(BX2,BX1,BX3)>(data, this->gamma, this->C2Iso,
                                               this->haveHall, this->xH);
          break;
        case KDIR:
          HllMHD<KDIR,ARG_EXPAND(MX3,MX1,MX2),
                 ARG_EXPAND(BX3,BX1,BX2)>(data, this->gamma, this->C2Iso,
                                               this->haveHall, this->xH);
          break;
      }
      break;
    case HLLD:
      switch(dir) {
        case IDIR:
          HlldMHD<IDIR,ARG_EXPAND(MX1,MX2,MX3),
                  ARG_EXPAND(BX1,BX2,BX3)>(data, this->gamma, this->C2Iso);
          break;
        case JDIR:
          HlldMHD<JDIR,ARG_EXPAND(MX2,MX1,MX3),
                  ARG_EXPAND(BX2,BX1,BX3)>(data, this->gamma, this->C2Iso);
          break;
        case KDIR:
          HlldMHD<KDIR,ARG_EXPAND(MX3,MX1,MX2),
                  ARG_EXPAND(BX3,BX1,BX2)>(data, this->gamma, this->C2Iso);
          break;
      }
      break;
    case ROE:
      switch(dir) {
        case IDIR:
          RoeMHD<IDIR,ARG_EXPAND(MX1,MX2,MX3),
                 ARG_EXPAND(BX1,BX2,BX3)>(data, this->gamma, this->C2Iso);
          break;
        case JDIR:
          RoeMHD<JDIR,ARG_EXPAND(MX2,MX1,MX3),
                 ARG_EXPAND(BX2,BX1,BX3)>(data, this->gamma, this->C2Iso);
          break;
        case KDIR:
          RoeMHD<KDIR,ARG_EXPAND(MX3,MX1,MX2),
                 ARG_EXPAND(BX3,BX1,BX2)>(data, this->gamma, this->C2Iso);
          break;
      }
      break;
#else
    case TVDLF: 
      switch(dir) {
        case IDIR:
          TvdlfHD<IDIR,MX1,MX2,MX3>(data, this->gamma, this->C2Iso);
          break;
        case JDIR:
          TvdlfHD<JDIR,MX2,MX1,MX3>(data, this->gamma, this->C2Iso);
          break;
        case KDIR:
          TvdlfHD<KDIR,MX3,MX1,MX2>(data, this->gamma, this->C2Iso);
          break;
      }
      break;
    case HLL:
      switch(dir) {
        case IDIR:
          HllHD<IDIR,MX1,MX2,MX3>(data, this->gamma, this->C2Iso);
          break;
        case JDIR:
          HllHD<JDIR,MX2,MX1,MX3>(data, this->gamma, this->C2Iso);
          break;
        case KDIR:
          HllHD<KDIR,MX3,MX1,MX2>(data, this->gamma, this->C2Iso);
          break;
      }
      break;
    case HLLC:
      switch(dir) {
        case IDIR:
          HllcHD<IDIR,MX1,MX2,MX3>(data, this->gamma, this->C2Iso);
          break;
        case JDIR:
          HllcHD<JDIR,MX2,MX1,MX3>(data, this->gamma, this->C2Iso);
          break;
        case KDIR:
          HllcHD<KDIR,MX3,MX1,MX2>(data, this->gamma, this->C2Iso);
          break;
      }
      break;
    case ROE:
      switch(dir) {
        case IDIR:
          RoeHD<IDIR,MX1,MX2,MX3>(data, this->gamma, this->C2Iso);
          break;
        case JDIR:
          RoeHD<JDIR,MX2,MX1,MX3>(data, this->gamma, this->C2Iso);
          break;
        case KDIR:
          RoeHD<KDIR,MX3,MX1,MX2>(data, this->gamma, this->C2Iso);
          break;
      }
      break;
#endif
    default: // do nothing
      break;
  }

  idfx::popRegion();
}
