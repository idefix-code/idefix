// ********************************************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ********************************************************************************************************

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
        if(hallDiffusivityFunc) hallDiffusivityFunc(data, t, data.xHall);
        else IDEFIX_ERROR("No user-defined Hall diffusivity function has been enrolled");
    }

    switch (mySolver) {
    #if MHD == YES
        case TVDLF: TvdlfMHD(data, dir, this->gamma, this->C2Iso);
            break;
        case HLL:   HllMHD(data, dir, this->gamma, this->C2Iso, this->haveHall, this->xH);
            break;
        case HLLD:  HlldMHD(data, dir, this->gamma, this->C2Iso);
            break;
        case ROE:   RoeMHD(data, dir, this->gamma, this->C2Iso);
            break;
    #else
        case TVDLF: TvdlfHD(data, dir, this->gamma, this->C2Iso);
            break;
        case HLL:   HllHD(data, dir, this->gamma, this->C2Iso);
            break;
        case HLLC:  HllcHD(data, dir, this->gamma, this->C2Iso);
            break;
        case ROE:   RoeHD(data, dir, this->gamma, this->C2Iso);
            break;
    #endif
        default: // do nothing
            break;
    }

    idfx::popRegion();

}