#include "../idefix.hpp"
#include "hydro.hpp"

#if MHD == YES
#include "solversMHD.hpp"
#else
#include "solversHD.hpp"
#endif

// Compute Riemann fluxes from states
void Hydro::CalcRiemannFlux(DataBlock & data, int dir) {

    idfx::pushRegion("Hydro::CalcRiemannFlux");
    
    switch (mySolver) {
    #if MHD == YES
        case TVDLF: TvdlfMHD(data, dir, this->gamma, this->C2Iso);
            break;
        case HLL:   HllMHD(data, dir, this->gamma, this->C2Iso);
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