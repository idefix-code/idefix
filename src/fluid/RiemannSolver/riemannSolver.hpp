// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef FLUID_RIEMANNSOLVER_RIEMANNSOLVER_HPP_
#define FLUID_RIEMANNSOLVER_RIEMANNSOLVER_HPP_

#include <string>
#include <memory>

#include "fluid.hpp"
#include "input.hpp"

// Forward declaration
template<typename Phys>
class ShockFlattening;

#include "extrapolateToFaces.hpp"

template <typename Phys>
class RiemannSolver {
 public:
  // Riemann Solver type

  enum Solver {TVDLF_MHD, HLL_MHD, HLLD_MHD, ROE_MHD, TVDLF, HLL, HLLC, ROE, HLL_DUST};

  RiemannSolver(Input &input, Fluid<Phys>* hydro);

  template <int> void CalcFlux(IdefixArray4D<real> &);

  Solver GetSolver() {
    return(mySolver);
  }

  void ShowConfig();

  // Riemann Solvers
  template<const int>
    void HlldMHD(IdefixArray4D<real> &);
  template<const int>
    void HllMHD(IdefixArray4D<real> &);
  template<const int>
    void RoeMHD(IdefixArray4D<real> &);
  template<const int>
    void TvdlfMHD(IdefixArray4D<real> &);

  template<const int>
    void HllcHD(IdefixArray4D<real> &);
  template<const int>
    void HllHD(IdefixArray4D<real> &);
  template<const int>
    void RoeHD(IdefixArray4D<real> &);
  template<const int>
    void TvdlfHD(IdefixArray4D<real> &);

  template<const int>
    void HllDust(IdefixArray4D<real> &);
  // Get the right slope limiter
  template<int dir>
  ExtrapolateToFaces<Phys, dir>* GetExtrapolator();

  std::unique_ptr<ShockFlattening<Phys>> shockFlattening;

 private:
  template <typename P, int dir, PLMLimiter L, int O>
  friend class ExtrapolateToFaces;

  IdefixArray4D<real> Vc;
  IdefixArray4D<real> Vs;
  IdefixArray4D<real> Flux;
  IdefixArray3D<real> cMax;
  Fluid<Phys>* hydro;
  DataBlock *data;

  Solver mySolver;

  // Because each direction is a different template, we can't use
  std::unique_ptr<ExtrapolateToFaces<Phys,IDIR>> slopeLimIDIR;
  std::unique_ptr<ExtrapolateToFaces<Phys,JDIR>> slopeLimJDIR;
  std::unique_ptr<ExtrapolateToFaces<Phys,KDIR>> slopeLimKDIR;

  bool haveShockFlattening;
};

#include "shockFlattening.hpp"

template <typename Phys>
RiemannSolver<Phys>::RiemannSolver(Input &input, Fluid<Phys>* hydro) : Vc{hydro->Vc},
                                      Vs{hydro->Vs},
                                      Flux{hydro->FluxRiemann},
                                      cMax{hydro->cMax},
                                      hydro{hydro},
                                      data{hydro->data}
                                      {
  // read Solver from input file
  if(!Phys::dust) {
    std::string solverString = input.Get<std::string>(std::string(Phys::prefix),"solver",0);
    if (solverString.compare("tvdlf") == 0) {
      if constexpr(Phys::mhd) {
        mySolver = TVDLF_MHD;
      } else {
        mySolver = TVDLF;
      }
    } else if (solverString.compare("hll") == 0) {
      if constexpr(Phys::mhd) {
        mySolver = HLL_MHD;
      } else {
        mySolver = HLL;
      }
    } else if (solverString.compare("hlld") == 0) {
      if constexpr(Phys::mhd) {
      mySolver = HLLD_MHD;
      } else {
        IDEFIX_ERROR("hlld Riemann solver requires a MHD fluid");
      }
    } else if (solverString.compare("hllc") == 0) {
      if constexpr(Phys::mhd) {
        IDEFIX_ERROR("hllc Riemann solver requires a HD fluid");
      } else {
          mySolver = HLLC;
      }
    } else if (solverString.compare("roe") == 0) {
      if constexpr(Phys::mhd) {
        mySolver = ROE_MHD;
      } else {
        mySolver = ROE;
      }
    } else {
      std::stringstream msg;
      if constexpr(Phys::mhd) {
        msg << "Unknown MHD solver type " << solverString;
      } else {
        msg << "Unknown HD solver type " << solverString;
      }
      IDEFIX_ERROR(msg);
    }
    // Check if Hall is enabled
    if(input.CheckEntry(std::string(Phys::prefix),"hall")>=0) {
        // Check consistency
        if(mySolver != HLL_MHD )
          IDEFIX_ERROR("Hall effect is only compatible with HLL Riemann solver.");
    }
  } else {
    // We're dealing with dust grains
    mySolver = HLL_DUST;
  }



  // Shock flattening
  this->haveShockFlattening = input.CheckEntry(std::string(Phys::prefix),"shockFlattening")>=0;
  // Init shock flattening
  if(haveShockFlattening) {
    this->shockFlattening = std::make_unique<ShockFlattening<Phys>>(
                              hydro,input.Get<real>(std::string(Phys::prefix),"shockFlattening",0));
  }

  // init slope limiters
  slopeLimIDIR = std::make_unique<ExtrapolateToFaces<Phys,IDIR>>(this);
  #if DIMENSIONS >= 2
  slopeLimJDIR = std::make_unique<ExtrapolateToFaces<Phys,JDIR>>(this);
  #endif
  #if DIMENSIONS == 3
  slopeLimKDIR = std::make_unique<ExtrapolateToFaces<Phys,KDIR>>(this);
  #endif
}

template <typename Phys>
void RiemannSolver<Phys>::ShowConfig() {
  idfx::cout << "RiemannSolver: ";
  switch(mySolver) {
    case TVDLF:
      idfx::cout << "tvdlf (HD)." << std::endl;
      break;
    case TVDLF_MHD:
      idfx::cout << "tvdlf (MHD)." << std::endl;
      break;
    case HLL:
      idfx::cout << "hll (HD)." << std::endl;
      break;
    case HLL_MHD:
      idfx::cout << "hll (MHD)." << std::endl;
      break;
    case HLLD_MHD:
        idfx::cout << "hlld (MHD)." << std::endl;
        break;
    case HLLC:
      idfx::cout << "hllc (HD)." << std::endl;
      break;
    case ROE:
      idfx::cout << "roe (HD)." << std::endl;
      break;
    case ROE_MHD:
      idfx::cout << "roe (MHD)." << std::endl;
      break;
    case HLL_DUST:
      idfx::cout << "HLL (Dust)." << std::endl;
      break;
    default:
      IDEFIX_ERROR("Unknown Riemann solver");
  }

  if(haveShockFlattening) {
    idfx::cout << Phys::prefix << ": Shock Flattening ENABLED." << std::endl;
  }
}

template <typename Phys>
template<const int dir>
ExtrapolateToFaces<Phys, dir>* RiemannSolver<Phys>::GetExtrapolator() {
  if constexpr(dir==IDIR) {
    return(this->slopeLimIDIR.get());
  } else if constexpr(dir==JDIR) {
    return(this->slopeLimJDIR.get());
  } else if constexpr(dir==KDIR) {
     return(this->slopeLimKDIR.get());
  }
}


#include "calcFlux.hpp"

#endif //FLUID_RIEMANNSOLVER_RIEMANNSOLVER_HPP_
