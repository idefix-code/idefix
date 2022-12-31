// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2022 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef FLUID_RIEMANNSOLVER_RIEMANNSOLVER_HPP_
#define FLUID_RIEMANNSOLVER_RIEMANNSOLVER_HPP_

#include "fluid.hpp"
#include "input.hpp"
#include "shockFlattening.hpp"
#include "../physics.hpp"
using Hydro = Fluid<Physics>;



template <typename Phys>
class RiemannSolver {
 public:
  // Riemann Solver type
  #if MHD == YES
  enum Solver {TVDLF=1, HLL, HLLD, ROE};
  #else
  enum Solver {TVDLF=1, HLL, HLLC, ROE};
  #endif

  RiemannSolver(Input &input, Fluid<Phys>* hydro);
  
  template <int> void CalcFlux(IdefixArray4D<real> &);

  Solver GetSolver() {
    return(mySolver);
  };

  void ShowConfig();

 // Riemann Solvers
#if MHD == YES
  template<const int>
    void HlldMHD(IdefixArray4D<real> &);
  template<const int>
    void HllMHD(IdefixArray4D<real> &);
  template<const int>
    void RoeMHD(IdefixArray4D<real> &);
  template<const int>
    void TvdlfMHD(IdefixArray4D<real> &);
#else
  template<const int>
    void HllcHD(IdefixArray4D<real> &);
  template<const int>
    void HllHD(IdefixArray4D<real> &);
  template<const int>
    void RoeHD(IdefixArray4D<real> &);
  template<const int>
    void TvdlfHD(IdefixArray4D<real> &);
#endif

 private:
  IdefixArray4D<real> Vc;
  IdefixArray4D<real> Vs;
  IdefixArray4D<real> Flux;
  IdefixArray3D<real> cMax;
  Fluid<Phys>* hydro;
  DataBlock *data;

  Solver mySolver;

  ShockFlattening shockFlattening;
  bool haveShockFlattening;
  
};

template <typename Phys>
RiemannSolver<Phys>::RiemannSolver(Input &input, Fluid<Phys>* hydro) : Vc{hydro->Vc},
                                      Vs{hydro->Vs}, 
                                      Flux{hydro->FluxRiemann}, 
                                      cMax{hydro->cMax},
                                      hydro{hydro},
                                      data{hydro->data}
                                      {

  // read Solver from input file
  std::string solverString = input.Get<std::string>(std::string(Phys::prefix),"solver",0);

  if (solverString.compare("tvdlf") == 0) {
    mySolver = TVDLF;
  } else if (solverString.compare("hll") == 0) {
    mySolver = HLL;
#if MHD == YES
  } else if (solverString.compare("hlld") == 0) {
    mySolver = HLLD;
#else
  } else if (solverString.compare("hllc") == 0) {
    mySolver = HLLC;
#endif
  } else if (solverString.compare("roe") == 0) {
    mySolver = ROE;
  } else {
    std::stringstream msg;
#if MHD == YES
    msg << "Unknown MHD solver type " << solverString;
#else
    msg << "Unknown HD solver type " << solverString;
#endif
    IDEFIX_ERROR(msg);
  }

  // Check if Hall is enabled
  if(input.CheckEntry(std::string(Phys::prefix),"hall")>=0) {
      // Check consistency
      if(mySolver != HLL )
        IDEFIX_ERROR("Hall effect is only compatible with HLL Riemann solver.");
  }

  // Shock flattening
  this->haveShockFlattening = input.CheckEntry(std::string(Phys::prefix),"shockFlattening")>=0;
  // Init shock flattening
  if(haveShockFlattening) {
    this->shockFlattening = ShockFlattening(hydro,input.Get<real>(std::string(Phys::prefix),"shockFlattening",0));
  }
}

template <typename Phys>
void RiemannSolver<Phys>::ShowConfig() {
  idfx::cout << "RiemannSolver: ";
  switch(mySolver) {
    case TVDLF:
      idfx::cout << "tvdlf." << std::endl;
      break;
    case HLL:
      idfx::cout << "hll." << std::endl;
      break;
    #if MHD==YES
      case HLLD:
        idfx::cout << "hlld." << std::endl;
        break;
    #else
      case HLLC:
        idfx::cout << "hllc." << std::endl;
        break;
    #endif
    case ROE:
      idfx::cout << "roe." << std::endl;
      break;
    default:
      IDEFIX_ERROR("Unknown Riemann solver");
  }

  if(haveShockFlattening) {
    idfx::cout << "Fluid: Shock Flattening ENABLED." << std::endl;
  }

}

#include "calcFlux.hpp"

#endif //FLUID_RIEMANNSOLVER_RIEMANNSOLVER_HPP_
