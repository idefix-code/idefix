// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2021 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include "tuner.hpp"
#include <string>

#include "idefix.hpp"
#include "dataBlock.hpp"
#include "input.hpp"
#include "timeIntegrator.hpp"
#include "setup.hpp"



namespace Tuner {
double testLoopType(DataBlock &data, Setup &setup, Input &input,
                    int numLoops, LoopPattern loopType) {
  idfx::defaultLoopPattern=loopType;

  TimeIntegrator tint(input,data);


  // Init the flow
  setup.InitFlow(data);
  data.SetBoundaries();

  // Launch a bunch of integrations
  int nint=0;
  tint.isSilent = true; // make the integration silent

  Kokkos::Timer timer;

  while(nint<numLoops) {
    tint.Cycle(data);
    nint++;
  }

  double tintegration = timer.seconds() / data.mygrid->np_int[IDIR] / data.mygrid->np_int[JDIR]
                            / data.mygrid->np_int[KDIR] / tint.getNcycles();

  //idfx::cout << "Perfs are " << 1/tintegration << " cell updates/second" << std::endl;

  return(1/tintegration);
}

// Convert lp integer into text
std::string LoopText(int lp) {
  std::string sout;

  switch(static_cast<LoopPattern>(lp)) {
    case LoopPattern::SIMDFOR:
      sout = std::string("SIMD");
      break;
    case LoopPattern::RANGE:
      sout = std::string("RANGE");
      break;
    case LoopPattern::MDRANGE:
      sout = std::string("MDRANGE");
      break;
    case LoopPattern::TPX:
      sout = std::string("TPX");
      break;
    case LoopPattern::TPTTRTVR:
      sout = std::string("TPTTRTVR");
      break;
    default:
      sout = std::string("UNKNOWN");
      break;
  }
  return(sout);
}
void tuneLoops(DataBlock &data, Setup &setup, Input &input, int numLoops) {
  LoopPattern oldLoop = idfx::defaultLoopPattern;
  real perfs = 0;
  LoopPattern bestLoop = oldLoop;

  // Do one dummy loop so that all is initialized
  idfx::cout << "Tuner: tuning idefix_loop... (this can take a while)" << std::endl;
  testLoopType(data,setup,input,numLoops,oldLoop);
  for(int i = 0 ; i != static_cast<int>(LoopPattern::UNDEFINED) ; i++) {
    //idfx::cout << "*************************************************" << std::endl;
    //idfx::cout << "Loop pattern " << LoopText(i) << "(" << i << ")" << std::endl;
    // Avoid SIMD for when using CUDA as this is not implemented
    #ifdef KOKKOS_ENABLE_CUDA
    if (static_cast<LoopPattern>(i) != LoopPattern::SIMDFOR) {
    #endif
      real thisPerfs = testLoopType(data,setup,input,numLoops,static_cast<LoopPattern>(i));
    #ifdef KOKKOS_ENABLE_CUDA
    } else {
      //idfx::cout << "Not implemented in Cuda, skipping" << std::endl;
    }
    #endif
    // Check if we have the best performances
    if(thisPerfs > perfs) {
      perfs=thisPerfs;
      bestLoop = static_cast<LoopPattern>(i);
    }
    // idfx::cout << "*************************************************" << std::endl;
  }
  idfx::cout << "Tuner: will be using " << LoopText(static_cast<int> (bestLoop))
             << " loops." << std::endl;

  idfx::defaultLoopPattern = bestLoop;

  // This is needed to re-init the datablock to its initial state
  TimeIntegrator tint(input,data);
}
} // namespace Tuner
