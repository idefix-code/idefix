// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2021 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************


#include "idefix.hpp"
#include "dataBlock.hpp"
#include "input.hpp"
#include "timeIntegrator.hpp"
#include "setup.hpp"

void testLoopType(DataBlock &data, Setup &setup, Input &input, int numLoops, LoopPattern loopType) {
  idfx::defaultLoopPattern=loopType;

  TimeIntegrator tint(input,data);

  // Init the flow
  setup.InitFlow(data);
  data.SetBoundaries();

  // Launch a bunch of integrations
  int nint=0;

  Kokkos::Timer timer;

  while(nint<numLoops) {
    tint.Cycle(data);
    nint++;
  }

  double tintegration = timer.seconds() / data.mygrid->np_int[IDIR] / data.mygrid->np_int[JDIR]
                            / data.mygrid->np_int[KDIR] / tint.getNcycles();

  idfx::cout << "Perfs are " << 1/tintegration << " cell updates/second" << std::endl;

}
void testLoops(DataBlock &data, Setup &setup, Input &input, int numLoops) {
  LoopPattern oldLoop = idfx::defaultLoopPattern;
  for(int i = 0 ; i != static_cast<int>(LoopPattern::UNDEFINED) ; i++) {
    idfx::cout << "*************************************************" << std::endl;
    idfx::cout << "Loop pattern " << i << std::endl;
    testLoopType(data,setup,input,numLoops,static_cast<LoopPattern>(i));
    idfx::cout << "*************************************************" << std::endl;
  }
  idfx::defaultLoopPattern = oldLoop;
}
