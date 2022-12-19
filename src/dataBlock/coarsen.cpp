// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************


#include "../idefix.hpp"
#include "dataBlock.hpp"
#include "dataBlockHost.hpp"
#include "fluid.hpp"

void DataBlock::Coarsen() {
  if(!haveGridCoarsening)  {
    IDEFIX_ERROR("DataBlock:Coarsen was called but grid coarsening is not enabled.");
  }
  ComputeGridCoarseningLevels();
  // This routine coarsen the *conservative* variables
  hydro->CoarsenFlow(hydro->Uc);
  #if MHD==YES
    hydro->CoarsenMagField(hydro->Vs);
  #endif
}

void DataBlock::EnrollGridCoarseningLevels(GridCoarseningFunc func) {
  if(!haveGridCoarsening) {
    IDEFIX_WARNING("DataBlock:EnrollCoarseningLevels was called but grid "
                    "coarsening is not enabled.");
  }
  this->gridCoarseningFunc = func;
}

void DataBlock::ComputeGridCoarseningLevels() {
  idfx::pushRegion("DataBlock::ComputeGridCoarseningLevels");
  static bool levelsHaveBeenComputedOnce = false;
  if((gridCoarseningFunc == NULL) && (haveGridCoarsening == GridCoarsening::dynamic)) {
    IDEFIX_ERROR("Dynamic grid Coarsening is enabled, "
                 "but no function has been enrolled to compute coarsening levels");
  }
  // if grid coarsening is enabled(=static), we compute the levels once
  // levels can be either initialised with the initial conditions, or with a dedicated
  // Coarsening function (if Enrollment has been called)
  if((haveGridCoarsening == GridCoarsening::enabled) && (!levelsHaveBeenComputedOnce)) {
    if(gridCoarseningFunc != NULL) {
      idfx::pushRegion("User-defined Coarsening function");
        gridCoarseningFunc(*this);
      idfx::popRegion();
      // We check the levels the first time this function is called. After that, there is no check!
      CheckCoarseningLevels();
    } else {
      IDEFIX_ERROR("Grid coarsening requires the enrollment of a grid coarsening function");
    }
    levelsHaveBeenComputedOnce = true;
  }
  if(haveGridCoarsening == GridCoarsening::dynamic) {
    idfx::pushRegion("User-defined Coarsening function");
      gridCoarseningFunc(*this);
    idfx::popRegion();
    levelsHaveBeenComputedOnce = true;
  }
  idfx::popRegion();
}

void DataBlock::CheckCoarseningLevels() {
  idfx::pushRegion("DataBlock::CheckCoarseningLevels()");
  // Check that the coarsening levels we have are valid
  // NB: this is a costly procedure, we can't repeat it at each loop!
  DataBlockHost d(*this);
  d.SyncFromDevice();
  for(int dir = 0 ; dir < DIMENSIONS ; dir++) {
    if(mygrid->coarseningDirection[dir]) {
      IdefixHostArray2D<int> arr = d.coarseningLevel[dir];
      const int Xt = (dir == IDIR ? JDIR : IDIR);
      const int Xb = (dir == KDIR ? JDIR : KDIR);
      for(int i = beg[Xt] ; i < end[Xt] ; i++) {
        for(int j = beg[Xb] ; j < end[Xb] ; j++) {
          if(arr(j,i) < 1) {
            std::stringstream str;
            str << "Incorrect grid coarsening levels" << std::endl;
            str << "at (i,j)=("<< i << "," << j << "): Coarsening level < 1!" << std::endl;
            IDEFIX_ERROR(str);
          }
          const int factor = 1 << (arr(j,i) - 1);
          if(np_int[dir] % factor != 0) {
            std::stringstream str;
            str << "local grid size not divisible by coarsening level" << std::endl;
            str << "at (i,j)=("<< i << "," << j << "): Coarsening level: ";
            str <<  arr(j,i) << std::endl;
            IDEFIX_ERROR(str);
          }
        }
      }
    }
  }

  idfx::popRegion();
}
