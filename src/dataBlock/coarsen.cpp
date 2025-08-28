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
  for(int dir = 0 ; dir < DIMENSIONS ; dir++) {
    if(mygrid->coarseningDirection[dir]) {
      IdefixHostArray2D<int> arr = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),
                                                                       coarseningLevel[dir]);
      for(int j = 0 ; j < arr.extent(0) ; j++) {
        for(int i = 0 ; i < arr.extent(1) ; i++) {
          if(std::isnan(arr(j,i))) {
            std::stringstream str;
            str << "Nan in grid coarsening levels" << std::endl;
            str << "at (i,j)=("<< i << "," << j << "): Coarsening level is NaN!" << std::endl;
            IDEFIX_ERROR(str);
          }
          if(arr(j,i) < 1) {
            std::stringstream str;
            str << "Coarsening level < 1!" << std::endl;
            str << "at (i,j)=("<< i << "," << j << "): ";
            str << "coarsening level= " << arr(j,i) << std::endl;
            IDEFIX_ERROR(str);
          }
          const int factor = 1 << (arr(j,i) - 1);
          if(np_int[dir] % factor != 0) {
            std::stringstream str;
            str << "Local grid size not divisible by coarsening level." << std::endl;
            str << "at (i,j)=("<< i << "," << j << "): ";
            str << "coarsening level= " << arr(j,i) << std::endl;
            str << np_int[dir] << " cannot be divided by 2^" << arr(j,i)-1;
            str << " = " << factor << std::endl;
            IDEFIX_ERROR(str);
          }
        }
      }
    }
  }

  idfx::popRegion();
}
