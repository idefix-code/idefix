// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2021 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include "../idefix.hpp"
#include "dataBlock.hpp"

DataBlock::DataBlock() {
  // Do nothing
}

void DataBlock::InitFromGrid(Grid &grid, Input &input) {
  // This initialisation is only valid for *serial*
  // MPI initialisation will involve domain decomposition of grids into DataBlocks

  idfx::pushRegion("DataBlock::InitFromGrid");

  this->mygrid=&grid;

  // Make a local copy of the grid for future usage.
  GridHost gridHost(grid);
  gridHost.SyncFromDevice();

  // Get the number of points from the parent grid object
  for(int dir = 0 ; dir < 3 ; dir++) {
    nghost[dir] = grid.nghost[dir];
    // Domain decomposition: decompose the full domain size in grid by the number of processes
    // in that direction
    np_int[dir] = grid.np_int[dir]/grid.nproc[dir];
    np_tot[dir] = np_int[dir]+2*nghost[dir];

    // Boundary conditions
    if (grid.xproc[dir]==0) {
      lbound[dir] = grid.lbound[dir];
      if(lbound[dir]==axis) this->haveAxis = true;
    } else {
      lbound[dir] = internal;
    }

    if (grid.xproc[dir] == grid.nproc[dir]-1) {
      rbound[dir] = grid.rbound[dir];
      if(rbound[dir]==axis) this->haveAxis = true;
    } else {
      rbound[dir] = internal;
    }

    beg[dir] = grid.nghost[dir];
    end[dir] = grid.nghost[dir]+np_int[dir];

    // Where does this datablock starts and end in the grid?
    // This assumes even distribution of points between procs
    gbeg[dir] = grid.nghost[dir] + grid.xproc[dir]*np_int[dir];
    gend[dir] = grid.nghost[dir] + (grid.xproc[dir]+1)*np_int[dir];

    // Local start and end of current datablock
    xbeg[dir] = gridHost.xl[dir](gbeg[dir]);
    xend[dir] = gridHost.xl[dir](gend[dir]);
  }

  if(idfx::psize>1) {
    idfx::cout << "DataBlock::InitFromGrid: local size is " << std::endl;

    for(int dir = 0 ; dir < DIMENSIONS ; dir++) {
      idfx::cout << "\t Direction X" << (dir+1) << ": " << xbeg[dir] << "...." << np_int[dir]
        << "...." << xend[dir] << std::endl;
    }
  }

  // Allocate the required fields
  for(int dir = 0 ; dir < 3 ; dir++) {
    x[dir] = IdefixArray1D<real>("DataBlock_x",np_tot[dir]);
    xr[dir] = IdefixArray1D<real>("DataBlock_xr",np_tot[dir]);
    xl[dir] = IdefixArray1D<real>("DataBlock_xl",np_tot[dir]);
    dx[dir] = IdefixArray1D<real>("DataBlock_dx",np_tot[dir]);
    xgc[dir] = IdefixArray1D<real>("DataBlock_xgc",np_tot[dir]);

    A[dir] = IdefixArray3D<real>("DataBlock_A",
                                 np_tot[KDIR]+KOFFSET, np_tot[JDIR]+JOFFSET, np_tot[IDIR]+IOFFSET);
  }

  dV = IdefixArray3D<real>("DataBlock_dV",np_tot[KDIR],np_tot[JDIR],np_tot[IDIR]);

#if GEOMETRY == SPHERICAL
  rt = IdefixArray1D<real>("DataBlock_rt",np_tot[IDIR]);
  sinx2m = IdefixArray1D<real>("DataBlock_sinx2m",np_tot[JDIR]);
  tanx2m = IdefixArray1D<real>("DataBlock_tanx2m",np_tot[JDIR]);
  sinx2 = IdefixArray1D<real>("DataBlock_sinx2",np_tot[JDIR]);
  tanx2 = IdefixArray1D<real>("DataBlock_tanx2",np_tot[JDIR]);
  dmu = IdefixArray1D<real>("DataBlock_dmu",np_tot[JDIR]);
#endif



  // Copy the relevant part of the coordinate system to the datablock
  for(int dir = 0 ; dir < 3 ; dir++) {
    int offset=gbeg[dir]-beg[dir];

    IdefixArray1D<real> x_input = grid.x[dir];
    IdefixArray1D<real> x_output= x[dir];
    IdefixArray1D<real> xr_input = grid.xr[dir];
    IdefixArray1D<real> xr_output= xr[dir];
    IdefixArray1D<real> xl_input = grid.xl[dir];
    IdefixArray1D<real> xl_output= xl[dir];
    IdefixArray1D<real> dx_input = grid.dx[dir];
    IdefixArray1D<real> dx_output= dx[dir];

    idefix_for("coordinates",0,np_tot[dir],
      KOKKOS_LAMBDA (int i) {
        x_output(i)  = x_input(i+offset);
        xr_output(i) = xr_input(i+offset);
        xl_output(i) = xl_input(i+offset);
        dx_output(i) = dx_input(i+offset);
      }
    );
  }

  // Iniaitlize the geometry
  this->MakeGeometry();

  // Initialize the hydro object attached to this datablock
  this->hydro.Init(input, grid, this);

#if RKL_ENABLED == YES
  // Initialize the RKL object
  rkl.Init(this);
#endif

  // Init MPI stack when needed
#ifdef WITH_MPI
  ////////////////////////////////////////////////////////////////////////////
  // Init variable mappers
  // The variable mapper list all of the variable which are exchanged in MPI boundary calls
  // This is required since we skip some of the variables in Vc to limit the amount of data
  // being exchanged
  #if MHD == YES
  int mapNVars = NVAR - DIMENSIONS; // We will not send magnetic field components which are in Vs
  #else
  int mapNVars = NVAR;
  #endif

  IdefixArray1D<int> mapVars("mapVars",mapNVars);

  // Create a host mirror
  IdefixArray1D<int>::HostMirror mapVarsHost = Kokkos::create_mirror_view(mapVars);
  // Init the list of variables we will exchange in MPI routines
  int ntarget = 0;
  for(int n = 0 ; n < mapNVars ; n++) {
    mapVarsHost(n) = ntarget;
    ntarget++;
    #if MHD == YES
      // Skip centered field components if they are also defined in Vs
      #if DIMENSIONS >= 1
        if(ntarget==BX1) ntarget++;
      #endif
      #if DIMENSIONS >= 2
        if(ntarget==BX2) ntarget++;
      #endif
      #if DIMENSIONS == 3
        if(ntarget==BX3) ntarget++;
      #endif
    #endif
  }
  // Synchronize the mapVars
  Kokkos::deep_copy(mapVars,mapVarsHost);
  this->mpi = new Mpi(this, mapVars, mapNVars, true);
#endif // MPI
  idfx::popRegion();
}

void DataBlock::ResetStage() {
  this->hydro.ResetStage();
}

DataBlock::~DataBlock() {
  #ifdef WITH_MPI
    if(mpi != NULL) {
      delete mpi;
    }
  #endif
}