// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2022 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include "../idefix.hpp"
#include "dataBlock.hpp"

void DataBlock::InitFromGrid(Grid &grid, Input &input) {
  idfx::pushRegion("DataBlock::InitFromGrid");

  this->mygrid=&grid;

  // Make a local copy of the grid for future usage.
  GridHost gridHost(grid);
  gridHost.SyncFromDevice();

  nghost = std::vector<int>(3);
  np_int = std::vector<int>(3);
  np_tot = std::vector<int>(3);
  lbound = std::vector<BoundaryType>(3);
  rbound = std::vector<BoundaryType>(3);
  beg = std::vector<int>(3);
  end = std::vector<int>(3);
  gbeg = std::vector<int>(3);
  gend = std::vector<int>(3);

  xbeg = std::vector<real>(3);
  xend = std::vector<real>(3);

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
    xend[dir] = gridHost.xr[dir](gend[dir]-1);
  }

  // Allocate the required fields
  for(int dir = 0 ; dir < 3 ; dir++) {
    x.push_back(IdefixArray1D<real>("DataBlock_x",np_tot[dir]));
    xr.push_back(IdefixArray1D<real>("DataBlock_xr",np_tot[dir]));
    xl.push_back(IdefixArray1D<real>("DataBlock_xl",np_tot[dir]));
    dx.push_back(IdefixArray1D<real>("DataBlock_dx",np_tot[dir]));
    xgc.push_back(IdefixArray1D<real>("DataBlock_xgc",np_tot[dir]));

    A.push_back(IdefixArray3D<real>("DataBlock_A",
                                 np_tot[KDIR]+KOFFSET, np_tot[JDIR]+JOFFSET, np_tot[IDIR]+IOFFSET));
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

  // Initialize grid coarsening if needed
  if(grid.haveGridCoarsening != GridCoarsening::disabled) {
    this->haveGridCoarsening = grid.haveGridCoarsening;
    this->coarseningDirection = grid.coarseningDirection;
    this->coarseningLevel = std::vector<IdefixArray2D<int>>(3);

    for(int dir = 0 ; dir < 3 ; dir++) {
      if(coarseningDirection[dir]) {
        const int Xt = (dir == IDIR ? JDIR : IDIR);
        const int Xb = (dir == KDIR ? JDIR : KDIR);

        // Allocate coarsening level arrays
        coarseningLevel[dir] = IdefixArray2D<int>(
                                  "DataBlock_corseLevel",
                                  np_tot[Xb],
                                  np_tot[Xt]);
        // Make a local reference
        IdefixArray2D<int> coarseInit = coarseningLevel[dir];
        // Init coarsening level array to one everywhere
        idefix_for("init_coarsening", 0, np_tot[Xb], 0, np_tot[Xt],
                KOKKOS_LAMBDA(int j, int i) {
                  coarseInit(j,i) = 1;
                });
      }
    }
  }

  // Iniaitlize the geometry
  this->MakeGeometry();

  // Initialise the state containers
  // (by default, datablock only initialise the current state, which is a reference
  // to arrays in the daughter object

  this->states["current"] = StateContainer();

  // Initialize the hydro object attached to this datablock
  this->hydro.Init(input, grid, this);

  // Initialise Fargo if needed
  if(input.CheckBlock("Fargo")) {
    fargo.Init(input, this);
    this->haveFargo = true;
  }

  // Initialise gravity if needed
  if(input.CheckBlock("Gravity")) {
    gravity.Init(input, this);
    this->haveGravity = true;
  }

  idfx::popRegion();
}

void DataBlock::ResetStage() {
  this->hydro.ResetStage();
}

// Set the boundaries of the data structures in this datablock
void DataBlock::SetBoundaries() {
  if(haveGridCoarsening) {
    ComputeGridCoarseningLevels();
    hydro.CoarsenFlow(hydro.Vc);
    #if MHD==YES
      hydro.CoarsenMagField(hydro.Vs);
    #endif
  }
  hydro.boundary.SetBoundaries(t);
}



void DataBlock::ShowConfig() {
  if(idfx::psize>1) {
    idfx::cout << "DataBlock: this process grid size is " << std::endl;

    for(int dir = 0 ; dir < DIMENSIONS ; dir++) {
      idfx::cout << "\t Direction X" << (dir+1) << ": " << xbeg[dir] << "...." << np_int[dir]
        << "...." << xend[dir] << std::endl;
    }
  }
  hydro.ShowConfig();
  if(haveFargo) fargo.ShowConfig();
  if(haveGravity) gravity.ShowConfig();
}


real DataBlock::ComputeTimestep() {
  // Compute the timestep using all of the enabled modules in the current dataBlock

  // First with the hydro block
  auto InvDt = hydro.InvDt;
  real dt;
  idefix_reduce("Timestep_reduction",
          beg[KDIR], end[KDIR],
          beg[JDIR], end[JDIR],
          beg[IDIR], end[IDIR],
          KOKKOS_LAMBDA (int k, int j, int i, real &dtmin) {
                  dtmin=FMIN(ONE_F/InvDt(k,j,i),dtmin);
              },
          Kokkos::Min<real>(dt));
  Kokkos::fence();
  return(dt);
}
