// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include "../idefix.hpp"
#include "dataBlock.hpp"
#include "fluid.hpp"
#include "vtk.hpp"
#include "dump.hpp"

DataBlock::DataBlock(Grid &grid, Input &input) {
  idfx::pushRegion("DataBlock::DataBlock");

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
  std::string label;
  for(int dir = 0 ; dir < 3 ; dir++) {
    label = "DataBlock_x" + std::to_string(dir);
    x.push_back(IdefixArray1D<real>(label, np_tot[dir]));

    label = "DataBlock_xr" + std::to_string(dir);
    xr.push_back(IdefixArray1D<real>(label,np_tot[dir]));

    label = "DataBlock_xl" + std::to_string(dir);
    xl.push_back(IdefixArray1D<real>(label,np_tot[dir]));

    label = "DataBlock_dx" + std::to_string(dir);
    dx.push_back(IdefixArray1D<real>(label,np_tot[dir]));

    label = "DataBlock_xgc" + std::to_string(dir);
    xgc.push_back(IdefixArray1D<real>(label,np_tot[dir]));

    label = "DataBlock_A" + std::to_string(dir);
    A.push_back(IdefixArray3D<real>(label,
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

  // Iniaitlize the geometry
  this->MakeGeometry();

  // Initialise the state containers
  // (by default, datablock only initialise the current state, which is a reference
  // to arrays in the daughter object

  this->states["current"] = StateContainer();

  // Initialize the Dump object
  this->dump = std::make_unique<Dump>(this);

  // Initialize the VTK object
  this->vtk = std::make_unique<Vtk>(input, this);

  // Initialize the hydro object attached to this datablock
  this->hydro = std::make_unique<Fluid<Physics>>(grid, input, this);

  // Initialise Fargo if needed
  if(input.CheckBlock("Fargo")) {
    this->fargo = std::make_unique<Fargo>(input, Physics::nvar, this);
    this->haveFargo = true;
  }

  // initialise planets if needed
  if(input.CheckBlock("Planet")) {
    planetarySystem.Init(input, this);
    this->haveplanetarySystem = true;
  }

  // Initialise gravity if needed (automatically if planets are present)
  if(input.CheckBlock("Gravity") || haveplanetarySystem) {
    gravity.Init(input, this);
    this->haveGravity = true; // TODO(mauxionj): why do it here and in init gravity ?
  }

  // Register variables that need to be saved in case of restart dump
  dump->RegisterVariable(&t, "time");
  dump->RegisterVariable(&dt, "dt");

  idfx::popRegion();
}

void DataBlock::ResetStage() {
  this->hydro->ResetStage();
}

// Set the boundaries of the data structures in this datablock
void DataBlock::SetBoundaries() {
  if(haveGridCoarsening) {
    ComputeGridCoarseningLevels();
    hydro->CoarsenFlow(hydro->Vc);
    #if MHD==YES
      hydro->CoarsenMagField(hydro->Vs);
    #endif
  }
  hydro->boundary->SetBoundaries(t);
}



void DataBlock::ShowConfig() {
  if(idfx::psize>1) {
    idfx::cout << "DataBlock: this process grid size is " << std::endl;

    for(int dir = 0 ; dir < DIMENSIONS ; dir++) {
      idfx::cout << "\t Direction X" << (dir+1) << ": " << xbeg[dir] << "...." << np_int[dir]
        << "...." << xend[dir] << std::endl;
    }
  }
  hydro->ShowConfig();
  if(haveFargo) fargo->ShowConfig();
  if(haveplanetarySystem) planetarySystem.ShowConfig();
  if(haveGravity) gravity.ShowConfig();
}


real DataBlock::ComputeTimestep() {
  // Compute the timestep using all of the enabled modules in the current dataBlock

  // First with the hydro block
  auto InvDt = hydro->InvDt;
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

// Recompute magnetic fields from vector potential in dedicated fluids
void DataBlock::DeriveVectorPotential() {
  if constexpr(Physics::mhd) {
    #ifdef EVOLVE_VECTOR_POTENTIAL
      hydro->emf->ComputeMagFieldFromA(hydro->Ve, hydro->Vs);
    #endif
  }
}
