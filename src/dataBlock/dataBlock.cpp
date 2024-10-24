// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include <algorithm>
#include <memory>
#include <string>
#include "idefix.hpp"
#include "dataBlock.hpp"
#include "fluid.hpp"
#include "gravity.hpp"
#include "planetarySystem.hpp"
#include "vtk.hpp"
#include "dump.hpp"
#ifdef WITH_HDF5
#include "xdmf.hpp"
#endif

DataBlock::DataBlock(Grid &grid, Input &input) {
  idfx::pushRegion("DataBlock::DataBlock");

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
    xend[dir] = gridHost.xr[dir](gend[dir]-1);
  }

  // Allocate the required fields
  std::string label;
  for(int dir = 0 ; dir < 3 ; dir++) {
    label = "DataBlock_x" + std::to_string(dir);
    x[dir] = IdefixArray1D<real>(label, np_tot[dir]);

    label = "DataBlock_xr" + std::to_string(dir);
    xr[dir] = IdefixArray1D<real>(label,np_tot[dir]);

    label = "DataBlock_xl" + std::to_string(dir);
    xl[dir] = IdefixArray1D<real>(label,np_tot[dir]);

    label = "DataBlock_dx" + std::to_string(dir);
    dx[dir] = IdefixArray1D<real>(label,np_tot[dir]);

    label = "DataBlock_xgc" + std::to_string(dir);
    xgc[dir] = IdefixArray1D<real>(label,np_tot[dir]);

    label = "DataBlock_A" + std::to_string(dir);
    A[dir] = IdefixArray3D<real>(label,
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

  // Initialize our sub-domain
  this->ExtractSubdomain();

  // Initialize the geometry
  this->MakeGeometry();

  // Initialise the state containers
  // (by default, datablock only initialise the current state, which is a reference
  // to arrays in the daughter object

  this->states["current"] = StateContainer();

  // Initialize the Dump object
  this->dump = std::make_unique<Dump>(input, this);

  // Initialize the VTK object
  this->vtk = std::make_unique<Vtk>(input, this, "data");

  // Init XDMF objects for HDF5 outputs
  #ifdef WITH_HDF5
    this->xdmf= std::make_unique<Xdmf>(input,this);
  #endif


  // Initialize the hydro object attached to this datablock
  this->hydro = std::make_unique<Fluid<DefaultPhysics>>(grid, input, this);

  // Initialise Fargo if needed
  if(input.CheckBlock("Fargo")) {
    this->fargo = std::make_unique<Fargo>(input, DefaultPhysics::nvar, this);
    this->haveFargo = true;
  }

  // initialise planets if needed
  if(input.CheckBlock("Planet")) {
    this->planetarySystem = std::make_unique<PlanetarySystem>(input, this);
    this->haveplanetarySystem = true;
  }

  // Initialise gravity if needed (automatically if planets are present)
  if(input.CheckBlock("Gravity") || haveplanetarySystem) {
    this->gravity = std::make_unique<Gravity>(input, this);
    this->haveGravity = true;
  }

  // Initialise dust grains if needed
  if(input.CheckBlock("Dust")) {
    haveDust = true;
    int nSpecies = input.Get<int>("Dust","nSpecies",0);
    for(int i = 0 ; i < nSpecies ; i++) {
      dust.emplace_back(std::make_unique<Fluid<DustPhysics>>(grid, input, this, i));
    }
  }
  // Register variables that need to be saved in case of restart dump
  dump->RegisterVariable(&t, "time");
  dump->RegisterVariable(&dt, "dt");

  idfx::popRegion();
}

/**
 * @brief Construct a new Data Block as a subgrid
 *
 * @param subgrid : subgrid from which the local datablock should be extracted from
 */
DataBlock::DataBlock(SubGrid *subgrid) {
  idfx::pushRegion("DataBlock:DataBlock(SubGrid)");
  Grid *grid = subgrid->parentGrid;
  this->mygrid = subgrid->grid.get();

  // Make a local copy of the grid for future usage.
  GridHost gridHost(*grid);
  gridHost.SyncFromDevice();


  // Get the number of points from the parent grid object
  for(int dir = 0 ; dir < 3 ; dir++) {
    nghost[dir] = grid->nghost[dir];
    // Domain decomposition: decompose the full domain size in grid by the number of processes
    // in that direction
    np_int[dir] = grid->np_int[dir]/grid->nproc[dir];
    np_tot[dir] = np_int[dir]+2*nghost[dir];

    // Boundary conditions
    if (grid->xproc[dir]==0) {
      lbound[dir] = grid->lbound[dir];
      if(lbound[dir]==axis) this->haveAxis = true;
    } else {
      lbound[dir] = internal;
    }

    if (grid->xproc[dir] == grid->nproc[dir]-1) {
      rbound[dir] = grid->rbound[dir];
      if(rbound[dir]==axis) this->haveAxis = true;
    } else {
      rbound[dir] = internal;
    }

    beg[dir] = grid->nghost[dir];
    end[dir] = grid->nghost[dir]+np_int[dir];

    // Where does this datablock starts and end in the grid?
    // This assumes even distribution of points between procs
    gbeg[dir] = grid->nghost[dir] + grid->xproc[dir]*np_int[dir];
    gend[dir] = grid->nghost[dir] + (grid->xproc[dir]+1)*np_int[dir];

    // Local start and end of current datablock
    xbeg[dir] = gridHost.xl[dir](gbeg[dir]);
    xend[dir] = gridHost.xr[dir](gend[dir]-1);
  }

  // Next we erase the info the slice direction
  int refdir = subgrid->direction;
  nghost[refdir] = 0;
  np_int[refdir] = 1;
  np_tot[refdir] = 1;
  beg[refdir] = 0;
  end[refdir] = 1;
  gbeg[refdir] = subgrid->index;
  gend[refdir] = subgrid->index+1;
  xbeg[refdir] = gridHost.x[refdir](subgrid->index);
  xend[refdir] = gridHost.x[refdir](subgrid->index);

  // Reset gbeg/gend so that they don't refer anymore to the parent grid
  gbeg[refdir] = 0;
  gend[refdir] = 1;


  // Allocate the required fields (only a limited set for a datablock from a subgrid)
  std::string label;
  for(int dir = 0 ; dir < 3 ; dir++) {
    label = "DataBlock_x" + std::to_string(dir);
    x[dir] = IdefixArray1D<real>(label, np_tot[dir]);

    label = "DataBlock_xr" + std::to_string(dir);
    xr[dir] = IdefixArray1D<real>(label,np_tot[dir]);

    label = "DataBlock_xl" + std::to_string(dir);
    xl[dir] = IdefixArray1D<real>(label,np_tot[dir]);

    label = "DataBlock_dx" + std::to_string(dir);
    dx[dir] = IdefixArray1D<real>(label,np_tot[dir]);
  }

  // Initialize our sub-domain
  this->ExtractSubdomain();


  idfx::popRegion();
}

void DataBlock::ResetStage() {
  this->hydro->ResetStage();
  if(haveDust) {
    for(int i = 0 ; i < dust.size() ; i++) {
      dust[i]->ResetStage();
    }
  }
}

void DataBlock::ConsToPrim() {
  this->hydro->ConvertConsToPrim();
  if(haveDust) {
    for(int i = 0 ; i < dust.size() ; i++) {
      dust[i]->ConvertConsToPrim();
    }
  }
}

void DataBlock::PrimToCons() {
  this->hydro->ConvertPrimToCons();
  if(haveDust) {
    for(int i = 0 ; i < dust.size() ; i++) {
      dust[i]->ConvertPrimToCons();
    }
  }
}

// Set the boundaries of the data structures in this datablock
void DataBlock::SetBoundaries() {
  if(haveGridCoarsening) {
    ComputeGridCoarseningLevels();
    hydro->CoarsenFlow(hydro->Vc);
    #if MHD==YES
      hydro->CoarsenMagField(hydro->Vs);
    #endif
    if(haveDust) {
      for(int i = 0 ; i < dust.size() ; i++) {
        dust[i]->CoarsenFlow(dust[i]->Vc);
      }
    }
  }
  if(haveDust) {
    for(int i = 0 ; i < dust.size() ; i++) {
      dust[i]->boundary->SetBoundaries(t);
    }
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
  if(haveplanetarySystem) planetarySystem->ShowConfig();
  if(haveGravity) gravity->ShowConfig();
  if(haveUserStepFirst) idfx::cout << "DataBlock: User's first step has been enrolled."
                                   << std::endl;
  if(haveUserStepLast) idfx::cout << "DataBlock: User's last step has been enrolled."
                                  << std::endl;
  if(haveDust) {
    idfx::cout << "DataBlock: evolving " << dust.size() << " dust species." << std::endl;
    // Only show the config the first dust specie
    dust[0]->ShowConfig();
    /*
    for(int i = 0 ; i < dust.size() ; i++) {
      dust[i]->ShowConfig();
    }*/
  }
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
  if(haveDust) {
    for(int n = 0 ; n < dust.size() ; n++) {
      real dtDust;
      auto InvDt = dust[n]->InvDt;
      idefix_reduce("Timestep_reduction_dust",
          beg[KDIR], end[KDIR],
          beg[JDIR], end[JDIR],
          beg[IDIR], end[IDIR],
          KOKKOS_LAMBDA (int k, int j, int i, real &dtmin) {
                  dtmin=FMIN(ONE_F/InvDt(k,j,i),dtmin);
              },
          Kokkos::Min<real>(dtDust));
      dt = std::min(dt,dtDust);
    }
  }
  Kokkos::fence();
  return(dt);
}

// Recompute magnetic fields from vector potential in dedicated fluids
void DataBlock::DeriveVectorPotential() {
  if constexpr(DefaultPhysics::mhd) {
    #ifdef EVOLVE_VECTOR_POTENTIAL
      hydro->emf->ComputeMagFieldFromA(hydro->Ve, hydro->Vs);
    #endif
  }
}

void DataBlock::LaunchUserStepLast() {
  if(haveUserStepLast) {
    idfx::pushRegion("User::UserStepLast");
    if(userStepLast != nullptr)
      userStepLast(*this, this->t, this->dt);
    else
      IDEFIX_ERROR("UserStepLast not properly initialized");
    idfx::popRegion();
  }
}

void DataBlock::LaunchUserStepFirst() {
  if(haveUserStepFirst) {
    idfx::pushRegion("User::UserStepFirst");
    if(userStepFirst != nullptr)
      userStepFirst(*this, this->t, this->dt);
    else
      IDEFIX_ERROR("UserStepLast not properly initialized");
    idfx::popRegion();
  }
}

void DataBlock::EnrollUserStepLast(StepFunc func) {
  haveUserStepLast = true;
  userStepLast = func;
}

void DataBlock::EnrollUserStepFirst(StepFunc func) {
  haveUserStepFirst = true;
  userStepFirst = func;
}
