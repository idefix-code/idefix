// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2022 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include <string>

#include "idefix.hpp"
#include "gridHost.hpp"
#include "grid.hpp"

Grid::Grid(Input &input) {
  idfx::pushRegion("Grid::Grid(Input)");

  xbeg = std::vector<real>(3);
  xend = std::vector<real>(3);

  nghost = std::vector<int>(3);

  // Get grid size from input file, block [Grid]
  int npoints[3];
  for(int dir = 0 ; dir < 3 ; dir++) {
    npoints[dir] = 1;
    nghost[dir] = 0;
    std::string label = std::string("X")+std::to_string(dir+1)+std::string("-grid");
    int numPatch = input.Get<int>("Grid",label,0);

    if(dir<DIMENSIONS) {
      #if ORDER < 4
        nghost[dir] = 2;
      #else
        nghost[dir] = 3;
      #endif
      npoints[dir] = 0;
      for(int patch = 0; patch < numPatch ; patch++) {
        npoints[dir] += input.Get<int>("Grid",label,2+3*patch );
      }
    }
  }

  np_tot = std::vector<int>(3);
  np_int = std::vector<int>(3);
  lbound = std::vector<BoundaryType>(3);
  rbound = std::vector<BoundaryType>(3);

  for(int dir = 0 ; dir < 3 ; dir++) {
    np_tot[dir] = npoints[dir] + 2*nghost[dir];
    np_int[dir] = npoints[dir];

    std::string label = std::string("X")+std::to_string(dir+1)+std::string("-beg");
    std::string boundary = input.Get<std::string>("Boundary",label,0);

    if(boundary.compare("outflow") == 0) {
      lbound[dir] = outflow;
    } else if(boundary.compare("periodic") == 0) {
      lbound[dir] = periodic;
    } else if(boundary.compare("reflective") == 0) {
      lbound[dir] = reflective;
    } else if(boundary.compare("internal") == 0) {
      lbound[dir] = internal;
    } else if(boundary.compare("shearingbox") == 0) {
      lbound[dir] = shearingbox;
    } else if(boundary.compare("axis") == 0) {
      if(dir!= JDIR) {
        IDEFIX_ERROR("Axis Boundaries are only applicable to X2");
      }
      lbound[dir] = axis;
      haveAxis = true;
    } else if(boundary.compare("userdef") == 0) {
      lbound[dir] = userdef;
    } else {
      std::stringstream msg;
      msg << "Unknown boundary type " << boundary;
      IDEFIX_ERROR(msg);
    }

    label = std::string("X")+std::to_string(dir+1)+std::string("-end");
    boundary = input.Get<std::string>("Boundary",label,0);
    if(boundary.compare("outflow") == 0) {
      rbound[dir] = outflow;
    } else if(boundary.compare("periodic") == 0) {
      rbound[dir] = periodic;
    } else if(boundary.compare("reflective") == 0) {
      rbound[dir] = reflective;
    } else if(boundary.compare("internal") == 0) {
      rbound[dir] = internal;
    } else if(boundary.compare("shearingbox") == 0) {
      rbound[dir] = shearingbox;
    } else if(boundary.compare("axis") == 0) {
      if(dir!= JDIR) {
        IDEFIX_ERROR("Axis Boundaries are only applicable to X2");
      }
      rbound[dir] = axis;
      haveAxis = true;
    } else if(boundary.compare("userdef") == 0) {
      rbound[dir] = userdef;
    } else {
      std::stringstream msg;
      msg << "Unknown boundary type " << boundary;
      IDEFIX_ERROR(msg);
    }
  }

  // Allocate the grid structure on device. Initialisation will come from GridHost
  x = std::vector<IdefixArray1D<real>>(3);
  xr = std::vector<IdefixArray1D<real>>(3);
  xl = std::vector<IdefixArray1D<real>>(3);
  dx = std::vector<IdefixArray1D<real>>(3);
  for(int dir = 0 ; dir < 3 ; dir++) {
    x[dir] = IdefixArray1D<real>("Grid_x",np_tot[dir]);
    xr[dir] = IdefixArray1D<real>("Grid_xr",np_tot[dir]);
    xl[dir] = IdefixArray1D<real>("Grid_xl",np_tot[dir]);
    dx[dir] = IdefixArray1D<real>("Grid_dx",np_tot[dir]);
  }

  // Allocate proc structure (by default: one proc in each direction, size one)
  nproc = std::vector<int> (3);
  xproc = std::vector<int> (3);
  for(int i=0 ; i < 3; i++) {
    nproc[i] = 1;
    xproc[i] = 0;
  }

#ifdef WITH_MPI
  // Domain decomposition required for the grid

  // Init periodicity array
  int period[3];
  for(int i=0 ; i < 3 ; i++)
    period[i] = 0;

  // Check if the dec option has been passed when number of procs > 1
  if(idfx::psize>1) {
    if(input.CheckEntry("CommandLine","dec")  != DIMENSIONS) {
      // No command line decomposition, make auto-decomposition if possible
      // (only when nproc and dimensions are powers of 2, and in 1D)
      if(DIMENSIONS == 1) {
        nproc[0] = idfx::psize;
      } else {
        if(!isPow2(idfx::psize))
          IDEFIX_ERROR(
            "Automatic domain decomposition requires the number of processes to be a power of 2. "
            "Alternatively, set a manual decomposition with -dec"
          );
        for(int dir = 0; dir < 3; dir++) {
          if(!isPow2(np_int[dir]))
            IDEFIX_ERROR(
              "Automatic domain decomposition requires nx1, nx2 and nx3 to be powers of 2. "
              "Alternatively, set a manual decomposition with -dec"
            );
        }
        makeDomainDecomposition();
      }
    } else {
      // Manual domain decomposition (with -dec option)
      int ntot=1;
      for(int dir=0 ; dir < DIMENSIONS; dir++) {
        nproc[dir] = input.Get<int>("CommandLine","dec",dir);
        // Check that the dimension is effectively divisible by number of procs
        if(np_int[dir] % nproc[dir])
          IDEFIX_ERROR("Grid size must be a multiple of the domain decomposition");
        // Count the total number of procs we'll need for the specified domain decomposition
        ntot = ntot * nproc[dir];
      }
      if(ntot != idfx::psize) {
        std::stringstream msg;
        msg << "The number of MPI process (" << idfx::psize
            << ") does not match your domain decomposition (";
        for(int dir=0 ; dir < DIMENSIONS; dir++) {
          msg << nproc[dir];
          if(dir<DIMENSIONS-1) msg << ", ";
        }
        msg << ").";
        IDEFIX_ERROR(msg);
      }
    }
  }

  // Add periodicity indications
  for(int dir=0 ; dir < DIMENSIONS; dir++) {
    if(rbound[dir] == periodic || rbound[dir] == shearingbox) period[dir] = 1;
  }

  // Create cartesian communicator along with cartesian coordinates.
  MPI_Cart_create(MPI_COMM_WORLD, 3, nproc.data(), period, 0, &CartComm);
  MPI_Cart_coords(CartComm, idfx::prank, 3, xproc.data());

  MPI_Barrier(MPI_COMM_WORLD);


  if(haveAxis) {
      // create axis communicator to be able to exchange data over the axis
      // (only retain the phi dimension)
      int remainDims[3] = {false, false, true};
      MPI_SAFE_CALL(MPI_Cart_sub(CartComm, remainDims, &AxisComm));
  }
#endif

  // init coarsening
  if(input.CheckEntry("Grid","coarsening")>=0) {
    std::string coarsenType = input.Get<std::string>("Grid","coarsening",0);
    if(coarsenType.compare("static")==0) {
      this->haveGridCoarsening = GridCoarsening::enabled;
    } else if(coarsenType.compare("dynamic")==0) {
      this->haveGridCoarsening = GridCoarsening::dynamic;
    } else {
      std::stringstream msg;
      msg << "Grid coarsening can only be static or dynamic. I got: " << coarsenType;
      IDEFIX_ERROR(msg);
    }
    this->coarseningDirection = std::vector<bool>(3, false);
    int directions = input.CheckEntry("Grid","coarsening");
    for(int i = 1 ; i < directions ; i++) {
      std::string dirname = input.Get<std::string>("Grid","coarsening",i);
      if(dirname.compare("X1")==0) {
        coarseningDirection[IDIR] = true;
      } else if(dirname.compare("X2")==0) {
        coarseningDirection[JDIR] = true;
      } else if(dirname.compare("X3")==0) {
        coarseningDirection[KDIR] = true;
      } else {
        std::stringstream msg;
        msg << "Grid coarsening direction can only be X1, X2 and/or X3. I got: " << dirname;
        IDEFIX_ERROR(msg);
      }
    }

    this->haveGridCoarsening = GridCoarsening::enabled;
  }
  idfx::popRegion();
}

bool Grid::isPow2(int n) {
  return( (n & (n-1)) == 0);
}

// Produce a domain decomposition in nProc assuming psize and np_int[] are powers of 2
void Grid::makeDomainDecomposition() {
  // initialize the routine
  int nleft=idfx::psize;
  int nlocal[3];
  for(int dir = 0; dir < 3; dir++) {
    nproc[dir] = 1;
    nlocal[dir] = np_int[dir];
  }

  // Loop
  while(nleft>1) {
    // Find the direction where there is a maximum of point
    int dirmax;
    int nmax=1;
    int ndir=2;

    for(int dir = ndir; dir >= 0; dir--) {
      // We do this loop backward so that we divide the domain first in the last dimension
      // (better for cache optimisation)
      if(nlocal[dir]>nmax ) {
        nmax=nlocal[dir];
        dirmax=dir;
      }
    }
    // At this point, we have nmax points in direction dirmax,
    // which is the direction we're going to divide by 2
    if(nmax==1)
      IDEFIX_ERROR("Your domain size is too small to be decomposed "
                   "on this number of MPI processes");
    nlocal[dirmax]=nlocal[dirmax]/2;
    nproc[dirmax]=nproc[dirmax]*2;
    nleft=nleft/2;
  }
}
/*
Grid& Grid::operator=(const Grid& grid) {
    for(int dir = 0 ; dir < 3 ; dir++) {
        x[dir] = grid.x[dir];
        xr[dir] = grid.xr[dir];
        xl[dir] = grid.xl[dir];
        dx[dir] = grid.dx[dir];
        xbeg[dir] = grid.xbeg[dir];
        xend[dir] = grid.xend[dir];
        np_tot[dir] = grid.np_tot[dir];
        np_int[dir] = grid.np_int[dir];
        nghost[dir] = grid.nghost[dir];
        lbound[dir] = grid.lbound[dir];
        rbound[dir] = grid.rbound[dir];
    }

    return *this;
}
*/

void Grid::ShowConfig() {
  idfx::cout << "Grid: full grid size is " << std::endl;
  for(int dir = 0 ; dir < DIMENSIONS ; dir++) {
    std::string lboundString, rboundString;
      switch(lbound[dir]) {
        case outflow:
          lboundString="outflow";
          break;
        case reflective:
          lboundString="reflective";
          break;
        case periodic:
          lboundString="periodic";
          break;
        case internal:
          lboundString="internal";
          break;
        case shearingbox:
          lboundString="shearingbox";
          break;
        case axis:
          lboundString="axis";
          break;
        case userdef:
          lboundString="userdef";
          break;
        default:
          lboundString="unknown";
      }
      switch(rbound[dir]) {
        case outflow:
          rboundString="outflow";
          break;
        case reflective:
          rboundString="reflective";
          break;
        case periodic:
          rboundString="periodic";
          break;
        case internal:
          rboundString="internal";
          break;
        case shearingbox:
          rboundString="shearingbox";
          break;
        case axis:
          rboundString="axis";
          break;
        case userdef:
          rboundString="userdef";
          break;
        default:
          rboundString="unknown";
      }

      idfx::cout << "\t Direction X" << (dir+1) << ": " << lboundString << "\t" << xbeg[dir]
                 << "...." << np_int[dir] << "...." << xend[dir] << "\t" << rboundString
                 << std::endl;
  }
  #ifdef WITH_MPI
    idfx::cout << "Grid: MPI domain decomposition is (";
    for(int dir = 0 ; dir < DIMENSIONS ; dir++) {
      idfx::cout << " " << nproc[dir] << " ";
    }
    idfx::cout << ")" << std::endl;
    idfx::cout << "Grid: Current MPI proc coordinates (";

    for(int dir = 0; dir < 3; dir++) {
      idfx::cout << xproc[dir];
      if(dir < 2) idfx::cout << ", ";
    }
    idfx::cout << ")" << std::endl;
  #endif
  if(haveGridCoarsening) {
    if(haveGridCoarsening == GridCoarsening::enabled ) {
      idfx::cout << "Grid: static grid coarsening enabled in direction(s) ";
    } else if (haveGridCoarsening == GridCoarsening::dynamic ) {
      idfx::cout << "Grid: dynamic grid coarsening enabled in direction(s) ";
    } else {
      IDEFIX_ERROR("Unknown grid coarsening");
    }
    for(int i = 0 ; i < 3 ; i++) {
      if(coarseningDirection[i]) {
        idfx::cout << "X" << i+1 << " ";
      }
    }
    idfx::cout << std::endl;
  }
}
