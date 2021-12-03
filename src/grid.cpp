// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2021 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include <string>

#include "idefix.hpp"
#include "gridHost.hpp"
#include "grid.hpp"

Grid::Grid() {
  // Do nothing
}

Grid::Grid(Input &input) {
  idfx::pushRegion("Grid::Grid(Input)");

  idfx::cout << "Grid: allocating Grid." << std::endl;

  // Get grid size from input file, block [Grid]
  int npoints[3];
  for(int dir = 0 ; dir < 3 ; dir++) {
    npoints[dir] = 1;
    nghost[dir] = 0;
    std::string label = std::string("X")+std::to_string(dir+1)+std::string("-grid");
    int numPatch = input.GetInt("Grid",label,0);

    if(dir<DIMENSIONS) {
      #if ORDER < 4
        nghost[dir] = 2;
      #else
        nghost[dir] = 3;
      #endif
      npoints[dir] = 0;
      for(int patch = 0; patch < numPatch ; patch++) {
        npoints[dir] += input.GetInt("Grid",label,2+3*patch );
      }
    }
  }


  for(int dir = 0 ; dir < 3 ; dir++) {
    np_tot[dir] = npoints[dir] + 2*nghost[dir];
    np_int[dir] = npoints[dir];

    std::string label = std::string("X")+std::to_string(dir+1)+std::string("-beg");
    std::string boundary = input.GetString("Boundary",label,0);

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
    boundary = input.GetString("Boundary",label,0);
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

  for(int dir = 0 ; dir < 3 ; dir++) {
    x[dir] = IdefixArray1D<real>("Grid_x",np_tot[dir]);
    xr[dir] = IdefixArray1D<real>("Grid_xr",np_tot[dir]);
    xl[dir] = IdefixArray1D<real>("Grid_xl",np_tot[dir]);
    dx[dir] = IdefixArray1D<real>("Grid_dx",np_tot[dir]);
  }

  // Allocate proc structure (by default: one proc in each direction, size one)
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
      // No command line decomposition, see if auto-decomposition is possible
      // (only when nproc and dimensions are powers of 2)
      bool autoDecomposition=true;
      if(!isPow2(idfx::psize)) autoDecomposition=false;
      for(int dir = 0; dir < 3; dir++) {
        if(np_int[dir]>1) {
          if(!isPow2(np_int[dir])) autoDecomposition=false;
        }
      }
      if(autoDecomposition)
        makeDomainDecomposition();
      else
        IDEFIX_ERROR("-dec option is mandatory when nproc or nx1, nx2, nx3 are not powers of 2.");
    } else {
      // Manual domain decomposition (with -dec option)
      int ntot=1;
      for(int dir=0 ; dir < DIMENSIONS; dir++) {
        nproc[dir] = input.GetInt("CommandLine","dec",dir);
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
  MPI_Cart_create(MPI_COMM_WORLD, 3, nproc, period, 0, &CartComm);
  MPI_Cart_coords(CartComm, idfx::prank, 3, xproc);

  MPI_Barrier(MPI_COMM_WORLD);
  idfx::cout << "Grid::Grid: Current MPI proc coordinates (";

  for(int dir = 0; dir < 3; dir++) {
    idfx::cout << xproc[dir];
    if(dir < 2) idfx::cout << ", ";
  }
  idfx::cout << ")" << std::endl;

#endif

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

    // If we have the axis, we should not decompose in phi
    // This is a bit too conservative since if we're not doing full two pi, domain decomposition
    // in phi is allowed. However, the grid bounds are only initialised later, so we don't have
    // this information yet. Hence, automatic domain decomposition is conservative in that case.
    // At least we know it works, even though it's probably sub-optimal in some cases.
    if(haveAxis) ndir = 1;

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

  idfx::cout << "Grid::makeDomainDecomposition: grid is (";
  for(int dir = 0 ; dir < DIMENSIONS ; dir++) {
    idfx::cout << " " << nproc[dir] << " ";
  }
  idfx::cout << ")" << std::endl;
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
