// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include <string>
#include <vector>

#include "idefix.hpp"
#include "grid.hpp"
#include "gridHost.hpp"

GridHost::GridHost(Grid &grid) {
  idfx::pushRegion("GridHost::GridHost(Grid)");

  this->grid=&grid;
  nghost = grid.nghost;
  np_tot = grid.np_tot;
  np_int = grid.np_int;

  lbound = grid.lbound;
  rbound = grid.rbound;

  xbeg = grid.xbeg;
  xend = grid.xend;

  haveAxis = grid.haveAxis;

  isRegularCartesian = grid.isRegularCartesian;

  // Create mirrors on host
  for(int dir = 0 ; dir < 3 ; dir++) {
    x[dir] = Kokkos::create_mirror_view(grid.x[dir]);
    xr[dir] = Kokkos::create_mirror_view(grid.xr[dir]);
    xl[dir] = Kokkos::create_mirror_view(grid.xl[dir]);
    dx[dir] = Kokkos::create_mirror_view(grid.dx[dir]);
  }

  idfx::popRegion();
}

void GridHost::MakeGrid(Input &input) {
  idfx::pushRegion("GridHost::MakeGrid");

  real xstart[3];
  real xend[3];
  // Create the grid

  // Get grid parameters from input file, block [Grid]
  for(int dir = 0 ; dir < 3 ; dir++) {
    // These are extra dimensions that are not being used.

    if(dir<DIMENSIONS) {
      std::string label = std::string("X")+std::to_string(dir+1)+std::string("-grid");
      int numPatch = input.Get<int>("Grid",label,0);

      xstart[dir] = input.Get<real>("Grid",label,1);
      xend[dir] = input.Get<real>("Grid",label,4+(numPatch-1)*3);

      this->xbeg[dir] = xstart[dir];
      this->xend[dir] = xend[dir];
      // First, we fill cells for any non strecthed patch
      // Loop on all the patches
      int idxstart = nghost[dir];
      for(int patch = 0 ; patch < numPatch ; patch++) {
        std::string patchType = input.Get<std::string>("Grid",label,3+patch*3);
        real patchStart = input.Get<real>("Grid",label,1+patch*3);
        real patchEnd = input.Get<real>("Grid",label,4+patch*3);
        int patchSize = input.Get<int>("Grid",label,2+patch*3);

        // If this is the first or last patch, also define ghost cells
        int ghostStart = 0;
        int ghostEnd = 0;
        if(patch == 0) ghostStart = nghost[dir];
        if(patch == numPatch-1) ghostEnd = nghost[dir];
        // Define the grid depending on patch type
        if(patchType.compare("u")==0) {
          // Uniform patch
          for(int i = 0 - ghostStart ; i < patchSize + ghostEnd ; i++) {
            dx[dir](i+idxstart) = (patchEnd-patchStart)/(patchSize);
            x[dir](i+idxstart)=patchStart + (i+HALF_F)*dx[dir](i+idxstart);
            xl[dir](i+idxstart)=patchStart  + i*dx[dir](i+idxstart);
            xr[dir](i+idxstart)=patchStart  + (i+1)*dx[dir](i+idxstart);
          }
        } else if(patchType.compare("l")==0) {
          // log-increasing patch
          isRegularCartesian = false;
          double alpha = (patchEnd + fabs(patchStart) - patchStart)/fabs(patchStart);

          for(int i = 0 - ghostStart ; i < patchSize + ghostEnd ; i++) {
            xl[dir](i+idxstart) = patchStart * pow(alpha,
                                    static_cast<double>(i) / (static_cast<double>(patchSize)));
            xr[dir](i+idxstart) = patchStart * pow(alpha,
                                    static_cast<double>(i+1) / (static_cast<double>(patchSize)));
            dx[dir](i+idxstart) = xr[dir](i+idxstart) - xl[dir](i+idxstart);
            x[dir](i+idxstart)= 0.5*(xr[dir](i+idxstart) + xl[dir](i+idxstart));
          }
        } else if((patchType.compare("s+"))&&(patchType.compare("s-"))) {
          std::stringstream msg;
          msg << "GridHost::MakeGrid: Unknown grid type :" << patchType << std::endl;
          IDEFIX_ERROR(msg);
        }

        // Increment offset
        idxstart += patchSize;
      }

      // Second, we fill cells for stretched patch
      // Loop on all the patches
      idxstart = nghost[dir];
      for(int patch = 0 ; patch < numPatch ; patch++) {
        std::string patchType = input.Get<std::string>("Grid",label,3+patch*3);
        real patchStart = input.Get<real>("Grid",label,1+patch*3);
        real patchEnd = input.Get<real>("Grid",label,4+patch*3);
        int patchSize = input.Get<int>("Grid",label,2+patch*3);

        // If this is the first or last patch, also define ghost cells
        int ghostStart = 0;
        int ghostEnd = 0;
        if(patch == 0) ghostStart = nghost[dir];
        if(patch == numPatch-1) ghostEnd = nghost[dir];
        // Define the grid depending on patch type
        if((patchType.compare("s+")==0)||(patchType.compare("s-")==0)) {
          isRegularCartesian = false;
          // Stretched grid
          // - means we take the initial dx on the left side, + on the right side
          // refPatch corresponds to the patch from which we compute the initial dx
          // of the stretched grid

          int refPatch=patch;
          if(patchType.compare("s+")==0) {
            refPatch=patch+1;
          } else {
            refPatch=patch-1;
          }
          // Sanity check
          // Check that the reference patch actually exist
          if(refPatch<0 || refPatch >= numPatch) {
            IDEFIX_ERROR("You're attempting to construct a stretched patch "
                         "from a non-existent patch");
          }
          // Check that we build from a uniform patch or after a logarithmic one
          std::string refPatchType = input.Get<std::string>("Grid",label,3+3*refPatch);
          if(refPatchType.compare("u")) {
            if((refPatchType.compare("l")) || (refPatch==patch+1)) {
              IDEFIX_ERROR("You can only construct a stretched patch "
                           "from a uniform grid "
                           "or AFTER a logarithmic grid");
            }
          }

          // Retrieve dx from the reference patch
          int offset = (refPatch==patch+1) ? patchSize : -1;
          double delta = dx[dir](idxstart+offset);
          double logdelta = log((patchEnd-patchStart)/delta);

          // Check that it is possible to make a stretch grid (bug report #28)
          if(std::fabs((patchEnd-patchStart)/patchSize - delta) < 1e-10) {
            IDEFIX_ERROR("A Stretch grid can be defined only if the stretched domain has a mean\n"
                         "spacing different from the reference uniform grid.\n"
                         "Try changing the number of points in your stretched grid.");
          }
          // Next we have to compute the stretch factor q. Let's start with a guess
          double q=1.05;
          // Use Newton method
          for(int iter=0; iter <= 50; iter++) {
            double f = log((pow(q,patchSize+1)-q)/(q-1))-logdelta;
            double fp = ((patchSize+1)*pow(q,patchSize)-1)/(pow(q,patchSize+1)-q)-1/(q-1);
            double dq = f/fp;
            // advance the guess
            q = q - dq;
            // Check whether we have converged
            if(fabs(dq)<1e-14*q) break;
            if(iter==50) IDEFIX_ERROR("Failed to create the stretched grid");
          }

          // once we know q, we can make the grid
          if(patchType.compare("s-")==0) {
            for(int i = 0 - ghostStart ; i < patchSize + ghostEnd ; i++) {
              xl[dir](i+idxstart) = patchStart + q*(pow(q,i)-1)/(q-1)*delta;
              xr[dir](i+idxstart) = patchStart + q*(pow(q,i+1)-1)/(q-1)*delta;
              dx[dir](i+idxstart) = pow(q,i+1)*delta;
              x[dir](i+idxstart)= 0.5*(xr[dir](i+idxstart) + xl[dir](i+idxstart));
            }
          } else {
            for(int i = 0 - ghostStart ; i < patchSize + ghostEnd ; i++) {
              xl[dir](i+idxstart) = patchEnd - q*(pow(q,patchSize-i)-1)/(q-1)*delta;
              xr[dir](i+idxstart) = patchEnd - q*(pow(q,patchSize-i-1)-1)/(q-1)*delta;
              dx[dir](i+idxstart) = pow(q,patchSize-i)*delta;
              x[dir](i+idxstart)= 0.5*(xr[dir](i+idxstart) + xl[dir](i+idxstart));
            }
          }
        }

        // Increment offset
        idxstart += patchSize;
      }
    } else {
      // dir >= DIMENSIONS/ Init simple uniform grid
      for(int i = 0 ; i < np_tot[dir] ; i++) {
        // Initialize to default values
        #if GEOMETRY != SPHERICAL
          xstart[dir] = -0.5;
          xend[dir] = 0.5;
        #elif GEOMETRY == SPHERICAL
          if(dir == JDIR) {
            xstart[dir] = 0;
            xend[dir] = M_PI;
          }
          if(dir == KDIR) {
            xstart[dir] = -0.5;
            xend[dir] = 0.5;
          }
        #endif
        this->xbeg[dir] = xstart[dir];
        this->xend[dir] = xend[dir];
        dx[dir](i) = xend[dir] - xstart[dir];
        x[dir](i)=0.5*(xstart[dir]+xend[dir]);
        xl[dir](i)=xstart[dir];
        xr[dir](i)=xend[dir];
      }
    }
  }

  // Check that axis treatment is compatible with the domain
if(haveAxis) {
  #if GEOMETRY != SPHERICAL
    IDEFIX_ERROR("Axis boundaries only compatible with Spherical boundary conditions");
  #endif
  #if DIMENSIONS < 2
    IDEFIX_ERROR("Axis Boundaries requires at least two dimenions");
  #endif
  #ifdef SINGLE_PRECISION
    const real smallNumber = 1e-5;
  #else
    const real smallNumber = 1e-10;
  #endif
  if((fabs(xbeg[JDIR])>smallNumber) && (lbound[JDIR] == axis)) {
    IDEFIX_ERROR("Axis Boundaries requires your X2 domain to start at exactly x2=0.0");
  }
  if((fabs(xend[JDIR]-M_PI)>smallNumber) && (rbound[JDIR] == axis )) {
    IDEFIX_ERROR("Axis Boundaries requires your X2 domain to end at exactly x2=Pi");
  }
  // Enforce symmetry of theta grid spacing when we cross the axis
  if(lbound[JDIR] == axis) {
    int jref = nghost[JDIR];
    for(int j = jref - 1 ; j>=0 ; j -- ) {
      dx[JDIR](j) = dx[JDIR](2*jref - j - 1);
      xl[JDIR](j) = xl[JDIR](j+1)- dx[JDIR](j);
      xr[JDIR](j) = xl[JDIR](j+1);
    }
  }
  if(rbound[JDIR] == axis) {
    int jref = nghost[JDIR]+np_int[JDIR] - 1;
    for(int j = jref + 1 ; j<np_tot[JDIR] ; j++ ) {
      dx[JDIR](j) = dx[JDIR](2*jref - j + 1);
      xr[JDIR](j) = xr[JDIR](j-1) + dx[JDIR](j);
      xl[JDIR](j) = xr[JDIR](j-1);
    }
  }
}

  idfx::popRegion();
}

void GridHost::SyncFromDevice() {
  idfx::pushRegion("GridHost::SyncFromDevice");

  for(int dir = 0 ; dir < 3 ; dir++) {
    Kokkos::deep_copy(x[dir],grid->x[dir]);
    Kokkos::deep_copy(xr[dir],grid->xr[dir]);
    Kokkos::deep_copy(xl[dir],grid->xl[dir]);
    Kokkos::deep_copy(dx[dir],grid->dx[dir]);
  }

  xbeg = grid->xbeg;
  xend = grid->xend;
  isRegularCartesian = grid->isRegularCartesian;

  idfx::popRegion();
}

void GridHost::SyncToDevice() {
  idfx::pushRegion("GridHost::SyncToDevice");

  // Sync with the device
  for(int dir = 0 ; dir < 3 ; dir++) {
    Kokkos::deep_copy(grid->x[dir],x[dir]);
    Kokkos::deep_copy(grid->xr[dir],xr[dir]);
    Kokkos::deep_copy(grid->xl[dir],xl[dir]);
    Kokkos::deep_copy(grid->dx[dir],dx[dir]);
  }

  grid->xbeg = xbeg;
  grid->xend = xend;
  grid->isRegularCartesian = isRegularCartesian;

  idfx::popRegion();
}
