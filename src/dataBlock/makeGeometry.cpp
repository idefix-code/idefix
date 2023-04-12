// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include "../idefix.hpp"
#include "dataBlock.hpp"

/**
 * Compute the geometrical terms for the grid.
 * This function computes all of the geometrical terms
 * (cell volumes, faces surfaces, trigonometric functions)
*/
void DataBlock::MakeGeometry() {
  idfx::pushRegion("DataBlock::MakeGeometry()");


  // Initialize grid coarsening if needed
  if(mygrid->haveGridCoarsening != GridCoarsening::disabled) {
    this->haveGridCoarsening = mygrid->haveGridCoarsening;
    this->coarseningDirection = mygrid->coarseningDirection;

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

  // Compute Volumes
  IdefixArray3D<real> dV  = this->dV;
  IdefixArray1D<real> dx1 = this->dx[IDIR];
  IdefixArray1D<real> dx2 = this->dx[JDIR];
  IdefixArray1D<real> dx3 = this->dx[KDIR];
  IdefixArray1D<real> x1  = this->x[IDIR];
  IdefixArray1D<real> x2  = this->x[JDIR];
  IdefixArray1D<real> x3  = this->x[KDIR];
  IdefixArray1D<real> x1p = this->xr[IDIR];
  IdefixArray1D<real> x2p = this->xr[JDIR];
  IdefixArray1D<real> x3p = this->xr[KDIR];
  IdefixArray1D<real> x1m = this->xl[IDIR];
  IdefixArray1D<real> x2m = this->xl[JDIR];
  IdefixArray1D<real> x3m = this->xl[KDIR];
  IdefixArray1D<real> rt  = this->rt;
  IdefixArray1D<real> sinx2m  = this->sinx2m;
  IdefixArray1D<real> tanx2m  = this->tanx2m;
  IdefixArray1D<real> sinx2   = this->sinx2;
  IdefixArray1D<real> tanx2   = this->tanx2;
  IdefixArray1D<real> dmu = this->dmu;

  idefix_for("Volumes",0,this->np_tot[KDIR],0,this->np_tot[JDIR],0,this->np_tot[IDIR],
    KOKKOS_LAMBDA (int k, int j, int i) {
#if GEOMETRY == CARTESIAN
      dV(k,j,i) = D_EXPAND(dx1(i), *dx2(j), *dx3(k));   // = dx*dy*dz

#elif GEOMETRY == CYLINDRICAL
      real dVr = FABS(x1p(i)*x1p(i) - x1m(i)*x1m(i))/2.0;
      dV(k,j,i) = D_EXPAND( dVr, *dx2(j), *ONE_F);
             // = |r|*dr*dz  (more accurately (x1p**2-x1m**2)/2*dphi*dz)

#elif GEOMETRY == POLAR
      real dVr = FABS(x1p(i)*x1p(i) - x1m(i)*x1m(i))/2.0;
      dV(k,j,i) = D_EXPAND( dVr, *dx2(j), *dx3(k));    // = |r|*dr*dphi*dz

#elif GEOMETRY == SPHERICAL
      real dVr = FABS(x1p(i)*x1p(i)*x1p(i) - x1m(i)*x1m(i)*x1m(i))/3.0;
      real dmu = FABS(cos(x2m(j)) - cos(x2p(j)));
      dV(k,j,i) = D_EXPAND( dVr, *dmu, *dx3(k));
#endif
    }
  );

  // Compute Geometrical cell centers
  IdefixArray1D<real> x1gc = this->xgc[IDIR];
  IdefixArray1D<real> x2gc = this->xgc[JDIR];
  IdefixArray1D<real> x3gc = this->xgc[KDIR];

  // X1 direction
  idefix_for("GeometricalCentersX1",0,np_tot[IDIR],
    KOKKOS_LAMBDA (int i) {
#if GEOMETRY == CARTESIAN
      x1gc(i) = x1(i);

#elif GEOMETRY == CYLINDRICAL || GEOMETRY == POLAR
      x1gc(i) = x1(i) + dx1(i)*dx1(i)/(12.0*x1(i));

#elif GEOMETRY == SPHERICAL
      x1gc(i) = x1(i) + 2.0*x1(i)*dx1(i)*dx1(i)
                / (12.0*x1(i)*x1(i) + dx1(i)*dx1(i));
      rt(i)   = (x1p(i)*x1p(i)*x1p(i) - x1m(i)*x1m(i)*x1m(i))
                / (x1p(i)*x1p(i)-x1m(i)*x1m(i)) / 1.5;
#endif
    }
  );

  // X2
  idefix_for("GeometricalCentersX2",0,np_tot[JDIR],
    KOKKOS_LAMBDA (int j) {
#if GEOMETRY != SPHERICAL
      x2gc(j) = x2(j);
#else
      real xL = x2m(j);
      real xR = x2p(j);
      x2gc(j)  = (sin(xR) - sin(xL)+ xL*cos(xL) - xR*cos(xR)) / (cos(xL)-cos(xR));
      sinx2m(j) = sin(xL);
      tanx2m(j) = tan(xL);
      sinx2(j) = sin(x2(j));
      tanx2(j) = tan(x2(j));
      dmu(j) = FABS(cos(xL)-cos(xR));
#endif
    }
  );

  // X3
  idefix_for("GeometricalCentersX3",0,np_tot[KDIR],
    KOKKOS_LAMBDA (int k) {
      x3gc(k) = x3(k);
    }
  );

  // Compute Areas
  IdefixArray3D<real> Ax1 = this->A[IDIR];
  IdefixArray3D<real> Ax2 = this->A[JDIR];
  IdefixArray3D<real> Ax3 = this->A[KDIR];

  // X1 direction
  int end = this->np_tot[IDIR];
  idefix_for("AreaX1",0,this->np_tot[KDIR],0,this->np_tot[JDIR],0,this->np_tot[IDIR]+IOFFSET,
    KOKKOS_LAMBDA (int k, int j, int i) {
#if GEOMETRY == CARTESIAN
      Ax1(k,j,i) = D_EXPAND(1.0, *dx2(j), *dx3(k));         // = dy*dz

#elif GEOMETRY == CYLINDRICAL
      if(i == end) {
        Ax1(k,j,i) = D_EXPAND(FABS(x1p(i-1)), *dx2(j), *ONE_F);
      } else {
        Ax1(k,j,i) = D_EXPAND(FABS(x1m(i)), *dx2(j), *ONE_F); // r*dz
      }

#elif GEOMETRY == POLAR
      if(i == end) {
        Ax1(k,j,i) = D_EXPAND(FABS(x1p(i-1)), *dx2(j), *dx3(k));
      } else {
        Ax1(k,j,i) = D_EXPAND(FABS(x1m(i)), *dx2(j), *dx3(k)); // r*dphi*dz
      }

#elif GEOMETRY == SPHERICAL
      real dmu = FABS(cos(x2m(j)) - cos(x2p(j)));
      if(i == end) {
        Ax1(k,j,i) = D_EXPAND(x1p(i-1)*x1p(i-1), *dmu, *dx3(k));
      } else {
        Ax1(k,j,i) = D_EXPAND(x1m(i)*x1m(i), *dmu, *dx3(k)); // r^2*dmu*dphi
      }
#endif
    }
  );

  // X2 direction
  end = this->np_tot[JDIR];
  idefix_for("AreaX2",0,this->np_tot[KDIR],0,this->np_tot[JDIR]+JOFFSET,0,this->np_tot[IDIR],
    KOKKOS_LAMBDA (int k, int j, int i) {
#if GEOMETRY == CARTESIAN
      Ax2(k,j,i) = D_EXPAND(dx1(i), *ONE_F, *dx3(k));        // = dx*dz

#elif GEOMETRY == CYLINDRICAL
      Ax2(k,j,i) = D_EXPAND(FABS(x1(i)), *dx1(i), *ONE_F);   // = r*dr

#elif GEOMETRY == POLAR
      Ax2(k,j,i) = D_EXPAND(dx1(i), *ONE_F, *dx3(k));        // = dr*dz

#elif GEOMETRY == SPHERICAL
      if (j == end) {
        Ax2(k,j,i) = D_EXPAND(x1(i)*dx1(i), *FABS(sin(x2p(j-1))), *dx3(k));
      } else {
        Ax2(k,j,i) = D_EXPAND(x1(i)*dx1(i), *FABS(sin(x2m(j))), *dx3(k)); // = r*dr*sin(thp)*dphi
      }
#endif
    }
  );

  // X3 direction
  end = this->np_tot[KDIR];
  idefix_for("AreaX3",0,this->np_tot[KDIR]+KOFFSET,0,this->np_tot[JDIR],0,this->np_tot[IDIR],
    KOKKOS_LAMBDA (int k, int j, int i) {
  #if GEOMETRY == CARTESIAN
      Ax3(k,j,i) = D_EXPAND(dx1(i), *dx2(j), *ONE_F);          // = dx*dy

  #elif GEOMETRY == CYLINDRICAL
      Ax3(k,j,i) = ONE_F;   // No 3rd direction in cylindrical coords

  #elif GEOMETRY == POLAR
      Ax3(k,j,i) = D_EXPAND(x1(i)*dx1(i), *dx2(j), *ONE_F);   // = r*dr*dphi

  #elif GEOMETRY == SPHERICAL
      Ax3(k,j,i) = D_EXPAND(x1(i)*dx1(i), *dx2(j), *ONE_F);   // = r*dr*dth
  #endif
    }
  );

  idfx::popRegion();
}

void DataBlock::ExtractSubdomain() {
  idfx::pushRegion("DataBlock::ExtractSubdomain");
  // Copy the relevant part of the coordinate system to the datablock
  for(int dir = 0 ; dir < 3 ; dir++) {
    int offset=gbeg[dir]-beg[dir];

    IdefixArray1D<real> x_input = mygrid->x[dir];
    IdefixArray1D<real> x_output= x[dir];
    IdefixArray1D<real> xr_input = mygrid->xr[dir];
    IdefixArray1D<real> xr_output= xr[dir];
    IdefixArray1D<real> xl_input = mygrid->xl[dir];
    IdefixArray1D<real> xl_output= xl[dir];
    IdefixArray1D<real> dx_input = mygrid->dx[dir];
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

  idfx::popRegion();
}
