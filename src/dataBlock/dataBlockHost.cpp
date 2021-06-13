// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2021 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include "idefix.hpp"
#include "dataBlockHost.hpp"

DataBlockHost::DataBlockHost() {
  // Do nothing
}

DataBlockHost::DataBlockHost(DataBlock& datain) {
  idfx::pushRegion("DataBlockHost::DataBlockHost(DataBlock)");

  // copy the dataBlock object for later use
  this->data=&datain;

  // By default, no current
  this->haveCurrent = false;

  // Create mirrors (should be mirror_view)
  for(int dir = 0 ; dir < 3 ; dir++) {
    x[dir] = Kokkos::create_mirror_view(data->x[dir]);
    xr[dir] = Kokkos::create_mirror_view(data->xr[dir]);
    xl[dir] = Kokkos::create_mirror_view(data->xl[dir]);
    dx[dir] = Kokkos::create_mirror_view(data->dx[dir]);
    A[dir] = Kokkos::create_mirror_view(data->A[dir]);
    np_tot[dir] = data->np_tot[dir];
    np_int[dir] = data->np_int[dir];
    nghost[dir] = data->nghost[dir];

    xbeg[dir] = data->xbeg[dir];
    xend[dir] = data->xend[dir];
    beg[dir] = data->beg[dir];
    end[dir] = data->end[dir];
    gbeg[dir] = data->gbeg[dir];
    gend[dir] = data->gend[dir];

    // TO BE COMPLETED...
  }

  dV = Kokkos::create_mirror_view(data->dV);
  Vc = Kokkos::create_mirror_view(data->hydro.Vc);
  Uc = Kokkos::create_mirror_view(data->hydro.Uc);
  InvDt = Kokkos::create_mirror_view(data->hydro.InvDt);

#if MHD == YES
  Vs = Kokkos::create_mirror_view(data->hydro.Vs);
  this->haveCurrent = data->hydro.haveCurrent;
  if(data->hydro.haveCurrent) {
    J = Kokkos::create_mirror_view(data->hydro.J);
  }

  D_EXPAND( Ex3 = Kokkos::create_mirror_view(data->hydro.emf.ez);  ,
                                                             ,
            Ex1 = Kokkos::create_mirror_view(data->hydro.emf.ex);
            Ex2 = Kokkos::create_mirror_view(data->hydro.emf.ey);  )
#endif

  // Store the grid informations from the dataBlock
  for(int dir = 0 ; dir < 3 ; dir++) {
    Kokkos::deep_copy(x[dir],data->x[dir]);
    Kokkos::deep_copy(xr[dir],data->xr[dir]);
    Kokkos::deep_copy(xl[dir],data->xl[dir]);
    Kokkos::deep_copy(dx[dir],data->dx[dir]);
    Kokkos::deep_copy(A[dir],data->A[dir]);
  }

  Kokkos::deep_copy(dV,data->dV);

  idfx::popRegion();
}

// Synchronisation routines of Data (*Only*)
void DataBlockHost::SyncToDevice() {
  idfx::pushRegion("DataBlockHost::SyncToDevice()");

  Kokkos::deep_copy(data->hydro.Vc,Vc);
  Kokkos::deep_copy(data->hydro.InvDt,InvDt);

#if MHD == YES
  Kokkos::deep_copy(data->hydro.Vs,Vs);
  if(this->haveCurrent && data->hydro.haveCurrent) Kokkos::deep_copy(data->hydro.J,J);
  D_EXPAND( Kokkos::deep_copy(data->hydro.emf.ez,Ex3);  ,
                                                  ,
            Kokkos::deep_copy(data->hydro.emf.ex,Ex1);
            Kokkos::deep_copy(data->hydro.emf.ey,Ex2);  )
#endif

  Kokkos::deep_copy(data->hydro.Uc,Uc);

  idfx::popRegion();
}

void DataBlockHost::SyncFromDevice() {
  idfx::pushRegion("DataBlockHost::SyncFromDevice()");
  Kokkos::deep_copy(Vc,data->hydro.Vc);
  Kokkos::deep_copy(InvDt,data->hydro.InvDt);

#if MHD == YES
  Kokkos::deep_copy(Vs,data->hydro.Vs);
  if(this->haveCurrent && data->hydro.haveCurrent) Kokkos::deep_copy(J,data->hydro.J);

  D_EXPAND( Kokkos::deep_copy(Ex3,data->hydro.emf.ez);  ,
                                                  ,
            Kokkos::deep_copy(Ex1,data->hydro.emf.ex);
            Kokkos::deep_copy(Ex2,data->hydro.emf.ey);  )
#endif

  Kokkos::deep_copy(Uc,data->hydro.Uc);

  idfx::popRegion();
}

void DataBlockHost::MakeVsFromAmag(IdefixHostArray4D<real> &Ain) {
  IdefixHostArray1D<real> dx1 = this->dx[IDIR];
  IdefixHostArray1D<real> dx2 = this->dx[JDIR];
  IdefixHostArray1D<real> dx3 = this->dx[KDIR];

  IdefixHostArray1D<real> x1m = this->xl[IDIR];
  IdefixHostArray1D<real> x2m = this->xl[JDIR];
  IdefixHostArray1D<real> x1 = this->x[IDIR];
  IdefixHostArray1D<real> x2 = this->x[JDIR];

#if MHD == YES

  for(int k = data->beg[KDIR] ; k < data->end[KDIR] + KOFFSET ; k++) {
    for(int j = data->beg[JDIR] ; j < data->end[JDIR] + JOFFSET ; j++) {
      for(int i = data->beg[IDIR] ; i < data->end[IDIR] + IOFFSET; i++) {
  #if GEOMETRY == CARTESIAN
        Vs(BX1s,k,j,i) = D_EXPAND( ZERO_F                                               ,
                                   + 1/dx2(j) * (Ain(KDIR,k,j+1,i) - Ain(KDIR,k,j,i) )  ,
                                   - 1/dx3(k) * (Ain(JDIR,k+1,j,i) - Ain(JDIR,k,j,i) )  );
    #if DIMENSIONS >= 2
        Vs(BX2s,k,j,i) = D_EXPAND( - 1/dx1(i) * (Ain(KDIR,k,j,i+1) - Ain(KDIR,k,j,i) )  ,
                                                                                        ,
                                   + 1/dx3(k) * (Ain(IDIR,k+1,j,i) - Ain(IDIR,k,j,i) )  );
    #endif
    #if DIMENSIONS == 3
        Vs(BX3s,k,j,i) = 1/dx1(i) * (Ain(JDIR,k,j,i+1) - Ain(JDIR,k,j,i) )
                         - 1/dx2(j) * (Ain(IDIR,k,j+1,i) - Ain(IDIR,k,j,i) );
    #endif
  #endif
  #if GEOMETRY == CYLINDRICAL
        IDEFIX_ERROR("Not yet defined");
  #endif
  #if GEOMETRY == POLAR
        Vs(BX1s,k,j,i) = D_EXPAND( ZERO_F                                                       ,
                                   + 1/(x1m(i)*dx2(j)) * (Ain(KDIR,k,j+1,i) - Ain(KDIR,k,j,i) ) ,
                                   - 1/dx3(k) * (Ain(JDIR,k+1,j,i) - Ain(JDIR,k,j,i) )          );

        Vs(BX2s,k,j,i) = D_EXPAND( - 1/dx1(i) * (Ain(KDIR,k,j,i+1) - Ain(KDIR,k,j,i) )  ,
                                                                                        ,
                                   + 1/dx3(k) * (Ain(IDIR,k+1,j,i) - Ain(IDIR,k,j,i) )  );

    #if DIMENSIONS == 3
        Vs(BX3s,k,j,i) = 1/(x1(i)*dx1(i)) * (x1m(i+1)*Ain(JDIR,k,j,i+1) - x1m(i)*Ain(JDIR,k,j,i) )
                         - 1/(x1(i)*dx2(j)) * (Ain(IDIR,k,j+1,i) - Ain(IDIR,k,j,i) );
    #endif
  #endif
  #if GEOMETRY == SPHERICAL
        Vs(BX1s,k,j,i) = D_EXPAND( ZERO_F                                                ,
                                   + 1/(x1m(i)*(cos(x2m(j))
                                   - cos(x2m(j+1)))) * (sin(x2m(j+1))*Ain(KDIR,k,j+1,i)
                                   - sin(x2m(j))*Ain(KDIR,k,j,i) )                       ,
                                   - 1/(x1m(i)*sin(x2(j))*dx3(k)) * (Ain(JDIR,k+1,j,i)
                                   - Ain(JDIR,k,j,i) )                                   );

        real Ax2m = fabs(sin(x2m(j)));
        // Regularisation along the axis
        if(FABS(Ax2m)<1e-12) Ax2m = ONE_F;
        Vs(BX2s,k,j,i) = D_EXPAND( - 1/(x1(i)*dx1(i)) * (x1m(i+1)*Ain(KDIR,k,j,i+1)
                                   - x1m(i)*Ain(KDIR,k,j,i) )                           ,
                                                                                        ,
                                   + 1/(x1m(i)*Ax2m*dx3(k)) * (Ain(IDIR,k+1,j,i)
                                   - Ain(IDIR,k,j,i) )                                  );

    #if DIMENSIONS == 3
        Vs(BX3s,k,j,i) = 1/(x1(i)*dx1(i)) * (x1m(i+1)*Ain(JDIR,k,j,i+1) - x1m(i)*Ain(JDIR,k,j,i) )
                        - 1/(x1(i)*dx2(j)) * (Ain(IDIR,k,j+1,i) - Ain(IDIR,k,j,i) );
    #endif
  #endif
      }
    }
  }
#else
  IDEFIX_ERROR("This function cannot be used without MHD enabled");
#endif
}
