// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2022 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include "hydroboundary.hpp"
#include "hydro.hpp"
#include "dataBlock.hpp"
#include "boundaryloop.hpp"


void HydroBoundary::Init(Input & input, Grid &grid, Hydro* hydro) {
  idfx::pushRegion("HydroBoundary::Init");
  this->hydro = hydro;
  this->data = hydro->data;

  // This should be required only when shearing box is on
  sBArray = IdefixArray4D<real>("ShearingBoxArray",
                                NVAR,
                                data->np_tot[KDIR]+1,
                                data->np_tot[JDIR]+1,
                                data->nghost[IDIR]);

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

  std::vector<int> mapVars;

  // Init the list of variables we will exchange in MPI routines
  int ntarget = 0;
  for(int n = 0 ; n < mapNVars ; n++) {
    mapVars.push_back(ntarget);
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
  #if MHD == YES
    mpi.Init(data->mygrid, mapVars, data->nghost.data(), data->np_int.data(), true);
  #else
    mpi.Init(data->mygrid, mapVars, data->nghost.data(), data->np_int.data());
  #endif
#endif // MPI
  idfx::popRegion();
}

void HydroBoundary::EnrollFluxBoundary(UserDefBoundaryFunc myFunc) {
  this->haveFluxBoundary = true;
  this->fluxBoundaryFunc = myFunc;
}

void HydroBoundary::EnforceFluxBoundaries(int dir) {
  idfx::pushRegion("HydroBoundary::EnforceFluxBoundaries");
  if(haveFluxBoundary) {
    if(data->lbound[dir] != internal) {
      fluxBoundaryFunc(*data, dir, left, data->t);
    }
    if(data->rbound[dir] != internal) {
      fluxBoundaryFunc(*data, dir, right, data->t);
    }
  } else {
    IDEFIX_ERROR("Cannot enforce flux boundary conditions without enrolling a specific function");
  }
  idfx::popRegion();
}

void HydroBoundary::SetBoundaries(real t) {
  idfx::pushRegion("HydroBoundary::SetBoundaries");
  // set internal boundary conditions
  if(haveInternalBoundary) internalBoundaryFunc(*data, t);
  for(int dir=0 ; dir < DIMENSIONS ; dir++ ) {
      // MPI Exchange data when needed
    #ifdef WITH_MPI
    if(data->mygrid->nproc[dir]>1) {
      switch(dir) {
        case 0:
          mpi.ExchangeX1(hydro->Vc, hydro->Vs);
          break;
        case 1:
          mpi.ExchangeX2(hydro->Vc, hydro->Vs);
          break;
        case 2:
          mpi.ExchangeX3(hydro->Vc, hydro->Vs);
          break;
      }
    }
    #endif
    EnforceBoundaryDir(t, dir);
    #if MHD == YES
      // Reconstruct the normal field component when using CT
      ReconstructNormalField(dir);
    #endif
  } // Loop on dimension ends

#if MHD == YES
  // Remake the cell-centered field.
  ReconstructVcField(hydro->Vc);
#endif
  idfx::popRegion();
}


// Enforce boundary conditions by writing into ghost zones
void HydroBoundary::EnforceBoundaryDir(real t, int dir) {
  idfx::pushRegion("Hydro::EnforceBoundaryDir");

  // left boundary

  switch(data->lbound[dir]) {
    case internal:
      // internal is used for MPI-enforced boundary conditions. Nothing to be done here.
      break;

    case periodic:
      if(data->mygrid->nproc[dir] > 1) break; // Periodicity already enforced by MPI calls
      EnforcePeriodic(dir,left);
      break;

    case reflective:
      EnforceReflective(dir,left);
      break;

    case outflow:
      EnforceOutflow(dir,left);
      break;

    case shearingbox:
      EnforceShearingBox(t,dir,left);
      break;
    case axis:
      hydro->myAxis.EnforceAxisBoundary(left);
      break;
    case userdef:
      if(this->haveUserDefBoundary)
        this->userDefBoundaryFunc(*data, dir, left, t);
      else
        IDEFIX_ERROR("No function has been enrolled to define your own boundary conditions");
      break;

    default:
      std::stringstream msg ("Boundary condition type is not yet implemented");
      IDEFIX_ERROR(msg);
  }

  // right boundary

  switch(data->rbound[dir]) {
    case internal:
      // internal is used for MPI-enforced boundary conditions. Nothing to be done here.
      break;

    case periodic:
      if(data->mygrid->nproc[dir] > 1) break; // Periodicity already enforced by MPI calls
      EnforcePeriodic(dir,right);
      break;
    case reflective:
      EnforceReflective(dir,right);
      break;
    case outflow:
      EnforceOutflow(dir,right);
      break;
    case shearingbox:
      EnforceShearingBox(t,dir,right);
      break;
    case axis:
      hydro->myAxis.EnforceAxisBoundary(right);
      break;
    case userdef:
      if(this->haveUserDefBoundary)
        this->userDefBoundaryFunc(*data, dir, right, t);
      else
        IDEFIX_ERROR("No function has been enrolled to define your own boundary conditions");
      break;
    default:
      std::stringstream msg("Boundary condition type is not yet implemented");
      IDEFIX_ERROR(msg);
  }

  idfx::popRegion();
}


void HydroBoundary::ReconstructVcField(IdefixArray4D<real> &Vc) {
  idfx::pushRegion("Hydro::ReconstructVcField");

  IdefixArray4D<real> Vs=hydro->Vs;

  // Reconstruct cell average field when using CT
  idefix_for("ReconstructVcMagField",0,data->np_tot[KDIR],0,data->np_tot[JDIR],0,data->np_tot[IDIR],
    KOKKOS_LAMBDA (int k, int j, int i) {
      D_EXPAND( Vc(BX1,k,j,i) = HALF_F * (Vs(BX1s,k,j,i) + Vs(BX1s,k,j,i+1)) ;  ,
                Vc(BX2,k,j,i) = HALF_F * (Vs(BX2s,k,j,i) + Vs(BX2s,k,j+1,i)) ;  ,
                Vc(BX3,k,j,i) = HALF_F * (Vs(BX3s,k,j,i) + Vs(BX3s,k+1,j,i)) ;  )
    }
  );

  idfx::popRegion();
}


void HydroBoundary::ReconstructNormalField(int dir) {
  idfx::pushRegion("Hydro::ReconstructNormalField");

  // Reconstruct the field
  IdefixArray4D<real> Vs = hydro->Vs;
  // Coordinates
  IdefixArray1D<real> x1=data->x[IDIR];
  IdefixArray1D<real> x2=data->x[JDIR];
  IdefixArray1D<real> x3=data->x[KDIR];
  IdefixArray1D<real> dx1=data->dx[IDIR];
  IdefixArray1D<real> dx2=data->dx[JDIR];
  IdefixArray1D<real> dx3=data->dx[KDIR];

  IdefixArray3D<real> Ax1=data->A[IDIR];
  IdefixArray3D<real> Ax2=data->A[JDIR];
  IdefixArray3D<real> Ax3=data->A[KDIR];

  int nstart, nend;
  int nx1,nx2,nx3;

  // reconstruct BX1s
  nstart = data->beg[IDIR]-1;
  nend = data->end[IDIR];

  nx1=data->np_tot[IDIR];
  nx2=data->np_tot[JDIR];
  nx3=data->np_tot[KDIR];

  if(dir==IDIR) {
    idefix_for("ReconstructBX1s",0,nx3,0,nx2,
      KOKKOS_LAMBDA (int k, int j) {
        for(int i = nstart ; i>=0 ; i-- ) {
          Vs(BX1s,k,j,i) = 1.0 / Ax1(k,j,i) * ( Ax1(k,j,i+1)*Vs(BX1s,k,j,i+1)
                          +(D_EXPAND( ZERO_F                                                 ,
                            + Ax2(k,j+1,i) * Vs(BX2s,k,j+1,i) - Ax2(k,j,i) * Vs(BX2s,k,j,i) ,
                            + Ax3(k+1,j,i) * Vs(BX3s,k+1,j,i) - Ax3(k,j,i) * Vs(BX3s,k,j,i)  )));
        }

        for(int i = nend ; i<nx1 ; i++ ) {
          Vs(BX1s,k,j,i+1) = 1.0 / Ax1(k,j,i+1) * ( Ax1(k,j,i)*Vs(BX1s,k,j,i)
                            -(D_EXPAND(      ZERO_F                                           ,
                             + Ax2(k,j+1,i) * Vs(BX2s,k,j+1,i) - Ax2(k,j,i) * Vs(BX2s,k,j,i)  ,
                             + Ax3(k+1,j,i) * Vs(BX3s,k+1,j,i) - Ax3(k,j,i) * Vs(BX3s,k,j,i)  )));
        }
      }
    );
  }

#if DIMENSIONS >=2
  if(dir==JDIR) {
    nstart = data->beg[JDIR]-1;
    nend = data->end[JDIR];
    if(!hydro->haveAxis) {
      idefix_for("ReconstructBX2s",0,data->np_tot[KDIR],0,data->np_tot[IDIR],
        KOKKOS_LAMBDA (int k, int i) {
          for(int j = nstart ; j>=0 ; j-- ) {
            Vs(BX2s,k,j,i) = 1.0 / Ax2(k,j,i) * ( Ax2(k,j+1,i)*Vs(BX2s,k,j+1,i)
                        +(D_EXPAND( Ax1(k,j,i+1) * Vs(BX1s,k,j,i+1) - Ax1(k,j,i) * Vs(BX1s,k,j,i)  ,
                                                                                                  ,
                              + Ax3(k+1,j,i) * Vs(BX3s,k+1,j,i) - Ax3(k,j,i) * Vs(BX3s,k,j,i) )));
          }
          for(int j = nend ; j<nx2 ; j++ ) {
            Vs(BX2s,k,j+1,i) = 1.0 / Ax2(k,j+1,i) * ( Ax2(k,j,i)*Vs(BX2s,k,j,i)
                        -(D_EXPAND( Ax1(k,j,i+1) * Vs(BX1s,k,j,i+1) - Ax1(k,j,i) * Vs(BX1s,k,j,i)  ,
                                                                                                  ,
                              + Ax3(k+1,j,i) * Vs(BX3s,k+1,j,i) - Ax3(k,j,i) * Vs(BX3s,k,j,i) )));
          }
        }
      );
    } else {
      // We have an axis, ask myAxis to do that job for us
      hydro->myAxis.ReconstructBx2s();
    }
  }
#endif

#if DIMENSIONS == 3
  if(dir==KDIR) {
    nstart = data->beg[KDIR]-1;
    nend = data->end[KDIR];

    idefix_for("ReconstructBX3s",0,data->np_tot[JDIR],0,data->np_tot[IDIR],
      KOKKOS_LAMBDA (int j, int i) {
        for(int k = nstart ; k>=0 ; k-- ) {
          Vs(BX3s,k,j,i) = 1.0 / Ax3(k,j,i) * ( Ax3(k+1,j,i)*Vs(BX3s,k+1,j,i)
                          + ( Ax1(k,j,i+1) * Vs(BX1s,k,j,i+1) - Ax1(k,j,i) * Vs(BX1s,k,j,i)
                          + Ax2(k,j+1,i) * Vs(BX2s,k,j+1,i) - Ax2(k,j,i) * Vs(BX2s,k,j,i) ));
        }
        for(int k = nend ; k<nx3 ; k++ ) {
          Vs(BX3s,k+1,j,i) = 1.0 / Ax3(k+1,j,i) * ( Ax3(k,j,i)*Vs(BX3s,k,j,i)
                            - ( Ax1(k,j,i+1) * Vs(BX1s,k,j,i+1) - Ax1(k,j,i) * Vs(BX1s,k,j,i)
                            +  Ax2(k,j+1,i) * Vs(BX2s,k,j+1,i) - Ax2(k,j,i) * Vs(BX2s,k,j,i) ));
        }
      }
    );
  }
#endif

  idfx::popRegion();
}


void HydroBoundary::EnrollUserDefBoundary(UserDefBoundaryFunc myFunc) {
  this->userDefBoundaryFunc = myFunc;
  this->haveUserDefBoundary = true;
}

void HydroBoundary::EnrollInternalBoundary(InternalBoundaryFunc myFunc) {
  this->internalBoundaryFunc = myFunc;
  this->haveInternalBoundary = true;
}

void HydroBoundary::EnforcePeriodic(int dir, BoundarySide side ) {
  IdefixArray4D<real> Vc = hydro->Vc;
  int nxi = data->np_int[IDIR];
  int nxj = data->np_int[JDIR];
  int nxk = data->np_int[KDIR];

  const int ighost = data->nghost[IDIR];
  const int jghost = data->nghost[JDIR];
  const int kghost = data->nghost[KDIR];

  BoundaryForAll("BoundaryPeriodic", dir, side,
        KOKKOS_LAMBDA (int n, int k, int j, int i) {
          int iref, jref, kref;
          // This hack takes care of cases where we have more ghost zones than active zones
          if(dir==IDIR)
            iref = ighost + (i+ighost*(nxi-1))%nxi;
          else
            iref = i;
          if(dir==JDIR)
            jref = jghost + (j+jghost*(nxj-1))%nxj;
          else
            jref = j;
          if(dir==KDIR)
            kref = kghost + (k+kghost*(nxk-1))%nxk;
          else
            kref = k;

          Vc(n,k,j,i) = Vc(n,kref,jref,iref);
        });

  #if MHD==YES
    IdefixArray4D<real> Vs = hydro->Vs;
    if(dir==JDIR || dir==KDIR) {
      nxi = data->np_int[IDIR]+1;
      nxj = data->np_int[JDIR];
      nxk = data->np_int[KDIR];
      BoundaryForX1s("BoundaryPeriodicX1s",dir,side,
      KOKKOS_LAMBDA (int k, int j, int i) {
        int iref, jref, kref;
          // This hack takes care of cases where we have more ghost zones than active zones
          if(dir==IDIR)
            iref = ighost + (i+ighost*(nxi-1))%nxi;
          else
            iref = i;
          if(dir==JDIR)
            jref = jghost + (j+jghost*(nxj-1))%nxj;
          else
            jref = j;
          if(dir==KDIR)
            kref = kghost + (k+kghost*(nxk-1))%nxk;
          else
            kref = k;

          Vs(BX1s,k,j,i) = Vs(BX1s,kref,jref,iref);
      });
    }
    #if DIMENSIONS >=2
      if(dir==IDIR || dir==KDIR) {
        nxi = data->np_int[IDIR];
        nxj = data->np_int[JDIR]+1;
        nxk = data->np_int[KDIR];
        BoundaryForX2s("BoundaryPeriodicX2s",dir,side,
        KOKKOS_LAMBDA (int k, int j, int i) {
          int iref, jref, kref;
            // This hack takes care of cases where we have more ghost zones than active zones
            if(dir==IDIR)
              iref = ighost + (i+ighost*(nxi-1))%nxi;
            else
              iref = i;
            if(dir==JDIR)
              jref = jghost + (j+jghost*(nxj-1))%nxj;
            else
              jref = j;
            if(dir==KDIR)
              kref = kghost + (k+kghost*(nxk-1))%nxk;
            else
              kref = k;

            Vs(BX2s,k,j,i) = Vs(BX2s,kref,jref,iref);
        });
      }
    #endif
    #if DIMENSIONS == 3
      nxi = data->np_int[IDIR];
      nxj = data->np_int[JDIR];
      nxk = data->np_int[KDIR]+1;
      if(dir==IDIR || dir==JDIR) {
        BoundaryForX3s("BoundaryPeriodicX3s",dir,side,
        KOKKOS_LAMBDA (int k, int j, int i) {
          int iref, jref, kref;
            // This hack takes care of cases where we have more ghost zones than active zones
            if(dir==IDIR)
              iref = ighost + (i+ighost*(nxi-1))%nxi;
            else
              iref = i;
            if(dir==JDIR)
              jref = jghost + (j+jghost*(nxj-1))%nxj;
            else
              jref = j;
            if(dir==KDIR)
              kref = kghost + (k+kghost*(nxk-1))%nxk;
            else
              kref = k;

            Vs(BX3s,k,j,i) = Vs(BX3s,kref,jref,iref);
        });
      }
    #endif
  #endif// MHD
}


void HydroBoundary::EnforceReflective(int dir, BoundarySide side ) {
  IdefixArray4D<real> Vc = hydro->Vc;
  const int nxi = data->np_int[IDIR];
  const int nxj = data->np_int[JDIR];
  const int nxk = data->np_int[KDIR];

  const int ighost = data->nghost[IDIR];
  const int jghost = data->nghost[JDIR];
  const int kghost = data->nghost[KDIR];

  BoundaryForAll("BoundaryReflective", dir, side,
        KOKKOS_LAMBDA (int n, int k, int j, int i) {
          // ref= 2*ibound -i -1
          // with ibound = nghost on the left and ibount = nghost + nx -1 on the right
          const int iref = (dir==IDIR) ? 2*(ighost + side*(nxi-1)) - i - 1 : i;
          const int jref = (dir==JDIR) ? 2*(jghost + side*(nxj-1)) - j - 1 : j;
          const int kref = (dir==KDIR) ? 2*(kghost + side*(nxk-1)) - k - 1 : k;

          const int sign = (n == VX1+dir) ? -1.0 : 1.0;

          Vc(n,k,j,i) = sign * Vc(n,kref,jref,iref);
        });

  #if MHD==YES
    IdefixArray4D<real> Vs = hydro->Vs;
    if(dir==JDIR || dir==KDIR) {
      BoundaryForX1s("BoundaryReflectiveX1s",dir,side,
        KOKKOS_LAMBDA (int k, int j, int i) {
          // ref= 2*ibound -i -1
          // with ibound = nghost on the left and ibount = nghost + nx -1 on the right
          //const int iref = (dir==IDIR) ? 2*(ighost + side*(nxi-1)) - i - 1 : i;
          const int jref = (dir==JDIR) ? 2*(jghost + side*(nxj-1)) - j - 1 : j;
          const int kref = (dir==KDIR) ? 2*(kghost + side*(nxk-1)) - k - 1 : k;

          Vs(BX1s,k,j,i) = -Vs(BX1s,kref,jref,i);
        });
    }
    #if DIMENSIONS >=2
      if(dir==IDIR || dir==KDIR) {
        BoundaryForX2s("BoundaryReflectiveX2s",dir,side,
          KOKKOS_LAMBDA (int k, int j, int i) {
            const int iref = (dir==IDIR) ? 2*(ighost + side*(nxi-1)) - i - 1 : i;
            //const int jref = (dir==JDIR) ? 2*(jghost + side*(nxj-1)) - j - 1 : j;
            const int kref = (dir==KDIR) ? 2*(kghost + side*(nxk-1)) - k - 1 : k;

              Vs(BX2s,k,j,i) = -Vs(BX2s,kref,j,iref);
          });
      }
    #endif
    #if DIMENSIONS == 3
      if(dir==IDIR || dir==JDIR) {
        BoundaryForX3s("BoundaryReflectiveX3s",dir,side,
          KOKKOS_LAMBDA (int k, int j, int i) {
            const int iref = (dir==IDIR) ? 2*(ighost + side*(nxi-1)) - i - 1 : i;
            const int jref = (dir==JDIR) ? 2*(jghost + side*(nxj-1)) - j - 1 : j;
            //const int kref = (dir==KDIR) ? 2*(kghost + side*(nxk-1)) - k - 1 : k;

            Vs(BX3s,k,j,i) = -Vs(BX3s,k,jref,iref);
          });
      }
    #endif
  #endif// MHD
}

void HydroBoundary::EnforceOutflow(int dir, BoundarySide side ) {
  IdefixArray4D<real> Vc = hydro->Vc;
  const int nxi = data->np_int[IDIR];
  const int nxj = data->np_int[JDIR];
  const int nxk = data->np_int[KDIR];

  const int ighost = data->nghost[IDIR];
  const int jghost = data->nghost[JDIR];
  const int kghost = data->nghost[KDIR];

  BoundaryForAll("BoundaryOutflow", dir, side,
        KOKKOS_LAMBDA (int n, int k, int j, int i) {
          // ref= ibound
          // with ibound = nghost on the left and ibound = nghost + nx -1 on the right
          const int iref = (dir==IDIR) ? ighost + side*(nxi-1) : i;
          const int jref = (dir==JDIR) ? jghost + side*(nxj-1) : j;
          const int kref = (dir==KDIR) ? kghost + side*(nxk-1) : k;

          // should it go inwards or outwards?
          // side = 1 on the left and =-1 on the right
          const int sign = 1-2*side;

          if( (n== VX1+dir) && (sign*Vc(n,kref,jref,iref) >= ZERO_F) ) {
            Vc(n,k,j,i) = ZERO_F;
          } else {
            Vc(n,k,j,i) = Vc(n,kref,jref,iref);
          }
        });

  #if MHD==YES
    IdefixArray4D<real> Vs = hydro->Vs;
    if(dir==JDIR || dir==KDIR) {
      BoundaryForX1s("BoundaryOutflowX1s",dir,side,
        KOKKOS_LAMBDA (int k, int j, int i) {
          // with ibound = nghost on the left and ibount = nghost + nx -1 on the right
          //const int iref = (dir==IDIR) ? ighost + side*(nxi-1) : i;
          const int jref = (dir==JDIR) ? jghost + side*(nxj-1) : j;
          const int kref = (dir==KDIR) ? kghost + side*(nxk-1) : k;

          Vs(BX1s,k,j,i) = Vs(BX1s,kref,jref,i);
        });
    }
    #if DIMENSIONS >=2
      if(dir==IDIR || dir==KDIR) {
        BoundaryForX2s("BoundaryOutflowX2s",dir,side,
          KOKKOS_LAMBDA (int k, int j, int i) {
            const int iref = (dir==IDIR) ? ighost + side*(nxi-1) : i;
            //const int jref = (dir==JDIR) ? jghost + side*(nxj-1) : j;
            const int kref = (dir==KDIR) ? kghost + side*(nxk-1) : k;

            Vs(BX2s,k,j,i) = Vs(BX2s,kref,j,iref);
          });
      }
    #endif
    #if DIMENSIONS == 3
      if(dir==IDIR || dir==JDIR) {
        BoundaryForX3s("BoundaryOutflowX3s",dir,side,
          KOKKOS_LAMBDA (int k, int j, int i) {
            const int iref = (dir==IDIR) ? ighost + side*(nxi-1) : i;
            const int jref = (dir==JDIR) ? jghost + side*(nxj-1) : j;
            //const int kref = (dir==KDIR) ? kghost + side*(nxk-1) : k;

            Vs(BX3s,k,j,i) = Vs(BX3s,k,jref,iref);
          });
      }
    #endif
  #endif// MHD
}

void HydroBoundary::EnforceShearingBox(real t, int dir, BoundarySide side) {
  if(dir != IDIR)
    IDEFIX_ERROR("Shearing box boundaries can only be applied along the X1 direction");
  if(data->mygrid->nproc[JDIR]>1)
    IDEFIX_ERROR("Shearing box is not yet compatible with domain decomposition in X2");

  // First thing is to enforce periodicity (already performed by MPI)
  if(data->mygrid->nproc[dir] == 1) EnforcePeriodic(dir, side);

  IdefixArray4D<real> scrh = sBArray;
  IdefixArray4D<real> Vc = hydro->Vc;

  const int nxi = data->np_int[IDIR];
  const int nxj = data->np_int[JDIR];
  const int nxk = data->np_int[KDIR];

  const int ighost = data->nghost[IDIR];
  const int jghost = data->nghost[JDIR];
  const int kghost = data->nghost[KDIR];

  // Where does the boundary starts along x1?
  const int istart = side*(ighost+nxi);

  // Shear rate
  const real S  = hydro->sbS;

  // Box size
  const real Lx = data->mygrid->xend[IDIR] - data->mygrid->xbeg[IDIR];
  const real Ly = data->mygrid->xend[JDIR] - data->mygrid->xbeg[JDIR];

  // total number of cells in y (active domain)
  const int ny = data->mygrid->np_int[JDIR];
  const real dy = Ly/ny;

  // Compute offset in y modulo the box size
  const int sign=2*side-1;
  const real sbVelocity = sign*S*Lx;
  real dL = std::fmod(sbVelocity*t,Ly);

  // translate this into # of cells
  const int m = static_cast<int> (std::floor(dL/dy+HALF_F));

  // remainding shift
  const real eps = dL / dy - m;


  // New we need to perform the shift
  BoundaryForAll("BoundaryShearingBox", dir, side,
        KOKKOS_LAMBDA ( int n, int k, int j, int i) {
          // jorigin
          const int jo = jghost + ((j-m-jghost)%nxj+nxj)%nxj;
          const int jop2 = jghost + ((jo+2-jghost)%nxj+nxj)%nxj;
          const int jop1 = jghost + ((jo+1-jghost)%nxj+nxj)%nxj;
          const int jom1 = jghost + ((jo-1-jghost)%nxj+nxj)%nxj;
          const int jom2 = jghost + ((jo-2-jghost)%nxj+nxj)%nxj;

          // Define Left and right fluxes
          // Fluxes are defined from slop-limited interpolation
          // Using Van-leer slope limiter (consistently with the main advection scheme)
          real Fl,Fr;
          real dqm, dqp, dq;

          if(eps>=ZERO_F) {
            // Compute Fl
            dqm = Vc(n,k,jom1,i) - Vc(n,k,jom2,i);
            dqp = Vc(n,k,jo,i) - Vc(n,k,jom1,i);
            dq = (dqp*dqm > ZERO_F ? TWO_F*dqp*dqm/(dqp + dqm) : ZERO_F);

            Fl = Vc(n,k,jom1,i) + 0.5*dq*(1.0-eps);
            //Compute Fr
            dqm=dqp;
            dqp = Vc(n,k,jop1,i) - Vc(n,k,jo,i);
            dq = (dqp*dqm > ZERO_F ? TWO_F*dqp*dqm/(dqp + dqm) : ZERO_F);

            Fr = Vc(n,k,jo,i) + 0.5*dq*(1.0-eps);
          } else {
            //Compute Fl
            dqm = Vc(n,k,jo,i) - Vc(n,k,jom1,i);
            dqp = Vc(n,k,jop1,i) - Vc(n,k,jo,i);
            dq = (dqp*dqm > ZERO_F ? TWO_F*dqp*dqm/(dqp + dqm) : ZERO_F);

            Fl = Vc(n,k,jo,i) - 0.5*dq*(1.0+eps);
            // Compute Fr
            dqm=dqp;
            dqp = Vc(n,k,jop2,i) - Vc(n,k,jop1,i);
            dq = (dqp*dqm > ZERO_F ? TWO_F*dqp*dqm/(dqp + dqm) : ZERO_F);

            Fr = Vc(n,k,jop1,i) - 0.5*dq*(1.0+eps);
          }
          scrh(n,k,j,i-istart) = Vc(n,k,jo,i) - eps*(Fr - Fl);
        });
  // Copy scrach back to our boundary
  BoundaryForAll("BoundaryShearingBoxCopy", dir, side,
        KOKKOS_LAMBDA ( int n, int k, int j, int i) {
          Vc(n,k,j,i) = scrh(n,k,j,i-istart);
          if(n==VX2) Vc(n,k,j,i) += sbVelocity;
        });

  // Magnetised version of the same thing
  #if MHD==YES
    IdefixArray4D<real> Vs = hydro->Vs;
    #if DIMENSIONS >= 2
      for(int component = BX2s ; component < DIMENSIONS ; component++) {
        BoundaryFor("BoundaryShearingBoxBXs", dir, side,
        KOKKOS_LAMBDA (int k, int j, int i) {
          // jorigin
          const int jo = jghost + ((j-m-jghost)%nxj+nxj)%nxj;
          const int jop2 = jghost + ((jo+2-jghost)%nxj+nxj)%nxj;
          const int jop1 = jghost + ((jo+1-jghost)%nxj+nxj)%nxj;
          const int jom1 = jghost + ((jo-1-jghost)%nxj+nxj)%nxj;
          const int jom2 = jghost + ((jo-2-jghost)%nxj+nxj)%nxj;

          // Define Left and right fluxes
          // Fluxes are defined from slop-limited interpolation
          // Using Van-leer slope limiter (consistently with the main advection scheme)
          real Fl,Fr;
          real dqm, dqp, dq;

          if(eps>=ZERO_F) {
            // Compute Fl
            dqm = Vs(component,k,jom1,i) - Vs(component,k,jom2,i);
            dqp = Vs(component,k,jo,i) - Vs(component,k,jom1,i);
            dq = (dqp*dqm > ZERO_F ? TWO_F*dqp*dqm/(dqp + dqm) : ZERO_F);

            Fl = Vs(component,k,jom1,i) + 0.5*dq*(1.0-eps);
            //Compute Fr
            dqm=dqp;
            dqp = Vs(component,k,jop1,i) - Vs(component,k,jo,i);
            dq = (dqp*dqm > ZERO_F ? TWO_F*dqp*dqm/(dqp + dqm) : ZERO_F);

            Fr = Vs(component,k,jo,i) + 0.5*dq*(1.0-eps);
          } else {
            //Compute Fl
            dqm = Vs(component,k,jo,i) - Vs(component,k,jom1,i);
            dqp = Vs(component,k,jop1,i) - Vs(component,k,jo,i);
            dq = (dqp*dqm > ZERO_F ? TWO_F*dqp*dqm/(dqp + dqm) : ZERO_F);

            Fl = Vs(component,k,jo,i) - 0.5*dq*(1.0+eps);
            // Compute Fr
            dqm=dqp;
            dqp = Vs(component,k,jop2,i) - Vs(component,k,jop1,i);
            dq = (dqp*dqm > ZERO_F ? TWO_F*dqp*dqm/(dqp + dqm) : ZERO_F);

            Fr = Vs(component,k,jop1,i) - 0.5*dq*(1.0+eps);
          }
          scrh(0,k,j,i-istart) = Vs(component,k,jo,i) - eps*(Fr - Fl);
        });
        // Copy scratch back to our boundary
        BoundaryFor("BoundaryShearingBoxCopyBXs", dir, side,
              KOKKOS_LAMBDA ( int k, int j, int i) {
                Vs(component,k,j,i) = scrh(0,k,j,i-istart);
              });
      }// loop on components
    #endif// DIMENSIONS
  #endif // MHD
}
