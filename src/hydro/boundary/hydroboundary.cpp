// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2021 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
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
  mpi.Init(this->data, mapVars, mapNVars, true);
#endif // MPI
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
          mpi.ExchangeX1();
          break;
        case 1:
          mpi.ExchangeX2();
          break;
        case 2:
          mpi.ExchangeX3();
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

  IdefixArray4D<real> Vc = hydro->Vc;
  IdefixArray4D<real> Vs = hydro->Vs;

  int ibeg,iend,jbeg,jend,kbeg,kend;
  int ioffset,joffset,koffset;
  int ighost,jghost,kghost;

  ighost = data->nghost[IDIR];
  jghost = data->nghost[JDIR];
  kghost = data->nghost[KDIR];

  int nxi = data->np_int[IDIR];
  int nxj = data->np_int[JDIR];
  int nxk = data->np_int[KDIR];


  real sbLx = hydro->sbLx;
  real sbS  = hydro->sbS;

  ioffset = (dir == IDIR) ? data->np_int[IDIR] : 0;
  joffset = (dir == JDIR) ? data->np_int[JDIR] : 0;
  koffset = (dir == KDIR) ? data->np_int[KDIR] : 0;


  // left boundary
  ibeg=0;
  iend= (dir == IDIR) ? ighost : data->np_tot[IDIR];
  jbeg=0;
  jend= (dir == JDIR) ? jghost : data->np_tot[JDIR];
  kbeg=0;
  kend= (dir == KDIR) ? kghost : data->np_tot[KDIR];

  switch(data->lbound[dir]) {
    case internal:
      // internal is used for MPI-enforced boundary conditions. Nothing to be done here.
      break;

    case periodic:
      if(data->mygrid->nproc[dir] > 1) break; // Periodicity already enforced by MPI calls
      EnforcePeriodic(dir,left);
      break;

    case reflective:
      idefix_for("BoundaryBegReflective",0,NVAR,kbeg,kend,jbeg,jend,ibeg,iend,
        KOKKOS_LAMBDA (int n, int k, int j, int i) {
          int iref= (dir==IDIR) ? ighost : i;
          int jref= (dir==JDIR) ? jghost : j;
          int kref= (dir==KDIR) ? kghost : k;

          if( n==VX1+dir)
            Vc(n,k,j,i) = ZERO_F;
          else
            Vc(n,k,j,i) = Vc(n,kref,jref,iref);
        }
      );

#if MHD == YES
      for(int component=0; component<DIMENSIONS; component++) {
        int ieb,jeb,keb;
        if(component == IDIR)
          ieb=iend+1;
        else
          ieb=iend;
        if(component == JDIR)
          jeb=jend+1;
        else
          jeb=jend;
        if(component == KDIR)
          keb=kend+1;
        else
          keb=kend;
        if(component != dir) { // skip normal component
          idefix_for("BoundaryBegOutflowVs",kbeg,keb,jbeg,jeb,ibeg,ieb,
            KOKKOS_LAMBDA (int k, int j, int i) {
              int iref= (dir==IDIR) ? ighost : i;
              int jref= (dir==JDIR) ? jghost : j;
              int kref= (dir==KDIR) ? kghost : k;

              // Don't touch the normal component !
              Vs(component,k,j,i) = Vs(component,kref,jref,iref);
            }
          );
        }
      }
#endif
      break;

    case outflow:
      idefix_for("BoundaryBegOutflow",0,NVAR,kbeg,kend,jbeg,jend,ibeg,iend,
        KOKKOS_LAMBDA (int n, int k, int j, int i) {
          int iref= (dir==IDIR) ? ighost : i;
          int jref= (dir==JDIR) ? jghost : j;
          int kref= (dir==KDIR) ? kghost : k;

          if( (n==VX1+dir) && (Vc(n,kref,jref,iref) >= ZERO_F))
            Vc(n,k,j,i) = ZERO_F;
          else
            Vc(n,k,j,i) = Vc(n,kref,jref,iref);
        }
      );

#if MHD == YES
      for(int component=0; component<DIMENSIONS; component++) {
        int ieb,jeb,keb;
        if(component == IDIR)
          ieb=iend+1;
        else
          ieb=iend;
        if(component == JDIR)
          jeb=jend+1;
        else
          jeb=jend;
        if(component == KDIR)
          keb=kend+1;
        else
          keb=kend;
        if(component != dir) { // skip normal component
          idefix_for("BoundaryBegOutflowVs",kbeg,keb,jbeg,jeb,ibeg,ieb,
            KOKKOS_LAMBDA (int k, int j, int i) {
              int iref= (dir==IDIR) ? ighost : i;
              int jref= (dir==JDIR) ? jghost : j;
              int kref= (dir==KDIR) ? kghost : k;

              // Don't touch the normal component !
              Vs(component,k,j,i) = Vs(component,kref,jref,iref);
            }
          );
        }
      }
#endif
      break;

    case shearingbox:
      if(data->mygrid->nproc[dir] > 1) {
        // if shearing box enabled, the MPI call has already enforced strict periodicicty,
        // so we just need to enforce the offset
        real voffset=-sbLx*sbS;

        idefix_for("BoundaryBegShearingBox",kbeg,kend,jbeg,jend,ibeg,iend,
          KOKKOS_LAMBDA (int k, int j, int i) {
            Vc(VX2,k,j,i) = Vc(VX2,k,j,i) + voffset;
          }
        );
      } else {
        idefix_for("BoundaryBegShearingBox",0,NVAR,kbeg,kend,jbeg,jend,ibeg,iend,
          KOKKOS_LAMBDA (int n, int k, int j, int i) {
            real voffset= (n == VX2) ? - sbLx * sbS : ZERO_F;
            Vc(n,k,j,i) = Vc(n,k+koffset,j+joffset,i+ioffset) + voffset;
          }
        );

#if MHD == YES
        for(int component=0; component<DIMENSIONS; component++) {
          int ieb,jeb,keb;
          if(component == IDIR)
            ieb=iend+1;
          else
            ieb=iend;
          if(component == JDIR)
            jeb=jend+1;
          else
            jeb=jend;
          if(component == KDIR)
            keb=kend+1;
          else
            keb=kend;
          if(component != dir) { // skip normal component
            idefix_for("BoundaryBegShearingBoxVs",kbeg,keb,jbeg,jeb,ibeg,ieb,
              KOKKOS_LAMBDA (int k, int j, int i) {
                Vs(component,k,j,i) = Vs(component,k+koffset,j+joffset,i+ioffset);
              }
            );
          }
        }
#endif
      }
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
  ibeg= (dir == IDIR) ? ioffset + ighost : 0;
  iend = data->np_tot[IDIR];
  jbeg= (dir == JDIR) ? joffset + jghost : 0;
  jend = data->np_tot[JDIR];
  kbeg= (dir == KDIR) ? koffset + kghost : 0;
  kend = data->np_tot[KDIR];

  switch(data->rbound[dir]) {
    case internal:
      // internal is used for MPI-enforced boundary conditions. Nothing to be done here.
      break;

    case periodic:
      if(data->mygrid->nproc[dir] > 1) break; // Periodicity already enforced by MPI calls
      EnforcePeriodic(dir,right);
      break;
    case reflective:
      idefix_for("BoundaryEndReflective",0,NVAR,kbeg,kend,jbeg,jend,ibeg,iend,
        KOKKOS_LAMBDA (int n, int k, int j, int i) {
          int iref= (dir==IDIR) ? ighost + ioffset - 1 : i;
          int jref= (dir==JDIR) ? jghost + joffset - 1 : j;
          int kref= (dir==KDIR) ? kghost + koffset - 1 : k;

          if( n==VX1+dir)
            Vc(n,k,j,i) = ZERO_F;
          else
            Vc(n,k,j,i) = Vc(n,kref,jref,iref);
        }
      );

#if MHD == YES
      for(int component=0; component<DIMENSIONS; component++) {
        int ieb,jeb,keb;
        if(component == IDIR)
          ieb=iend+1;
        else
          ieb=iend;
        if(component == JDIR)
          jeb=jend+1;
        else
          jeb=jend;
        if(component == KDIR)
          keb=kend+1;
        else
          keb=kend;
        if(component != dir) { // skip normal component
          idefix_for("BoundaryEndReflectiveVs",kbeg,keb,jbeg,jeb,ibeg,ieb,
            KOKKOS_LAMBDA (int k, int j, int i) {
              int iref= (dir==IDIR) ? ighost + ioffset - 1 : i;
              int jref= (dir==JDIR) ? jghost + joffset - 1 : j;
              int kref= (dir==KDIR) ? kghost + koffset - 1 : k;
              Vs(component,k,j,i) = Vs(component,kref,jref,iref);
            }
          );
        }
      }
#endif
      break;
    case outflow:
      idefix_for("BoundaryEndOutflow",0,NVAR,kbeg,kend,jbeg,jend,ibeg,iend,
        KOKKOS_LAMBDA (int n, int k, int j, int i) {
          int iref= (dir==IDIR) ? ighost + ioffset - 1 : i;
          int jref= (dir==JDIR) ? jghost + joffset - 1 : j;
          int kref= (dir==KDIR) ? kghost + koffset - 1 : k;

          if( (n==VX1+dir) && (Vc(n,kref,jref,iref) <= ZERO_F))
            Vc(n,k,j,i) = ZERO_F;
          else
            Vc(n,k,j,i) = Vc(n,kref,jref,iref);
        }
      );

#if MHD == YES
      for(int component=0; component<DIMENSIONS; component++) {
        int ieb,jeb,keb;
        if(component == IDIR)
          ieb=iend+1;
        else
          ieb=iend;
        if(component == JDIR)
          jeb=jend+1;
        else
          jeb=jend;
        if(component == KDIR)
          keb=kend+1;
        else
          keb=kend;
        if(component != dir) { // skip normal component
          idefix_for("BoundaryEndOutflowVs",kbeg,keb,jbeg,jeb,ibeg,ieb,
            KOKKOS_LAMBDA (int k, int j, int i) {
              int iref= (dir==IDIR) ? ighost + ioffset - 1 : i;
              int jref= (dir==JDIR) ? jghost + joffset - 1 : j;
              int kref= (dir==KDIR) ? kghost + koffset - 1 : k;
              Vs(component,k,j,i) = Vs(component,kref,jref,iref);
            }
          );
        }
      }
#endif
      break;

    case shearingbox:
      if(data->mygrid->nproc[dir] > 1) {
        // if shearing box enabled, the MPI call has already enforced strict periodicicty,
        // so we just need to enforce the offset
        real voffset=sbLx*sbS;

        idefix_for("BoundaryEndShearingBox",kbeg,kend,jbeg,jend,ibeg,iend,
          KOKKOS_LAMBDA (int k, int j, int i) {
            Vc(VX2,k,j,i) = Vc(VX2,k,j,i) + voffset;
          }
        );
      } else {
        idefix_for("BoundaryEndShearingBox",0,NVAR,kbeg,kend,jbeg,jend,ibeg,iend,
          KOKKOS_LAMBDA (int n, int k, int j, int i) {
            real voffset= (n == VX2) ? + sbLx * sbS : ZERO_F;
            Vc(n,k,j,i) = Vc(n,k-koffset,j-joffset,i-ioffset) + voffset;
          }
        );

#if MHD == YES
        for(int component=0; component<DIMENSIONS; component++) {
          int ieb,jeb,keb;
          if(component == IDIR)
            ieb=iend+1;
          else
            ieb=iend;
          if(component == JDIR)
            jeb=jend+1;
          else
            jeb=jend;
          if(component == KDIR)
            keb=kend+1;
          else
            keb=kend;
          if(component != dir) { // skip normal component
            idefix_for("BoundaryEndShearingBoxVs",kbeg,keb,jbeg,jeb,ibeg,ieb,
              KOKKOS_LAMBDA (int k, int j, int i) {
                Vs(component,k,j,i) = Vs(component,k-koffset,j-joffset,i-ioffset);
              }
            );
          }
        }
#endif
      }
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
  idfx::cout << "Hydro: User-defined boundary condition has been enrolled" << std::endl;
}

void HydroBoundary::EnrollInternalBoundary(InternalBoundaryFunc myFunc) {
  this->internalBoundaryFunc = myFunc;
  this->haveInternalBoundary = true;
  idfx::cout << "Hydro: User-defined internal boundary condition has been enrolled" << std::endl;
}

void HydroBoundary::EnforcePeriodic(int dir, BoundarySide side ) {
  IdefixArray4D<real> Vc = hydro->Vc;
  const int nxi = data->np_int[IDIR];
  const int nxj = data->np_int[JDIR];
  const int nxk = data->np_int[KDIR];

  const int ighost = data->nghost[IDIR];
  const int jghost = data->nghost[JDIR];
  const int kghost = data->nghost[KDIR];

  boundary_for_all("BoundaryPeriodic", dir, side,
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
      const int nxi = data->np_int[IDIR]+1;
      boundary_for_X1s("BoundaryPeriodicX1s",dir,side,
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
    #if COMPONENTS >=2
      if(dir==IDIR || dir==KDIR) {
        const int nxj = data->np_int[JDIR]+1;
        boundary_for_X2s("BoundaryPeriodicX2s",dir,side,
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
    #if COMPONENTS == 3
      const int nxk = data->np_int[JDIR]+1;
      if(dir==IDIR || dir==JDIR) {
        boundary_for_X3s("BoundaryPeriodicX3s",dir,side,
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
  #endif// MHD
}
