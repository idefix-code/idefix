// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2021 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include "hydro.hpp"
#include "dataBlock.hpp"

void Hydro::SetBoundary(real t) {
  // set internal boundary conditions
  if(haveInternalBoundary) internalBoundaryFunc(*data, t);
  for(int dir=0 ; dir < DIMENSIONS ; dir++ ) {
      // MPI Exchange data when needed
    #ifdef WITH_MPI
    if(data->mygrid->nproc[dir]>1) {
      switch(dir) {
        case 0:
          data->mpi.ExchangeX1();
          break;
        case 1:
          data->mpi.ExchangeX2();
          break;
        case 2:
          data->mpi.ExchangeX3();
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
  ReconstructVcField(this->Vc);
#endif
}


// Enforce boundary conditions by writing into ghost zones
void Hydro::EnforceBoundaryDir(real t, int dir) {
  idfx::pushRegion("Hydro::EnforceBoundaryDir");

  IdefixArray4D<real> Vc = this->Vc;
  IdefixArray4D<real> Vs = this->Vs;

  int ibeg,iend,jbeg,jend,kbeg,kend;
  int ioffset,joffset,koffset;
  int ighost,jghost,kghost;

  ighost = data->nghost[IDIR];
  jghost = data->nghost[JDIR];
  kghost = data->nghost[KDIR];

  int nxi = data->np_int[IDIR];
  int nxj = data->np_int[JDIR];
  int nxk = data->np_int[KDIR];


  real sbLx = this->sbLx;
  real sbS  = this->sbS;

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
      idefix_for("BoundaryBegPeriodic",0,NVAR,kbeg,kend,jbeg,jend,ibeg,iend,
        KOKKOS_LAMBDA (int n, int k, int j, int i) {
          int iref, jref, kref;
          // This hack takes care of cases where we have more ghost zones than active zones
          if(dir==IDIR)
            iref = ighost + nxi - 1 - (ighost-i-1)%nxi;
          else
            iref = i;
          if(dir==JDIR)
            jref = jghost + nxj - 1 - (jghost-j-1)%nxj;
          else
            jref = j;
          if(dir==KDIR)
            kref = kghost + nxk - 1 - (kghost-k-1)%nxk;
          else
            kref = k;

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
        if(component != dir) {
          // skip normal component
          idefix_for("BoundaryBegShearingBoxVs",kbeg,keb,jbeg,jeb,ibeg,ieb,
            KOKKOS_LAMBDA (int k, int j, int i) {
              int iref, jref, kref;
              if(dir==IDIR)
                iref = ighost + nxi - 1 - (ighost-i-1)%nxi;
              else
                iref = i;
              if(dir==JDIR)
                jref = jghost + nxj - 1 - (jghost-j-1)%nxj;
              else
                jref = j;
              if(dir==KDIR)
                kref = kghost + nxk - 1 - (kghost-k-1)%nxk;
              else
                kref = k;
              Vs(component,k,j,i) = Vs(component,kref,jref,iref);
            }
          );
        }
      }
#endif
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
      this->myAxis.EnforceAxisBoundary(left);
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
      idefix_for("BoundaryEndPeriodic",0,NVAR,kbeg,kend,jbeg,jend,ibeg,iend,
        KOKKOS_LAMBDA (int n, int k, int j, int i) {
          int iref, jref, kref;
          // This hack takes care of cases where we have more ghost zones than active zones
          if(dir==IDIR)
            iref = ighost + (i-(ighost+nxi))%nxi;
          else
            iref = i;
          if(dir==JDIR)
            jref = jghost + (j-(jghost+nxj))%nxj;
          else
            jref = j;
          if(dir==KDIR)
            kref = kghost + (k-(kghost+nxk))%nxk;
          else
            kref = k;
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
          idefix_for("BoundaryEndPeriodicVs",kbeg,keb,jbeg,jeb,ibeg,ieb,
            KOKKOS_LAMBDA (int k, int j, int i) {
              int iref, jref, kref;
              if(dir==IDIR)
                iref = ighost + (i-(ighost+nxi))%nxi;
              else
                iref = i;
              if(dir==JDIR)
                jref = jghost + (j-(jghost+nxj))%nxj;
              else
                jref = j;
              if(dir==KDIR)
                kref = kghost + (k-(kghost+nxk))%nxk;
              else
                kref = k;
              Vs(component,k,j,i) = Vs(component,kref,jref,iref);
            }
          );
        }
      }
#endif
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
      this->myAxis.EnforceAxisBoundary(right);
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


void Hydro::ReconstructVcField(IdefixArray4D<real> &Vc) {
  idfx::pushRegion("Hydro::ReconstructVcField");

  IdefixArray4D<real> Vs=this->Vs;

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


void Hydro::ReconstructNormalField(int dir) {
  idfx::pushRegion("Hydro::ReconstructNormalField");

  // Reconstruct the field
  IdefixArray4D<real> Vc = this->Vc;
  IdefixArray4D<real> Vs = this->Vs;
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
    if(!haveAxis) {
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
      myAxis.ReconstructBx2s();
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
