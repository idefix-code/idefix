// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2021 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include "axis.hpp"
#include "hydro.hpp"
#include "dataBlock.hpp"


void Axis::Init(Grid &grid, Hydro *h) {
  this->hydro = h;
  this->data = this->hydro->data;
  this->emf = & this->hydro->emf;

  // Do we have a full circle?
  // Check that we have a fraction of 2PI:
  double should_be_integer = 2.0*M_PI/fabs(grid.xend[KDIR] - grid.xbeg[KDIR]);

  if(fabs(should_be_integer - round(should_be_integer))>1e-10) {
    IDEFIX_ERROR("The grid extent in X3 should be an integer fraction of 2Pi");
  }

  #if GEOMETRY != SPHERICAL
    IDEFIX_ERROR("Axis boundary conditions are only designed to handle spherical geometry");
  #endif
  idfx::cout << "Axis: Axis regularisation enabled ";

  if(fabs((grid.xend[KDIR] - grid.xbeg[KDIR] -2.0*M_PI)) < 1e-10) {
    this->isTwoPi = true;
    idfx::cout << "with full (2pi) azimuthal extension" << std::endl;;
    #ifdef WITH_MPI
      // Check that there is no domain decomposition in phi
      if(data->mygrid->nproc[KDIR]>1) {
        IDEFIX_ERROR("Axis boundaries are not compatible with MPI domain decomposition in X3");
      }
    #endif
  } else {
    this->isTwoPi = false;
    idfx::cout << "with partial (<2pi) azimuthal extension" << std::endl;
  }

  // Check where the axis is lying.
  if(hydro->data->lbound[JDIR] == axis) axisLeft = true;
  if(hydro->data->rbound[JDIR] == axis) axisRight = true;

  // Init the symmetry array (used to flip the signs of arrays accross the axis)
  symmetryVc = IdefixArray1D<int>("Axis:SymmetryVc",NVAR);
  IdefixArray1D<int>::HostMirror symmetryVcHost = Kokkos::create_mirror_view(symmetryVc);
  // Init the array
  for (int nv = 0; nv < NVAR; nv++) {
    symmetryVcHost(nv) = 1;
    if (nv == VX2)
      symmetryVcHost(nv) = -1;
    if (nv == VX3)
      symmetryVcHost(nv) = -1;
    if (nv == BX2)
      symmetryVcHost(nv) = -1;
    if (nv == BX3)
      symmetryVcHost(nv) = -1;
  }
  Kokkos::deep_copy(symmetryVc, symmetryVcHost);

#if MHD == YES
  symmetryVs = IdefixArray1D<int>("Axis:SymmetryVs",DIMENSIONS);
  IdefixArray1D<int>::HostMirror symmetryVsHost = Kokkos::create_mirror_view(symmetryVs);
  // Init the array
  for(int nv = 0; nv < DIMENSIONS; nv++) {
    symmetryVsHost(nv) = 1;
    if (nv == BX2s)
      symmetryVsHost(nv) = -1;
    if (nv == BX3s)
      symmetryVsHost(nv) = -1;
  }
  Kokkos::deep_copy(symmetryVs, symmetryVsHost);
#endif

  #if MHD == YES
  this->Ex1Avg = IdefixArray1D<real>("Axis:Ex1Avg",hydro->data->np_tot[IDIR]);
  #endif
}

void Axis::SymmetrizeEx1Side(int jref) {
#if DIMENSIONS == 3
  IdefixArray3D<real> Ex1 = emf->ex;
  IdefixAtomicArray1D<real> Ex1Avg = this->Ex1Avg;

  if(isTwoPi) {
    idefix_for("Ex1_ini",0,data->np_tot[IDIR],
        KOKKOS_LAMBDA(int i) {
          Ex1Avg(i) = ZERO_F;
        });

    idefix_for("Ex1_Symmetrize",data->beg[KDIR],data->end[KDIR],0,data->np_tot[IDIR],
      KOKKOS_LAMBDA(int k,int i) {
        Ex1Avg(i) += Ex1(k,jref,i);
      });

    int ncells=data->mygrid->np_int[KDIR];

    idefix_for("Ex1_Store",0,data->np_tot[KDIR]+KOFFSET,0,data->np_tot[IDIR],
    KOKKOS_LAMBDA(int k,int i) {
      Ex1(k,jref,i) = Ex1Avg(i)/((real) ncells);
    });
  } else {
    // if we're not doing full two pi, the flow is symmetric with respect to the axis, and the axis
    // EMF is simply zero
    idefix_for("Ex1_Store",0,data->np_tot[KDIR]+KOFFSET,0,data->np_tot[IDIR],
    KOKKOS_LAMBDA(int k,int i) {
      Ex1(k,jref,i) = ZERO_F;
    });
  }
#endif
}

// Average the Emf component along the axis
void Axis::SymmetrizeEx1() {
  idfx::pushRegion("Axis::SymmetrizeEx1");

  if(this->axisLeft) {
    int jref = hydro->data->beg[JDIR];
    this->SymmetrizeEx1Side(jref);
  }
  if(this->axisRight) {
    int jref = hydro->data->end[JDIR];
    this->SymmetrizeEx1Side(jref);
  }

  idfx::popRegion();
}

// enforce the boundary conditions on the ghost zone accross the axis
void Axis::EnforceAxisBoundary(int side) {
  idfx::pushRegion("Axis::EnforceAxisBoundary");
  IdefixArray4D<real> Vc = hydro->Vc;
  IdefixArray1D<int> sVc = this->symmetryVc;

  int ibeg = 0;
  int iend = data->np_tot[IDIR];
  int jref, jbeg,jend;
  int offset;
  if(side == left) {
    jref = data->beg[JDIR];
    jbeg = 0;
    jend = data->beg[JDIR];
    offset = -1;
  }
  if(side == right) {
    jref = data->end[JDIR]-1;
    jbeg = data->end[JDIR];
    jend = data->np_tot[JDIR];
    offset = 1;
  }

  int kbeg = 0;
  int kend = data->np_tot[KDIR];
  int np_int_k = data->mygrid->np_int[KDIR];
  int nghost_k = data->mygrid->nghost[KDIR];

  if(isTwoPi) {
    idefix_for("BoundaryEndOutflow",0,NVAR,kbeg,kend,jbeg,jend,ibeg,iend,
            KOKKOS_LAMBDA (int n, int k, int j, int i) {
              int kcomp = nghost_k + (( k - nghost_k + np_int_k/2) % np_int_k);

              Vc(n,k,j,i) = sVc(n)*Vc(n, kcomp, 2*jref-j+offset,i);
            });
  } else {
    idefix_for("BoundaryEndOutflow",0,NVAR,kbeg,kend,jbeg,jend,ibeg,iend,
            KOKKOS_LAMBDA (int n, int k, int j, int i) {
              // kcomp = k by construction since we're doing a fraction of twopi

              Vc(n,k,j,i) = sVc(n)*Vc(n, k, 2*jref-j+offset,i);
            });
  }

  #if MHD == YES
    IdefixArray4D<real> Vs = hydro->Vs;
    IdefixArray1D<int> sVs = this->symmetryVs;

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
      if(component != JDIR) { // skip normal component
        if(isTwoPi) {
          idefix_for("BoundaryEndOutflowVs",kbeg,keb,jbeg,jeb,ibeg,ieb,
            KOKKOS_LAMBDA (int k, int j, int i) {
              int kcomp = nghost_k + (( k - nghost_k + np_int_k/2) % np_int_k);
              Vs(component,k,j,i) = sVs(component)*Vs(component,kcomp, 2*jref-j+offset,i);
            }
          );
        } else {
          idefix_for("BoundaryEndOutflowVs",kbeg,keb,jbeg,jeb,ibeg,ieb,
            KOKKOS_LAMBDA (int k, int j, int i) {
              Vs(component,k,j,i) = sVs(component)*Vs(component,k, 2*jref-j+offset,i);
            }
          );
        }
      }
    }
  #endif


  idfx::popRegion();
}

// Reconstruct Bx2s taking care of the sides where an axis is lying
void Axis::ReconstructBx2s() {
  idfx::pushRegion("Axis::EnforceAxisBoundary");
  IdefixArray4D<real> Vs = hydro->Vs;
  IdefixArray3D<real> Ax1=data->A[IDIR];
  IdefixArray3D<real> Ax2=data->A[JDIR];
  IdefixArray3D<real> Ax3=data->A[KDIR];
  int nstart = data->beg[JDIR]-1;
  int nend = data->end[JDIR];
  int ntot = data->np_tot[JDIR];

  int signLeft = 1;
  int signRight = 1;
  if(axisLeft) signLeft = -1;
  if(axisRight) signRight = -1;

#if DIMENSIONS >= 2
  // This loop is a copy of ReconstructNormalField, with the proper sign when we cross the axis
  idefix_for("Axis::ReconstructBX2s",0,data->np_tot[KDIR],0,data->np_tot[IDIR],
        KOKKOS_LAMBDA (int k, int i) {
          for(int j = nstart ; j>=0 ; j-- ) {
            Vs(BX2s,k,j,i) = 1.0 / Ax2(k,j,i) * ( Ax2(k,j+1,i)*Vs(BX2s,k,j+1,i)
                        +(D_EXPAND( Ax1(k,j,i+1) * Vs(BX1s,k,j,i+1) - Ax1(k,j,i) * Vs(BX1s,k,j,i)  ,
                                                                                                  ,
                              + signLeft*( Ax3(k+1,j,i) * Vs(BX3s,k+1,j,i)
                                                              - Ax3(k,j,i) * Vs(BX3s,k,j,i) ))));
          }
          for(int j = nend ; j<ntot ; j++ ) {
            Vs(BX2s,k,j+1,i) = 1.0 / Ax2(k,j+1,i) * ( Ax2(k,j,i)*Vs(BX2s,k,j,i)
                        -(D_EXPAND( Ax1(k,j,i+1) * Vs(BX1s,k,j,i+1) - Ax1(k,j,i) * Vs(BX1s,k,j,i)  ,
                                                                                                  ,
                              + signRight*( Ax3(k+1,j,i) * Vs(BX3s,k+1,j,i)
                                                               - Ax3(k,j,i) * Vs(BX3s,k,j,i) ))));
          }
        }
      );

  // Set BX2s on the axis to the average of the two agacent cells
  // This is required since Bx2s on the axis is not evolved since
  // there is no circulation around it
    bool left = axisLeft;
    bool right = axisRight;

    int jright = data->end[JDIR];
    int jleft = data->beg[JDIR];
    if(isTwoPi) {
      idefix_for("Axis:BoundaryAvg",0,data->np_tot[KDIR],0,data->np_tot[IDIR],
            KOKKOS_LAMBDA (int k, int i) {
              if(left) {
                Vs(BX2s,k,jleft,i) = HALF_F*(Vs(BX2s,k,jleft-1,i)+Vs(BX2s,k,jleft+1,i));
              }
              if(right) {
                Vs(BX2s,k,jright,i) = HALF_F*(Vs(BX2s,k,jright-1,i)+Vs(BX2s,k,jright+1,i));
              }
            }
          );
    } else {
      idefix_for("Axis:BoundaryAvg",0,data->np_tot[KDIR],0,data->np_tot[IDIR],
            KOKKOS_LAMBDA (int k, int i) {
              if(left) {
                Vs(BX2s,k,jleft,i) = ZERO_F;
              }
              if(right) {
                Vs(BX2s,k,jright,i) = ZERO_F;
              }
            }
          );
    }
#endif


  idfx::popRegion();
}
