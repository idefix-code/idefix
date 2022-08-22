// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2022 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include <vector>
#include "axis.hpp"
#include "hydro.hpp"
#include "dataBlock.hpp"


void Axis::Init(Grid &grid, Hydro *h) {
  this->hydro = h;
  this->data = this->hydro->data;
  this->emf = & this->hydro->emf;

  #if GEOMETRY != SPHERICAL
    IDEFIX_ERROR("Axis boundary conditions are only designed to handle spherical geometry");
  #endif


  if(fabs((grid.xend[KDIR] - grid.xbeg[KDIR] -2.0*M_PI)) < 1e-10) {
    this->isTwoPi = true;
    #ifdef WITH_MPI
      // Check that there is a domain decomposition in phi
      if(data->mygrid->nproc[KDIR]>1) {
        if(data->mygrid->nproc[KDIR]%2==1) {
          IDEFIX_ERROR("The numbre of processes in the phi direction should"
                        " be even for axis decomposition");
        }
        needMPIExchange = true;
      }
    #endif
  } else {
    this->isTwoPi = false;
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
  this->BAvg = IdefixArray2D<real>("Axis:BxAvg",hydro->data->np_tot[IDIR],2);
  if(hydro->haveCurrent) {
    this->JAvg = IdefixArray2D<real>("Axis:JAvg",hydro->data->np_tot[IDIR],3);
  }
  #endif

  #ifdef WITH_MPI
    if(needMPIExchange) {
      // Make MPI exchange datatypes
      InitMPI();
    }
  #endif
}

void Axis::ShowConfig() {
  idfx::cout << "Axis: Axis regularisation ENABLED." << std::endl;
  if(isTwoPi) {
    idfx::cout << "Axis: Full 2pi regularisation around the axis." << std::endl;
    if(needMPIExchange) {
      idfx::cout << "Axis: Using MPI exchanges for axis regularisation" << std::endl;
    }
  } else {
    idfx::cout << "Axis: Fractional (2pi/N) regularisation around the axis." << std::endl;
  }
}

void Axis::SymmetrizeEx1Side(int jref) {
#if DIMENSIONS == 3
  IdefixArray3D<real> Ex1 = emf->ex;
  IdefixArray1D<real> Ex1Avg = this->Ex1Avg;

  if(isTwoPi) {
    idefix_for("Ex1_ini",0,data->np_tot[IDIR],
        KOKKOS_LAMBDA(int i) {
          Ex1Avg(i) = ZERO_F;
        });

    idefix_for("Ex1_Symmetrize",data->beg[KDIR],data->end[KDIR],0,data->np_tot[IDIR],
      KOKKOS_LAMBDA(int k,int i) {
        Kokkos::atomic_add(&Ex1Avg(i),  Ex1(k,jref,i));
      });
    if(needMPIExchange) {
      #ifdef WITH_MPI
        // sum along all of the processes on the same r
        MPI_Allreduce(MPI_IN_PLACE, Ex1Avg.data(), data->np_tot[IDIR], realMPI,
                      MPI_SUM, data->mygrid->AxisComm);
      #endif
    }

    int ncells=data->mygrid->np_int[KDIR];

    idefix_for("Ex1_Store",0,data->np_tot[KDIR],0,data->np_tot[IDIR],
    KOKKOS_LAMBDA(int k,int i) {
      Ex1(k,jref,i) = Ex1Avg(i)/((real) ncells);
    });
  } else {
    // if we're not doing full two pi, the flow is symmetric with respect to the axis, and the axis
    // EMF is simply zero
    idefix_for("Ex1_Store",0,data->np_tot[KDIR],0,data->np_tot[IDIR],
    KOKKOS_LAMBDA(int k,int i) {
      Ex1(k,jref,i) = ZERO_F;
    });
  }
#endif
}

// Ex3 (=Ephi) on the axis is ill-defined. However, the length of cell edges along the phi
// direction is zero on the axis, so this EMF is not relevent when using CT.
// Nevertheless, when using a vector potential formulation, Ex3 on the axis may pile up
// in Ve(AX3e...), leading potentially numerical instabilities in that region.
// Hence, we enforce a regularisation of Ex3 for consistancy.

void Axis::RegularizeEx3side(int jref) {
  IdefixArray3D<real> Ex3 = emf->ez;

  idefix_for("Ex3_Regularise",0,data->np_tot[KDIR],0,data->np_tot[IDIR],
    KOKKOS_LAMBDA(int k,int i) {
      Ex3(k,jref,i) = 0.0;
    });
}

void Axis::RegularizeCurrentSide(int side) {
  // Compute the values of Jx, Jy and Jz that are consistent for all cells touching the axis
  #if DIMENSIONS == 3
    IdefixArray4D<real> J = hydro->J;
    int jref = 0;
    int sign = 0;

    if(side == left) {
      jref = data->beg[JDIR];
      sign = 1;
    }
    if(side == right) {
      jref = data->end[JDIR];
      sign = -1;
    }
    IdefixArray2D<real> JAvg = this->JAvg;
    IdefixArray1D<real> phi = data->x[KDIR];



    idefix_for("J_ini",0,data->np_tot[IDIR],0,3,
          KOKKOS_LAMBDA(int i, int n) {
            JAvg(i,n) = ZERO_F;
    });
    idefix_for("JHorizontal_compute",data->beg[KDIR],data->end[KDIR],0,data->np_tot[IDIR],
        KOKKOS_LAMBDA(int k,int i) {
          real Jthmid = sign*(J(JDIR,k  ,jref-1,i) +
                              J(JDIR,k  ,jref  ,i) +
                              J(JDIR,k+1,jref-1,i) +
                              J(JDIR,k+1,jref  ,i)) / 4.0;

          real Jphimid = J(KDIR,k, jref, i);

          Kokkos::atomic_add(&JAvg(i,IDIR), Jthmid * cos(phi(k)) - Jphimid * sin(phi(k)));
          Kokkos::atomic_add(&JAvg(i,JDIR), Jthmid * sin(phi(k)) + Jphimid * cos(phi(k)));
          Kokkos::atomic_add(&JAvg(i,KDIR), J(IDIR,k,jref+sign,i)); // We pick up the radial current
                                                                    // in the active zones
    });

    if(needMPIExchange) {
      #ifdef WITH_MPI
        // sum along all of the processes on the same r
        MPI_Allreduce(MPI_IN_PLACE, JAvg.data(), 3*data->np_tot[IDIR], realMPI,
                      MPI_SUM, data->mygrid->AxisComm);
      #endif
    }

    const int ncells=data->mygrid->np_int[KDIR];

    idefix_for("fixJ",0,data->np_tot[KDIR],0,data->np_tot[IDIR],
        KOKKOS_LAMBDA(int k,int i) {
          real Jx = JAvg(i,IDIR) / ((real) ncells*sign);
          real Jy = JAvg(i,JDIR) / ((real) ncells*sign);
          real Jz = JAvg(i,KDIR) / ((real) ncells);

          J(IDIR, k,jref,i) = Jz;
          J(JDIR, k,jref,i) = cos(phi(k))*Jx + sin(phi(k))*Jy;
          // There is nothing along KDIR since Jphi is never localised on the axis.
        });

  #endif // DIMENSIONS
}

// Average the Emf component along the axis
void Axis::RegularizeEMFs() {
  idfx::pushRegion("Axis::RegularizeEMFs");

  if(this->axisLeft) {
    int jref = hydro->data->beg[JDIR];
    SymmetrizeEx1Side(jref);
    RegularizeEx3side(jref);
  }
  if(this->axisRight) {
    int jref = hydro->data->end[JDIR];
    SymmetrizeEx1Side(jref);
    RegularizeEx3side(jref);
  }

  idfx::popRegion();
}

// Average the Emf component along the axis
void Axis::RegularizeCurrent() {
  idfx::pushRegion("Axis::RegularizeCurrent");

  if(this->axisLeft) {
    RegularizeCurrentSide(left);
  }
  if(this->axisRight) {
    RegularizeCurrentSide(right);
  }

  idfx::popRegion();
}

void Axis::FixBx2sAxis(int side) {
  // Compute the values of Bx and By that are consistent with BX2 along the axis
  #if DIMENSIONS == 3
    IdefixArray4D<real> Vs = hydro->Vs;
    IdefixArray2D<real> BAvg = this->BAvg;
    IdefixArray1D<real> phi = data->x[KDIR];

    int jref = 0;
    int sign = 0;

    if(side == left) {
      jref = data->beg[JDIR];
      sign = 1;
    }
    if(side == right) {
      jref = data->end[JDIR];
      sign = -1;
    }

    idefix_for("B_ini",0,data->np_tot[IDIR],0,2,
          KOKKOS_LAMBDA(int i, int n) {
            BAvg(i,n) = ZERO_F;
    });
    idefix_for("BHorizontal_compute",data->beg[KDIR],data->end[KDIR],0,data->np_tot[IDIR],
        KOKKOS_LAMBDA(int k,int i) {
          real Bthmid = sign*HALF_F*(Vs(BX2s,k,jref-1,i) + Vs(BX2s,k,jref+1,i));
          real Bphimid = HALF_F*(Vs(BX3s,k,jref-1,i) + Vs(BX3s,k,jref,i));
          //Bthmid = 0.0;
          //Bphimid = 0.0;

          Kokkos::atomic_add(&BAvg(i,IDIR), Bthmid * cos(phi(k)) - Bphimid * sin(phi(k)));
          Kokkos::atomic_add(&BAvg(i,JDIR), Bthmid * sin(phi(k)) + Bphimid * cos(phi(k)));
    });
    if(needMPIExchange) {
      #ifdef WITH_MPI
        // sum along all of the processes on the same r
        MPI_Allreduce(MPI_IN_PLACE, BAvg.data(), 2*data->np_tot[IDIR], realMPI,
                      MPI_SUM, data->mygrid->AxisComm);
      #endif
    }
    int ncells=data->mygrid->np_int[KDIR];

    idefix_for("fixBX2s",data->beg[KDIR],data->end[KDIR],0,data->np_tot[IDIR],
        KOKKOS_LAMBDA(int k,int i) {
          real Bx = BAvg(i,IDIR) / ((real) ncells*sign);
          real By = BAvg(i,JDIR) / ((real) ncells*sign);

          Vs(BX2s,k,jref,i) = cos(phi(k))*Bx + sin(phi(k))*By;
        });
  #endif // DIMENSIONS
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
    if(needMPIExchange) {
      ExchangeMPI(side);
    } else { // no MPI exchange
      idefix_for("BoundaryAxis",0,NVAR,kbeg,kend,jbeg,jend,ibeg,iend,
              KOKKOS_LAMBDA (int n, int k, int j, int i) {
                int kcomp = nghost_k + (( k - nghost_k + np_int_k/2) % np_int_k);

                Vc(n,k,j,i) = sVc(n)*Vc(n, kcomp, 2*jref-j+offset,i);
              });
    }// MPI Exchange
  } else {  // not 2pi
    idefix_for("BoundaryAxis",0,NVAR,kbeg,kend,jbeg,jend,ibeg,iend,
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
          if(!needMPIExchange) { //no mpi exchange
            idefix_for("BoundaryAxisVs",kbeg,keb,jbeg,jeb,ibeg,ieb,
              KOKKOS_LAMBDA (int k, int j, int i) {
                int kcomp = nghost_k + (( k - nghost_k + np_int_k/2) % np_int_k);
                Vs(component,k,j,i) = sVs(component)*Vs(component,kcomp, 2*jref-j+offset,i);
              }
            );
          } // mpi exchange
        } else { // not 2pi
          idefix_for("BoundaryAxisVs",kbeg,keb,jbeg,jeb,ibeg,ieb,
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
  idfx::pushRegion("Axis::ReconstructBx2s");
#if DIMENSIONS >= 2 && MHD == YES
  IdefixArray4D<real> Vs = hydro->Vs;
  IdefixArray3D<real> Ax1=data->A[IDIR];
  IdefixArray3D<real> Ax2=data->A[JDIR];
  IdefixArray3D<real> Ax3=data->A[KDIR];
  int nstart = data->beg[JDIR]-1;
  int nend = data->end[JDIR];
  int ntot = data->np_tot[JDIR];

  [[maybe_unused]] int signLeft = 1;
  [[maybe_unused]] int signRight = 1;
  if(axisLeft) signLeft = -1;
  if(axisRight) signRight = -1;

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
    bool haveleft = axisLeft;
    bool haveright = axisRight;

    int jright = data->end[JDIR];
    int jleft = data->beg[JDIR];
    if(isTwoPi) {
      if(haveleft) FixBx2sAxis(left);
      if(haveright) FixBx2sAxis(right);

    } else {
      idefix_for("Axis:BoundaryAvg",0,data->np_tot[KDIR],0,data->np_tot[IDIR],
            KOKKOS_LAMBDA (int k, int i) {
              if(haveleft) {
                Vs(BX2s,k,jleft,i) = ZERO_F;
              }
              if(haveright) {
                Vs(BX2s,k,jright,i) = ZERO_F;
              }
            }
          );
    }
#endif
  idfx::popRegion();
}



void Axis::ExchangeMPI(int side) {
  idfx::pushRegion("Axis::ExchangeMPI");
  #ifdef WITH_MPI
  // Load  the buffers with data
  int ibeg,iend,jbeg,jend,kbeg,kend,offset;
  int nx,ny,nz;
  auto bufferSend = this->bufferSend;
  IdefixArray1D<int> map = this->mapVars;
  IdefixArray4D<real> Vc = hydro->Vc;
  IdefixArray4D<real> Vs = hydro->Vs;

// If MPI Persistent, start receiving even before the buffers are filled

  MPI_Status sendStatus;
  MPI_Status recvStatus;

  double tStart = MPI_Wtime();
  MPI_SAFE_CALL(MPI_Start(&recvRequest));
  idfx::mpiCallsTimer += MPI_Wtime() - tStart;

  // Coordinates of the ghost region which needs to be transfered
  ibeg   = 0;
  iend   = data->np_tot[IDIR];
  nx     = data->np_tot[IDIR];  // Number of points in x
  jbeg   = 0;
  jend   = data->nghost[JDIR];
  offset = data->end[JDIR];     // Distance between beginning of left and right ghosts
  ny     = data->nghost[JDIR];
  kbeg   = data->beg[KDIR];
  kend   = data->end[KDIR];
  nz     = kend - kbeg;
  if(side==left) {
    idefix_for("LoadBufferX2Vc",0,mapNVars,kbeg,kend,jbeg,jend,ibeg,iend,
      KOKKOS_LAMBDA (int n, int k, int j, int i) {
        bufferSend(i + (j-jbeg)*nx + (k-kbeg)*nx*ny + n*nx*ny*nz ) =
                                                          Vc(map(n),k,j+ny,i);
      }
    );
    #if MHD == YES
      int VsIndex = mapNVars*nx*ny*nz;
      idefix_for("LoadBufferX2IDIR",kbeg,kend,jbeg,jend,ibeg,iend+1,
        KOKKOS_LAMBDA (int k, int j, int i) {
          bufferSend(i + (j-jbeg)*(nx+1) + (k-kbeg)*(nx+1)*ny + VsIndex ) =
                                                          Vs(IDIR,k,j+ny,i);
        }
      );
      VsIndex = mapNVars*nx*ny*nz + (nx+1)*ny*nz;
      idefix_for("LoadBufferX2KDIR",kbeg,kend+1,jbeg,jend,ibeg,iend,
        KOKKOS_LAMBDA (int k, int j, int i) {
          bufferSend(i + (j-jbeg)*nx + (k-kbeg)*nx*ny + VsIndex ) =
                                                          Vs(KDIR,k,j+ny,i);
        }
      );
    #endif
  } else if(side==right) {
    idefix_for("LoadBufferX2Vc",0,mapNVars,kbeg,kend,jbeg,jend,ibeg,iend,
      KOKKOS_LAMBDA (int n, int k, int j, int i) {
        bufferSend(i + (j-jbeg)*nx + (k-kbeg)*nx*ny + n*nx*ny*nz ) =
                                                          Vc(map(n),k,j+offset-ny,i);
      }
    );

    // Load face-centered field in the buffer
     #if MHD == YES
      int VsIndex = mapNVars*nx*ny*nz;
      idefix_for("LoadBufferX2IDIR",kbeg,kend,jbeg,jend,ibeg,iend+1,
        KOKKOS_LAMBDA (int k, int j, int i) {
          bufferSend(i + (j-jbeg)*(nx+1) + (k-kbeg)*(nx+1)*ny + VsIndex ) =
                                                          Vs(IDIR,k,j+offset-ny,i);
        }
      );
      VsIndex = mapNVars*nx*ny*nz + (nx+1)*ny*nz;

      idefix_for("LoadBufferX2KDIR",kbeg,kend+1,jbeg,jend,ibeg,iend,
        KOKKOS_LAMBDA (int k, int j, int i) {
          bufferSend(i + (j-jbeg)*nx + (k-kbeg)*nx*ny + VsIndex ) =
                                                          Vs(KDIR,k,j+offset-ny,i);
        }
      );
    #endif // MHD
  } // side==right

  Kokkos::fence();

  tStart = MPI_Wtime();
  MPI_SAFE_CALL(MPI_Start(&sendRequest));
  MPI_Wait(&recvRequest,&recvStatus);
  idfx::mpiCallsTimer += MPI_Wtime() - tStart;

  // Unpack
  auto bufferRecv=this->bufferRecv;
  auto sVc = this->symmetryVc;

  if(side==left) {
    idefix_for("StoreBufferAxis",0,mapNVars,kbeg,kend,jbeg,jend,ibeg,iend,
      KOKKOS_LAMBDA (int n, int k, int j, int i) {
        Vc(map(n),k,jend-(j-jbeg)-1,i) =
                    sVc(map(n))*bufferRecv(i + (j-jbeg)*nx + (k-kbeg)*nx*ny + n*nx*ny*nz );
      }
    );

    // Load face-centered field in the buffer
     #if MHD == YES
      int VsIndex = mapNVars*nx*ny*nz;
      auto sVs = this->symmetryVs;
      idefix_for("StoreBufferX2IDIR",kbeg,kend,jbeg,jend,ibeg,iend+1,
        KOKKOS_LAMBDA (int k, int j, int i) {
          Vs(IDIR,k,jend-(j-jbeg)-1,i) =
                    sVs(IDIR)*bufferRecv(i + (j-jbeg)*(nx+1) + (k-kbeg)*(nx+1)*ny + VsIndex );
        }
      );

      VsIndex = mapNVars*nx*ny*nz + (nx+1)*ny*nz;

      idefix_for("StoreBufferX2KDIR",kbeg,kend+1,jbeg,jend,ibeg,iend,
        KOKKOS_LAMBDA ( int k, int j, int i) {
          Vs(KDIR,k,jend-(j-jbeg)-1,i) =
                      sVs(KDIR)*bufferRecv(i + (j-jbeg)*nx + (k-kbeg)*nx*ny + VsIndex );
        }
      );
    #endif //MHD
  } else if(side==right) {
    idefix_for("StoreBufferAxis",0,mapNVars,kbeg,kend,jbeg,jend,ibeg,iend,
      KOKKOS_LAMBDA (int n, int k, int j, int i) {
        Vc(map(n),k,jend-(j-jbeg)-1+offset,i) =
                    sVc(map(n))*bufferRecv(i + (j-jbeg)*nx + (k-kbeg)*nx*ny + n*nx*ny*nz );
      }
    );

    // Load face-centered field in the buffer
    #if MHD == YES
    int VsIndex = mapNVars*nx*ny*nz;
    auto sVs = this->symmetryVs;
    idefix_for("StoreBufferX2IDIR",kbeg,kend,jbeg,jend,ibeg,iend+1,
      KOKKOS_LAMBDA (int k, int j, int i) {
        Vs(IDIR,k,jend-(j-jbeg)-1+offset,i) =
                    sVs(IDIR)*bufferRecv(i + (j-jbeg)*(nx+1) + (k-kbeg)*(nx+1)*ny + VsIndex );
      }
    );
    VsIndex = mapNVars*nx*ny*nz + (nx+1)*ny*nz;

    idefix_for("StoreBufferX2KDIR",kbeg,kend+1,jbeg,jend,ibeg,iend,
      KOKKOS_LAMBDA ( int k, int j, int i) {
        Vs(KDIR,k,jend-(j-jbeg)-1+offset,i) =
                    sVs(KDIR)*bufferRecv(i + (j-jbeg)*nx + (k-kbeg)*nx*ny + VsIndex );
      }
    );
    #endif
  }

  MPI_Wait(&sendRequest, &sendStatus);

  idfx::mpiCallsTimer += MPI_Wtime() - tStart;


  #endif  //MPI
  idfx::popRegion();
}

void Axis::InitMPI() {
  idfx::pushRegion("Axis::InitMPI");
  #ifdef WITH_MPI

  ////////////////////////////////////////////////////////////////////////////
  // Init variable mappers
  // The variable mapper list all of the variable which are exchanged in MPI boundary calls
  // This is required since we skip some of the variables in Vc to limit the amount of data
  // being exchanged
  #if MHD == YES
  this->mapNVars = NVAR - DIMENSIONS; // We will not send magnetic field components which are in Vs
  #else
  this->mapNVars = NVAR;
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

  this->mapVars = idfx::ConvertVectorToIdefixArray(mapVars);

  this->bufferSize = data->np_tot[IDIR] * data->nghost[JDIR] * data->np_int[KDIR] * mapNVars;
  #if MHD == YES
    // IDIR
    bufferSize += (data->np_tot[IDIR]+1) * data->nghost[JDIR] * data->np_int[KDIR];
    #if DIMENSIONS==3
    bufferSize += data->np_tot[IDIR] * data->nghost[JDIR] * (data->np_int[KDIR]+1);
    #endif  // DIMENSIONS
  #endif

  this->bufferRecv = IdefixArray1D<real>("bufferRecvAxis", bufferSize);
  this->bufferSend = IdefixArray1D<real>("bufferSendAxis", bufferSize);

  // init persistent communications
  // We receive from procRecv, and we send to procSend
  int procSend, procRecv;

  // We receive from procRecv, and we send to procSend, send to the right. Shift by half the domain
  MPI_SAFE_CALL(MPI_Cart_shift(data->mygrid->AxisComm,0,data->mygrid->nproc[KDIR]/2,
                               &procRecv,&procSend ));

  MPI_SAFE_CALL(MPI_Send_init(bufferSend.data(), bufferSize, realMPI, procSend,
                650, data->mygrid->AxisComm, &sendRequest));

  MPI_SAFE_CALL(MPI_Recv_init(bufferRecv.data(), bufferSize, realMPI, procRecv,
                650, data->mygrid->AxisComm, &recvRequest));

  #endif
  idfx::popRegion();
}
