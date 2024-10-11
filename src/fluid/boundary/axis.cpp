// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include <vector>
#include "axis.hpp"
#include "boundary.hpp"

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

  IdefixArray3D<real> Ex1 = this->ex;
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
        Kokkos::fence();
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
  IdefixArray3D<real> Ex3 = this->ez;

  idefix_for("Ex3_Regularise",0,data->np_tot[KDIR],0,data->np_tot[IDIR],
    KOKKOS_LAMBDA(int k,int i) {
      Ex3(k,jref,i) = 0.0;
    });
}


void Axis::RegularizeCurrentSide(int side) {
  // Compute the values of Jx, Jy and Jz that are consistent for all cells touching the axis
  #if DIMENSIONS == 3
    IdefixArray4D<real> J = this->J;
    IdefixArray4D<real> Vs = this->Vs;
    int js = 0;
    int jc = 0;
    int sign = 0;
    if(side == left) {
      js = data->beg[JDIR];
      jc = data->beg[JDIR];
      sign = 1;
    }
    if(side == right) {
      js = data->end[JDIR];
      jc = data->end[JDIR]-1;
      sign = -1;
    }
    IdefixArray1D<real> BAvg = this->Ex1Avg;
    IdefixArray1D<real> x2 = data->x[JDIR];
    IdefixArray1D<real> x1 = data->x[IDIR];
    IdefixArray1D<real> dx3 = data->dx[KDIR];

    idefix_for("B_ini",0,data->np_tot[IDIR],
          KOKKOS_LAMBDA(int i) {
            BAvg(i) = ZERO_F;
    });
    idefix_for("Compute_Bcirculation",data->beg[KDIR],data->end[KDIR],0,data->np_tot[IDIR],
        KOKKOS_LAMBDA(int k,int i) {
          Kokkos::atomic_add(&BAvg(i), Vs(BX3s,k,jc,i)*dx3(k) ); // Compute the circulation of
                                                                 // Bphi around the pole
    });

    if(needMPIExchange) {
      #ifdef WITH_MPI
        Kokkos::fence();
        // sum along all of the processes on the same r
        MPI_Allreduce(MPI_IN_PLACE, BAvg.data(), data->np_tot[IDIR], realMPI,
                      MPI_SUM, data->mygrid->AxisComm);
      #endif
    }

    const int ncells=data->mygrid->np_int[KDIR];

    real deltaPhi = data->mygrid->xend[KDIR] - data->mygrid->xbeg[KDIR];

    // Use the circulation around the pole of Bphi to determine Jr on the pole:
    // Delta phi r^2(1-cos theta) Jr = int r sin(theta) Bphi dphi

    idefix_for("fixJ",0,data->np_tot[KDIR],0,data->np_tot[IDIR],
        KOKKOS_LAMBDA(int k,int i) {
          real th = x2(jc);
          real fact = sign*sin(th)/(deltaPhi*x1(i)*(1-cos(th)));
          J(IDIR, k,js,i) = BAvg(i)*fact;
        });

  #endif // DIMENSIONS
}

// Average the Emf component along the axis

void Axis::RegularizeEMFs() {
  idfx::pushRegion("Axis::RegularizeEMFs");

  if(this->axisLeft) {
    int jref = data->beg[JDIR];
    SymmetrizeEx1Side(jref);
    RegularizeEx3side(jref);
  }
  if(this->axisRight) {
    int jref = data->end[JDIR];
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
    IdefixArray4D<real> Vs = this->Vs;
    IdefixArray2D<real> BAvg = this->BAvg;
    IdefixArray1D<real> phi = data->x[KDIR];

    int jin = 0;
    int jout = 0;
    int jaxe = 0;
    int sign = 0;

    if(side == left) {
      jin = data->beg[JDIR];
      jout = jin-1;
      jaxe = jin;
      sign = 1;
    }
    if(side == right) {
      jin = data->end[JDIR]-1;
      jout = jin+1;
      jaxe = jout;
      sign = -1;
    }

    idefix_for("B_ini",0,data->np_tot[IDIR],0,2,
          KOKKOS_LAMBDA(int i, int n) {
            BAvg(i,n) = ZERO_F;
    });
    idefix_for("BHorizontal_compute",data->beg[KDIR],data->end[KDIR],0,data->np_tot[IDIR],
        KOKKOS_LAMBDA(int k,int i) {
          real Bthmid = sign*HALF_F*(Vs(BX2s,k,jaxe-1,i) + Vs(BX2s,k,jaxe+1,i));
          real Bphimid = HALF_F*(Vs(BX3s,k,jin,i) + Vs(BX3s,k,jout,i));
          //Bthmid = 0.0;
          //Bphimid = 0.0;

          Kokkos::atomic_add(&BAvg(i,IDIR), Bthmid * cos(phi(k)) - Bphimid * sin(phi(k)));
          Kokkos::atomic_add(&BAvg(i,JDIR), Bthmid * sin(phi(k)) + Bphimid * cos(phi(k)));
    });
    if(needMPIExchange) {
      Kokkos::fence();
      #ifdef WITH_MPI
        // sum along all of the processes on the same r
        MPI_Allreduce(MPI_IN_PLACE, BAvg.data(), 2*data->np_tot[IDIR], realMPI,
                      MPI_SUM, data->mygrid->AxisComm);
      #endif
    }
    int ncells=data->mygrid->np_int[KDIR];

    idefix_for("fixBX2s",data->beg[KDIR],data->end[KDIR],0,data->np_tot[IDIR],
        KOKKOS_LAMBDA(int k,int i) {
          real Bx = BAvg(i,IDIR) / ((real) ncells);
          real By = BAvg(i,JDIR) / ((real) ncells);

          Vs(BX2s,k,jaxe,i) = sign*(cos(phi(k))*Bx + sin(phi(k))*By);
        });
  #endif // DIMENSIONS
}


void Axis::FixBx2sAxisGhostAverage(int side) {
  // This uses the same method as Athena (Stone+????) to enforce the BX2 value
  // on the axis :
  // average of BX2s on the left face of the last ghost cell and right face of
  // the first active cell (for the left-side bondary, conversely for the
  // right-side boundary)

  #if DIMENSIONS == 3
    IdefixArray4D<real> Vs = this->Vs;

    int jaxis = 0;

    if(side == left) {
      jaxis = data->beg[JDIR];
    }
    if(side == right) {
      jaxis = data->end[JDIR];
    }

    idefix_for("BHorizontal_averaging",data->beg[KDIR],data->end[KDIR],0,data->np_tot[IDIR],
        KOKKOS_LAMBDA(int k,int i) {
          Vs(BX2s,k,jaxis,i) = HALF_F*(Vs(BX2s,k,jaxis-1,i) + Vs(BX2s,k,jaxis+1,i));
        });
  #endif // DIMENSIONS
}




// enforce the boundary conditions on the ghost zone accross the axis

void Axis::EnforceAxisBoundary(int side) {
  idfx::pushRegion("Axis::EnforceAxisBoundary");
  IdefixArray4D<real> Vc = this->Vc;
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
      idefix_for("BoundaryAxis",0,this->nVar,kbeg,kend,jbeg,jend,ibeg,iend,
              KOKKOS_LAMBDA (int n, int k, int j, int i) {
                int kcomp = nghost_k + (( k - nghost_k + np_int_k/2) % np_int_k);

                Vc(n,k,j,i) = sVc(n)*Vc(n, kcomp, 2*jref-j+offset,i);
              });
    }// MPI Exchange
  } else {  // not 2pi
    idefix_for("BoundaryAxis",0,this->nVar,kbeg,kend,jbeg,jend,ibeg,iend,
            KOKKOS_LAMBDA (int n, int k, int j, int i) {
              // kcomp = k by construction since we're doing a fraction of twopi

              Vc(n,k,j,i) = sVc(n)*Vc(n, k, 2*jref-j+offset,i);
            });
  }

  if(haveMHD) {
    IdefixArray4D<real> Vs = this->Vs;
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
  }


  idfx::popRegion();
}

// Reconstruct Bx2s taking care of the sides where an axis is lying

void Axis::ReconstructBx2s() {
  idfx::pushRegion("Axis::ReconstructBx2s");
#if DIMENSIONS >= 2 && MHD == YES
  IdefixArray4D<real> Vs = this->Vs;
  IdefixArray3D<real> Ax1=data->A[IDIR];
  IdefixArray3D<real> Ax2=data->A[JDIR];
  IdefixArray3D<real> Ax3=data->A[KDIR];
  int nstart = data->beg[JDIR]-1;
  int nend = data->end[JDIR];
  int ntot = data->np_tot[JDIR];

  bool haveleft = axisLeft;
  bool haveright = axisRight;

  if(haveleft) {
    // This loop is a copy of ReconstructNormalField, with the proper sign when we cross the axis
    idefix_for("Axis::ReconstructBX2sLeft",0,data->np_tot[KDIR],0,data->np_tot[IDIR],
      KOKKOS_LAMBDA (int k, int i) {
        for(int j = nstart ; j>=0 ; j-- ) {
          Vs(BX2s,k,j,i) = 1.0 / Ax2(k,j,i) * ( Ax2(k,j+1,i)*Vs(BX2s,k,j+1,i)
                      +(D_EXPAND( Ax1(k,j,i+1) * Vs(BX1s,k,j,i+1) - Ax1(k,j,i) * Vs(BX1s,k,j,i)  ,
                                                                                                ,
                            - ( Ax3(k+1,j,i) * Vs(BX3s,k+1,j,i)
                                                            - Ax3(k,j,i) * Vs(BX3s,k,j,i) ))));
        }
      }
    );
  }
  if(haveright) {
    // This loop is a copy of ReconstructNormalField, with the proper sign when we cross the axis
    idefix_for("Axis::ReconstructBX2sRight",0,data->np_tot[KDIR],0,data->np_tot[IDIR],
      KOKKOS_LAMBDA (int k, int i) {
        for(int j = nend ; j<ntot ; j++ ) {
          Vs(BX2s,k,j+1,i) = 1.0 / Ax2(k,j+1,i) * ( Ax2(k,j,i)*Vs(BX2s,k,j,i)
                      -(D_EXPAND( Ax1(k,j,i+1) * Vs(BX1s,k,j,i+1) - Ax1(k,j,i) * Vs(BX1s,k,j,i)  ,
                                                                                                ,
                            -  ( Ax3(k+1,j,i) * Vs(BX3s,k+1,j,i)
                                                            - Ax3(k,j,i) * Vs(BX3s,k,j,i) ))));
        }
      }
    );
  }

    // Set BX2s on the axis to the average of the two agacent cells
    // This is required since Bx2s on the axis is not evolved since
    // there is no circulation around it


    int jright = data->end[JDIR];
    int jleft = data->beg[JDIR];
    if(isTwoPi) {
      #ifdef AXIS_BX2S_USE_ATHENA_REGULARISATION
        if(haveleft) FixBx2sAxisGhostAverage(left);
        if(haveright) FixBx2sAxisGhostAverage(right);
      #else
        if(haveleft) FixBx2sAxis(left);
        if(haveright) FixBx2sAxis(right);
      #endif
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
  IdefixArray4D<real> Vc = this->Vc;
  IdefixArray4D<real> Vs = this->Vs;

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
    if (haveMHD) {
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
    } // MHD
  } else if(side==right) {
    idefix_for("LoadBufferX2Vc",0,mapNVars,kbeg,kend,jbeg,jend,ibeg,iend,
      KOKKOS_LAMBDA (int n, int k, int j, int i) {
        bufferSend(i + (j-jbeg)*nx + (k-kbeg)*nx*ny + n*nx*ny*nz ) =
                                                          Vc(map(n),k,j+offset-ny,i);
      }
    );

    // Load face-centered field in the buffer
     if (haveMHD) {
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
    } // MHD
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
    if (haveMHD) {
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
    }
  } else if(side==right) {
    idefix_for("StoreBufferAxis",0,mapNVars,kbeg,kend,jbeg,jend,ibeg,iend,
      KOKKOS_LAMBDA (int n, int k, int j, int i) {
        Vc(map(n),k,jend-(j-jbeg)-1+offset,i) =
                    sVc(map(n))*bufferRecv(i + (j-jbeg)*nx + (k-kbeg)*nx*ny + n*nx*ny*nz );
      }
    );

    // Load face-centered field in the buffer
    if (haveMHD) {
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
    } // MHD
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
  if (haveMHD) {
    this->mapNVars = this->nVar - DIMENSIONS; // We don't send B components (already in Vs)
  } else {
    this->mapNVars = this->nVar;
  }

  std::vector<int> mapVars;
  // Init the list of variables we will exchange in MPI routines
  int ntarget = 0;
  for(int n = 0 ; n < mapNVars ; n++) {
    mapVars.push_back(ntarget);
    ntarget++;
    if (haveMHD) {
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
    } // MHD
  }

  this->mapVars = idfx::ConvertVectorToIdefixArray(mapVars);

  this->bufferSize = data->np_tot[IDIR] * data->nghost[JDIR] * data->np_int[KDIR] * mapNVars;
  if (haveMHD) {
    // IDIR
    bufferSize += (data->np_tot[IDIR]+1) * data->nghost[JDIR] * data->np_int[KDIR];
    #if DIMENSIONS==3
    bufferSize += data->np_tot[IDIR] * data->nghost[JDIR] * (data->np_int[KDIR]+1);
    #endif  // DIMENSIONS
  }

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
