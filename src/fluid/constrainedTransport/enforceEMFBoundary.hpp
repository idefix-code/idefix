// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef FLUID_CONSTRAINEDTRANSPORT_ENFORCEEMFBOUNDARY_HPP_
#define FLUID_CONSTRAINEDTRANSPORT_ENFORCEEMFBOUNDARY_HPP_

#include "constrainedTransport.hpp"
#include "fluid.hpp"
#include "dataBlock.hpp"
#include "axis.hpp"

// In some cases, to avoid the divergence of the normal field of cells that should be identical
// in periodic boundary conditions or with domain decomposition, it is needed to average or
// periodise the EMFs of corresponding cells, as proposed by Stone et al. 2020 (p4, top left col.)
// This is because in some specific cases involving curvilinear coordinates, roundoff errors
// can accumulate, leading to a small drift of face-values that should be strictly equal.
// This behaviour is enabled using the flag below.

//#define ENFORCE_EMF_CONSISTENCY

template<typename Phys>
void ConstrainedTransport<Phys>::EnforceEMFBoundary() {
  idfx::pushRegion("Emf::EnforceEMFBoundary");
#if MHD == YES
  if(data->hydro->haveEmfBoundary)
    this->data->hydro->emfBoundaryFunc(*data, data->t);

  if(this->data->hydro->haveAxis) {
    this->data->hydro->boundary->axis->RegularizeEMFs();
  }

  #ifdef ENFORCE_EMF_CONSISTENCY
    #ifdef WITH_MPI
      // This average the EMFs at the domain surface with immediate neighbours
      // to ensure the EMFs exactly match
      this->ExchangeAll();
    #endif
  #endif

  IdefixArray3D<real> ex = this->ex;
  IdefixArray3D<real> ey = this->ey;
  IdefixArray3D<real> ez = this->ez;

  // Enforce specific EMF regularisation
  for(int dir=0 ; dir < DIMENSIONS ; dir++ ) {
    if(data->lbound[dir] == shearingbox || data->rbound[dir] == shearingbox) {
      SymmetrizeEMFShearingBox();
    }
    #ifdef ENFORCE_EMF_CONSISTENCY
      if(data->lbound[dir] == periodic && data->rbound[dir] == periodic) {
        // If domain decomposed, periodicity is already enforced by ExchangeAll
        if(data->mygrid->nproc[dir] == 1) {
          int ioffset = (dir == IDIR) ? data->np_int[IDIR] : 0;
          int joffset = (dir == JDIR) ? data->np_int[JDIR] : 0;
          int koffset = (dir == KDIR) ? data->np_int[KDIR] : 0;

          int ibeg = (dir == IDIR) ? data->beg[IDIR] : 0;
          int iend = (dir == IDIR) ? data->beg[IDIR]+1 : data->np_tot[IDIR];
          int jbeg = (dir == JDIR) ? data->beg[JDIR] : 0;
          int jend = (dir == JDIR) ? data->beg[JDIR]+1 : data->np_tot[JDIR];
          int kbeg = (dir == KDIR) ? data->beg[KDIR] : 0;
          int kend = (dir == KDIR) ? data->beg[KDIR]+1 : data->np_tot[KDIR];
          idefix_for("BoundaryEMFPeriodic",kbeg,kend,jbeg,jend,ibeg,iend,
            KOKKOS_LAMBDA (int k, int j, int i) {
              real em;

              if(dir==IDIR) {
                em = HALF_F*(ez(k,j,i)+ez(k,j,i+ioffset));
                ez(k,j,i) = em;
                ez(k,j,i+ioffset) = em;

                #if DIMENSIONS == 3
                em = HALF_F*(ey(k,j,i)+ey(k,j,i+ioffset));
                ey(k,j,i) = em;
                ey(k,j,i+ioffset) = em;
                #endif
              }

              if(dir==JDIR) {
                em = HALF_F*(ez(k,j,i)+ez(k,j+joffset,i));
                ez(k,j,i) = em;
                ez(k,j+joffset,i) = em;

                #if DIMENSIONS == 3
                em = HALF_F*(ex(k,j,i)+ex(k,j+joffset,i));
                ex(k,j,i) = em;
                ex(k,j+joffset,i) = em;
                #endif
              }

              if(dir==KDIR) {
                em = HALF_F*(ex(k,j,i)+ex(k+koffset,j,i));
                ex(k,j,i) = em;
                ex(k+koffset,j,i) = em;

                em = HALF_F*(ey(k,j,i)+ey(k+koffset,j,i));
                ey(k,j,i) = em;
                ey(k+koffset,j,i) = em;
              }
            });
        }
      }
    #endif //ENFORCE_EMF_CONSISTENCY
  }
#endif // MHD==YES
  idfx::popRegion();
}

template<typename Phys>
void ConstrainedTransport<Phys>::SymmetrizeEMFShearingBox() {
  idfx::pushRegion("Emf::EnforceEMFBoundary");
  #if MHD == YES

    IdefixArray2D<real> sbEyL = this->sbEyL;
    IdefixArray2D<real> sbEyR = this->sbEyR;
    IdefixArray2D<real> sbEyRL = this->sbEyRL;

    IdefixArray3D<real> ey = this->ey;

    // Store emf components on the left and right
    if(data->lbound[IDIR]==shearingbox) {
      int i = data->beg[IDIR];
      idefix_for("StoreEyL", 0, data->np_tot[KDIR], 0, data->np_tot[JDIR],
                            KOKKOS_LAMBDA(int k,int j) {
                              sbEyL(k,j) = ey(k,j,i);
                            });
    }
    if(data->rbound[IDIR]==shearingbox) {
      int i = data->end[IDIR];
      idefix_for("StoreEyL", 0, data->np_tot[KDIR], 0, data->np_tot[JDIR],
                            KOKKOS_LAMBDA(int k, int j) {
                              sbEyR(k,j) = ey(k,j,i);
                            });
    }

    // todo: exchange sbEyR and sbEyL when domain-decomposed along y
    #ifdef WITH_MPI
      if(data->mygrid->nproc[IDIR]>1) {
        int procLeft, procRight;
        const int size = data->np_tot[KDIR]*data->np_tot[JDIR];
        MPI_Status status;

        MPI_SAFE_CALL(MPI_Cart_shift(data->mygrid->CartComm,0,1,&procLeft,&procRight));
        if(data->lbound[IDIR]==shearingbox) {
          // We send to our left (which, by periodicity, is the right end of the domain)
          // our value of sbEyL and get
          MPI_Sendrecv(sbEyL.data(), size, realMPI, procLeft, 2001,
                       sbEyR.data(), size, realMPI, procLeft, 2002,
                       data->mygrid->CartComm, &status );
        }
        if(data->rbound[IDIR]==shearingbox) {
          // We send to our right (which, by periodicity, is the left end (=beginning)
          // of the domain) our value of sbEyR and get sbEyL
          MPI_Sendrecv(sbEyR.data(), size, realMPI, procRight, 2002,
                       sbEyL.data(), size, realMPI, procRight, 2001,
                       data->mygrid->CartComm, &status );
        }
      }
    #endif
    // Extrapolate on the left
    if(data->lbound[IDIR]==shearingbox) {
      int i = data->beg[IDIR];
      ExtrapolateEMFShearingBox(left, sbEyR, sbEyRL);

      idefix_for("ReplaceEyLeft", 0, data->np_tot[KDIR], 0, data->np_tot[JDIR],
                              KOKKOS_LAMBDA(int k,int j) {
                                ey(k,j,i) = 0.5*(sbEyL(k,j)+sbEyRL(k,j));
                              });
    }

    // Extrapolate on the right
    if(data->rbound[IDIR]==shearingbox) {
      int i = data->end[IDIR];
      ExtrapolateEMFShearingBox(right, sbEyL, sbEyRL);

      idefix_for("ReplaceEyRight", 0, data->np_tot[KDIR], 0, data->np_tot[JDIR],
                              KOKKOS_LAMBDA(int k,int j) {
                                ey(k,j,i) = 0.5*(sbEyR(k,j)+sbEyRL(k,j));
                              });
    }

  #endif
  idfx::popRegion();
}


template<typename Phys>
void ConstrainedTransport<Phys>::ExtrapolateEMFShearingBox(BoundarySide side,
                                                   IdefixArray2D<real> Ein,
                                                   IdefixArray2D<real> Eout) {
  const int nxi = data->np_int[IDIR];
  const int nxj = data->np_int[JDIR];

  const int ighost = data->nghost[IDIR];
  const int jghost = data->nghost[JDIR];

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
  real dL = std::fmod(sbVelocity*data->t,Ly);

  // translate this into # of cells
  const int m = static_cast<int> (std::floor(dL/dy+HALF_F));

  // remainding shift
  const real eps = dL / dy - m;

  // New we need to perform the shift
  idefix_for("BoundaryShearingBoxEMF", 0, data->np_tot[KDIR],
                                       0, data->np_tot[JDIR],
        KOKKOS_LAMBDA (int k, int j) {
          // jorigin
          const int jo = jghost + ((j-m-jghost)%nxj+nxj)%nxj;
          const int jop2 = jghost + ((jo+2-jghost)%nxj+nxj)%nxj;
          const int jop1 = jghost + ((jo+1-jghost)%nxj+nxj)%nxj;
          const int jom1 = jghost + ((jo-1-jghost)%nxj+nxj)%nxj;
          const int jom2 = jghost + ((jo-2-jghost)%nxj+nxj)%nxj;

          // Define Left and right fluxes
          // Fluxes are defined from slope-limited interpolation
          // Using Van-leer slope limiter (consistently with the main advection scheme)
          real Fl,Fr;
          real dqm, dqp, dq;

          if(eps>=ZERO_F) {
            // Compute Fl
            dqm = Ein(k,jom1) - Ein(k,jom2);
            dqp = Ein(k,jo) - Ein(k,jom1);
            dq = (dqp*dqm > ZERO_F ? TWO_F*dqp*dqm/(dqp + dqm) : ZERO_F);

            Fl = Ein(k,jom1) + 0.5*dq*(1.0-eps);
            //Compute Fr
            dqm=dqp;
            dqp = Ein(k,jop1) - Ein(k,jo);
            dq = (dqp*dqm > ZERO_F ? TWO_F*dqp*dqm/(dqp + dqm) : ZERO_F);

            Fr = Ein(k,jo) + 0.5*dq*(1.0-eps);
          } else {
            //Compute Fl
            dqm = Ein(k,jo) - Ein(k,jom1);
            dqp = Ein(k,jop1) - Ein(k,jo);
            dq = (dqp*dqm > ZERO_F ? TWO_F*dqp*dqm/(dqp + dqm) : ZERO_F);

            Fl = Ein(k,jo) - 0.5*dq*(1.0+eps);
            // Compute Fr
            dqm=dqp;
            dqp = Ein(k,jop2) - Ein(k,jop1);
            dq = (dqp*dqm > ZERO_F ? TWO_F*dqp*dqm/(dqp + dqm) : ZERO_F);

            Fr = Ein(k,jop1) - 0.5*dq*(1.0+eps);
          }
          Eout(k,j) = Ein(k,jo) - eps*(Fr - Fl);
        });
}
#endif // FLUID_CONSTRAINEDTRANSPORT_ENFORCEEMFBOUNDARY_HPP_
