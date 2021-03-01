// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2021 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include "electroMotiveForce.hpp"
#include "hydro.hpp"
#include "dataBlock.hpp"


void ElectroMotiveForce::EnforceEMFBoundary() {
  idfx::pushRegion("Emf::EnforceEMFBoundary");
#if MHD == YES
  if(data->hydro.haveEmfBoundary)
    this->data->hydro.emfBoundaryFunc(*data, data->t);

  if(this->data->hydro.haveAxis) {
    this->data->hydro.myAxis.SymmetrizeEx1();
  }

  #ifdef WITH_MPI
    // This average the EMFs at the domain surface with immediate neighbours
    // to ensure the EMFs exactly match
    this->ExchangeAll();
  #endif

  IdefixArray3D<real> ex = this->ex;
  IdefixArray3D<real> ey = this->ey;
  IdefixArray3D<real> ez = this->ez;

  // Ensure EMF are strictly periodic to avoid accumulation of roundoff errors
  for(int dir=0 ; dir < DIMENSIONS ; dir++ ) {
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
  }
#endif // MHD==YES
  idfx::popRegion();
}
