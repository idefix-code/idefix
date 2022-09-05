// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2022 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include "hydro.hpp"
#include "dataBlock.hpp"

KOKKOS_INLINE_FUNCTION int Kmax(int n, int m) {
  return( n>m ? n : m);
}
// This function coarsen the flow according to the grid coarsening array

void Hydro::CoarsenFlow(IdefixArray4D<real> &Vi) {
  idfx::pushRegion("Hydro::CoarsenFlow");

  IdefixArray3D<real> dV   = data->dV;
  for(int dir = 0 ; dir < DIMENSIONS ; dir++) {
    if(!data->coarseningDirection[dir]) continue;
    int begDir = data->beg[dir];
    int endDir = data->end[dir];

    IdefixArray2D<int> coarseningLevel = data->coarseningLevel[dir];

    idefix_for("Hydro_CoarsenFlow",
          0, NVAR,
          data->beg[KDIR],data->end[KDIR],
          data->beg[JDIR],data->end[JDIR],
          data->beg[IDIR],data->end[IDIR],
      KOKKOS_LAMBDA (int n, int k, int j, int i) {
        int factor, index;
        int ioffset = 0;
        int joffset = 0;
        int koffset = 0;
        //factor = 2^(coarsening-1)
        if(dir==IDIR) {
          factor = 1 << (coarseningLevel(k,j) - 1);
          index = i;
          ioffset = 1;
        }
        if(dir==JDIR) {
          factor = 1 << (coarseningLevel(k,i) - 1);
          index = j;
          joffset = 1;
        }
        if(dir==KDIR) {
          factor = 1 << (coarseningLevel(j,i) - 1);
          index = k;
          koffset = 1;
        }
        // check if coarsening is required in this region
        if(factor>1) {
          // We average the cells by groups of "factor" in direction "dir",
          // so only the first element of each group will do the job.
          if( (index-begDir)%factor == 0) {
            real q = 0.0;
            real V = 0.0;
            for(int shift = 0 ; shift < factor ; shift++) {
              q = q + Vi(n, k + shift*koffset, j + shift*joffset, i+shift*ioffset)
                    * dV(k + shift*koffset, j + shift*joffset, i+shift*ioffset);
              V = V + dV(k + shift*koffset, j + shift*joffset, i+shift*ioffset);
            }
            // Average
            q = q/V;
            // Write back the cell elements
            for(int shift = 0 ; shift < factor ; shift++) {
              Vi(n, k + shift*koffset, j + shift*joffset, i+shift*ioffset) = q;
            }
          }
        }
    });
  }
  idfx::popRegion();
}

void Hydro::CoarsenMagField(IdefixArray4D<real> &Vsin) {
#if MHD == YES
  idfx::pushRegion("Hydro::CoarsenMagField");
  #if DIMENSIONS >= 2
  /**********************************
   * MHD Part                       *
   * ********************************/
  for(int dir = 0 ; dir < DIMENSIONS ; dir++) {
    if(!data->coarseningDirection[dir]) continue;
    int begDir = data->beg[dir];
    int endDir = data->end[dir];

    IdefixArray2D<int> coarseningLevel = data->coarseningLevel[dir];

    const int BXn = dir;
    const int BXt = (dir == IDIR ? BX2s : BX1s);
    const int BXb = (dir == KDIR ? BX2s : BX3s);

    IdefixArray3D<real> An = data->A[dir];
    IdefixArray3D<real> At = data->A[BXt];
    IdefixArray3D<real> Ab = data->A[BXb];

    int it = 0, ib = 0;
    int jt = 0, jb = 0;
    int kt = 0, kb = 0;
    int endk = data->end[KDIR];
    int endj = data->end[JDIR];
    int endi = data->end[IDIR];

    // find the index over which the faces are off-centered
    // (trick so that we can write a single loop)
    if(dir==IDIR) {
      jt=1;
      kb=1;
      endj++;
      #if DIMENSIONS == 3
      endk++;
      #endif
    }
    if(dir==JDIR) {
      it=1;
      kb=1;
      endi++;
      #if DIMENSIONS == 3
      endk++;
      #endif
    }
    if(dir==KDIR) {
      it=1;
      jb=1;
      endi++;
      endj++;
    }
    #if DIMENSIONS < 3
    // force kb to 0
    kb = 0;
    #endif

    // Coarsen components normal to the coarsening direction first (BXt and BXb)
    idefix_for("Hydro_CoarsenFlow_BXsn",
        data->beg[KDIR],endk,
        data->beg[JDIR],endj,
        data->beg[IDIR],endi,
    KOKKOS_LAMBDA (int k, int j, int i) {
      int coarsening_t =1;
      [[maybe_unused]] int coarsening_b = 1;
      int index;
      int ioffset = 0;
      int joffset = 0;
      int koffset = 0;
      if(dir==IDIR) {
        coarsening_t = Kmax(coarseningLevel(k,j-1), coarseningLevel(k,j));
        #if DIMENSIONS == 3
        coarsening_b = Kmax(coarseningLevel(k-1,j), coarseningLevel(k,j));
        #endif
        ioffset = 1;
        index=i;
      }
      if(dir==JDIR) {
        coarsening_t = Kmax(coarseningLevel(k,i-1), coarseningLevel(k,i));
        #if DIMENSIONS == 3
        coarsening_b = Kmax(coarseningLevel(k-1,i), coarseningLevel(k,i));
        #endif
        joffset = 1;
        index=j;
      }
      if(dir==KDIR) {
        coarsening_t = Kmax(coarseningLevel(j,i-1), coarseningLevel(j,i));
        coarsening_b = Kmax(coarseningLevel(j-1,i), coarseningLevel(j,i));
        koffset = 1;
        index=k;
      }

      //coarsening_t=0;
      //coarsening_b=0;
      // Treat t component
      if(coarsening_t>1) {
        int factor_t = 1 << (coarsening_t-1);
        // We average the cells by groups of "factor" in direction "dir",
        // so only the first element of each group will do the job.
        if( (index-begDir)%factor_t == 0) {
          real q = 0.0;
          real A = 0.0;
          for(int shift = 0 ; shift < factor_t ; shift++) {
            q = q + Vsin(BXt, k + shift*koffset, j + shift*joffset, i+shift*ioffset)
                  * At(k + shift*koffset, j + shift*joffset, i+shift*ioffset);
            A = A + At(k + shift*koffset, j + shift*joffset, i+shift*ioffset);
          }

          // If the Area is zero, do a point average instead (this happens on the axis instance)
          if(FABS(A) < 1e-10) {
            // do a point average instead of a surface average
            q=0.0;
            A=0.0;
            for(int shift = 0 ; shift < factor_t ; shift++) {
              q = q + Vsin(BXt, k + shift*koffset, j + shift*joffset, i+shift*ioffset);
              A = A + 1.0;
            }
          }
          q = q/A;
          // Write back the cell elements
          for(int shift = 0 ; shift < factor_t ; shift++) {
            Vsin(BXt, k + shift*koffset, j + shift*joffset, i+shift*ioffset) = q;
          }
        }
      }
      // Treat b component
      #if DIMENSIONS == 3
      if(coarsening_b>1) {
        int factor_b = 1 << (coarsening_b-1);
        // We average the cells by groups of "factor" in direction "dir",
        // so only the first element of each group will do the job.
        if( (index-begDir)%factor_b == 0) {
          real q = 0.0;
          real A = 0.0;
          for(int shift = 0 ; shift < factor_b ; shift++) {
            q = q + Vsin(BXb, k + shift*koffset, j + shift*joffset, i+shift*ioffset)
                  * Ab(k + shift*koffset, j + shift*joffset, i+shift*ioffset);
            A = A + Ab(k + shift*koffset, j + shift*joffset, i+shift*ioffset);
          }
          // If the Area is zero, do a point average instead (this happens on the axis instance)
          if(FABS(A) < 1e-10) {
            // do a point average instead of a surface average
            q=0.0;
            A=0.0;
            for(int shift = 0 ; shift < factor_b ; shift++) {
              q = q + Vsin(BXb, k + shift*koffset, j + shift*joffset, i+shift*ioffset);
              A = A + 1.0;
            }
          }
          // Average
          q = q/A;
          // Write back the cell elements
          for(int shift = 0 ; shift < factor_b ; shift++) {
            Vsin(BXb, k + shift*koffset, j + shift*joffset, i+shift*ioffset) = q;
          }
        }
      }
      #endif
    });

      // Coarsen components parralel the coarsening direction (BXn)
    idefix_for("Hydro_CoarsenFlow_BXsn",
        data->beg[KDIR],data->end[KDIR],
        data->beg[JDIR],data->end[JDIR],
        data->beg[IDIR],data->end[IDIR],
    KOKKOS_LAMBDA (int k, int j, int i) {
      int coarsening;
      int index;
      int ioffset = 0;
      int joffset = 0;
      int koffset = 0;
      // Pick the highest coarsening level of the current coordinates and its immediate neighbours
      if(dir==IDIR) {
        coarsening = coarseningLevel(k,j);
        coarsening = Kmax(coarsening, coarseningLevel(k,j-1));
        coarsening = Kmax(coarsening, coarseningLevel(k,j+1));
        #if DIMENSIONS == 3
        coarsening = Kmax(coarsening, coarseningLevel(k-1,j));
        coarsening = Kmax(coarsening, coarseningLevel(k+1,j));
        #endif
        ioffset = 1;
        index=i;
      }
      if(dir==JDIR) {
        coarsening = coarseningLevel(k,i);
        coarsening = Kmax(coarsening, coarseningLevel(k,i-1));
        coarsening = Kmax(coarsening, coarseningLevel(k,i+1));
        #if DIMENSIONS == 3
        coarsening = Kmax(coarsening, coarseningLevel(k-1,i));
        coarsening = Kmax(coarsening, coarseningLevel(k+1,i));
        #endif
        joffset = 1;
        index=j;
      }
      if(dir==KDIR) {
        coarsening = coarseningLevel(j,i);
        coarsening = Kmax(coarsening, coarseningLevel(j,i-1));
        coarsening = Kmax(coarsening, coarseningLevel(j,i+1));
        coarsening = Kmax(coarsening, coarseningLevel(j-1,i));
        coarsening = Kmax(coarsening, coarseningLevel(j+1,i));
        koffset = 1;
        index=k;
      }

      if(coarsening > 1) {
        // the current cell had tangential field components which have been coarsened. Hence,
        // We reconstruct the parrallel field component with divB=0

        int factor = 1 << (coarsening-1);
        if( (index-begDir)%factor == 0) {
          for(int shift = 0 ; shift < factor-1 ; shift++) {
            real qt = Vsin(BXt,k+kt+shift*koffset, j+jt+shift*joffset, i+it+shift*ioffset)
                        * At(k+kt+shift*koffset, j+jt+shift*joffset, i+it+shift*ioffset)
                    - Vsin(BXt,k+shift*koffset, j+shift*joffset, i+shift*ioffset)
                        *  At(k+shift*koffset, j+shift*joffset, i+shift*ioffset);
            real qb = Vsin(BXb,k+kb+shift*koffset, j+jb+shift*joffset, i+ib+shift*ioffset)
                        * Ab(k+kb+shift*koffset, j+jb+shift*joffset, i+ib+shift*ioffset)
                    - Vsin(BXb,k+shift*koffset, j+shift*joffset, i+shift*ioffset)
                        *  Ab(k+shift*koffset, j+shift*joffset, i+shift*ioffset);
/*
            if(i==10 && j == 3 && k == 10) {
              idfx::cout << "shift=" << shift << " ; oldBxn=" << Vsin(BXn, k+(shift+1)*koffset, j+(shift+1)*joffset, i+(shift+1)*ioffset) << std::endl;
              idfx::cout << "qt=" << qt << " ; qb=" << qb << std::endl;
            }*/

            Vsin(BXn, k+(shift+1)*koffset, j+(shift+1)*joffset, i+(shift+1)*ioffset) =
                1.0/An(k+(shift+1)*koffset, j+(shift+1)*joffset, i+(shift+1)*ioffset) *
                (
                  Vsin(BXn, k+shift*koffset, j+shift*joffset, i+shift*ioffset) *
                    An(k+shift*koffset, j+shift*joffset, i+shift*ioffset)
                  - qt - qb
                );
                /*
            if(i==10 && j == 3 && k == 10) {
              idfx::cout << "shift=" << shift << " ; newsBxn=" << Vsin(BXn, k+(shift+1)*koffset, j+(shift+1)*joffset, i+(shift+1)*ioffset) << std::endl;
            }*/
          }
        }
      }
    });
  }
  idfx::popRegion();
#endif
}

void Hydro::CoarsenVectorPotential() {
#if MHD == YES && defined(EVOLVE_VECTOR_POTENTIAL)
  idfx::pushRegion("Hydro::VectorPotential");

  IdefixArray4D<real> Ve = this->Ve;

  #if DIMENSIONS < 3
  IDEFIX_ERROR("Grid Coarsening with vector potential is not implemented for DIMENSIONS < 3");
  #endif

  for(int dir = 0 ; dir < DIMENSIONS ; dir++) {
    if(!data->coarseningDirection[dir]) continue;
    int begDir = data->beg[dir];
    int endDir = data->end[dir];

    const int AXn = dir;
    const int AXt = (dir == IDIR ? AX2e : AX1e);
    const int AXb = (dir == KDIR ? AX2e : AX3e);

    IdefixArray2D<int> coarseningLevel = data->coarseningLevel[dir];

    idefix_for("Hydro_CoarsenFlow",
          data->beg[KDIR],data->end[KDIR],
          data->beg[JDIR],data->end[JDIR],
          data->beg[IDIR],data->end[IDIR],
      KOKKOS_LAMBDA (int k, int j, int i) {
        int factor, index;
        int ioffset = 0;
        int joffset = 0;
        int koffset = 0;
        //factor = 2^(coarsening-1)
        if(dir==IDIR) {
          factor = 1 << (coarseningLevel(k,j) - 1);
          index = i;
          ioffset = 1;
        }
        if(dir==JDIR) {
          factor = 1 << (coarseningLevel(k,i) - 1);
          index = j;
          joffset = 1;
        }
        if(dir==KDIR) {
          factor = 1 << (coarseningLevel(j,i) - 1);
          index = k;
          koffset = 1;
        }
        // check if coarsening is required in this region
        if(factor>1) {
          // We average the cells by groups of "factor" in direction "dir",
          // so only the first element of each group will do the job.
          if( (index-begDir)%factor == 0) {
            real q = 0.0;
            for(int shift = 0 ; shift < factor ; shift++) {
              q = q + Ve(AXn, k + shift*koffset, j + shift*joffset, i+shift*ioffset);
            }
            // Average
            q = q/factor;
            // Write back the cell elements
            for(int shift = 0 ; shift < factor ; shift++) {
              Ve(AXn, k + shift*koffset, j + shift*joffset, i+shift*ioffset) = q;
              real a = ((real) shift) / ((real) factor);

              Ve(AXt, k + shift*koffset, j + shift*joffset, i+shift*ioffset) =
                  a * Ve(AXt, k + factor*koffset, j + factor*joffset, i+factor*ioffset)
                  + (1-a) * Ve(AXt, k, j, i);

              Ve(AXb, k + shift*koffset, j + shift*joffset, i+shift*ioffset) =
                  a * Ve(AXb, k + factor*koffset, j + factor*joffset, i+factor*ioffset)
                  + (1-a) * Ve(AXb, k, j, i);
            }
          }
        }
    });
  }

  #endif // DIMENSIONS >= 2
  idfx::popRegion();

  #endif // MHD
}
