// ***********************************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************************

#include "../idefix.hpp"
#include "hydro.hpp"

// Compute Corner EMFs from the one stored in the Riemann step
void Hydro::CalcCornerEMF(real t) {
  idfx::pushRegion("Hydro::CalcCornerEMF");

  // Corned EMFs
  IdefixArray3D<real> ex = this->emf.ex;
  IdefixArray3D<real> ey = this->emf.ey;
  IdefixArray3D<real> ez = this->emf.ez;

  // Face-centered EMFs
  IdefixArray3D<real> exj = this->emf.exj;
  IdefixArray3D<real> exk = this->emf.exk;
  IdefixArray3D<real> eyi = this->emf.eyi;
  IdefixArray3D<real> eyk = this->emf.eyk;
  IdefixArray3D<real> ezi = this->emf.ezi;
  IdefixArray3D<real> ezj = this->emf.ezj;

#if MHD == YES && DIMENSIONS >= 2

  #if EMF_AVERAGE == ARITHMETIC
  idefix_for("CalcCornerEMF",
             data->beg[KDIR],data->end[KDIR]+KOFFSET,
             data->beg[JDIR],data->end[JDIR]+JOFFSET,
             data->beg[IDIR],data->end[IDIR]+IOFFSET,
    KOKKOS_LAMBDA (int k, int j, int i) {
                  // CT_EMF_ArithmeticAverage (emf, 0.25);
      real w = ONE_FOURTH_F;
    #if DIMENSIONS == 3
      ex(k,j,i) = w * (exj(k,j,i) + exj(k-1,j,i) + exk(k,j,i) + exk(k,j-1,i));
      ey(k,j,i) = w * (eyi(k,j,i) + eyi(k-1,j,i) + eyk(k,j,i) + eyk(k,j,i-1));
    #endif
    #if DIMENSIONS >= 2
      ez(k,j,i) = w * (ezi(k,j,i) + ezi(k,j-1,i) + ezj(k,j,i) + ezj(k,j,i-1));
    #else
      ez(k,j,i) = w * (TWO_F*ezi(k,j,i) + ezj(k,j,i) + ezj(k,j,i-1));
    #endif
    }
  );
  #endif

  #if EMF_AVERAGE == UCT_CONTACT || EMF_AVERAGE == UCT0
  // 0. Compute cell-centered emf.
  IdefixArray4D<real> Vc = this->Vc;

  IdefixArray3D<real> Ex1 = this->emf.Ex1;
  IdefixArray3D<real> Ex2 = this->emf.Ex2;
  IdefixArray3D<real> Ex3 = this->emf.Ex3;

  // Required by Hall effect
  IdefixArray4D<real> J = this->J;
  real xHConstant = this->xH;
  ParabolicType haveHall = this->haveHall;
  IdefixArray3D<real> xHallArr = this->xHall;

  idefix_for("CalcCenterEMF",
             0,data->np_tot[KDIR],
             0,data->np_tot[JDIR],
             0,data->np_tot[IDIR],
    KOKKOS_LAMBDA (int k, int j, int i) {
      real vx1, vx2, vx3;
      real Bx1, Bx2, Bx3;

      vx1 = vx2 = vx3 = ZERO_F;
      Bx1 = Bx2 = Bx3 = ZERO_F;

      EXPAND( vx1 = Vc(VX1,k,j,i);  ,
              vx2 = Vc(VX2,k,j,i);  ,
              vx3 = Vc(VX3,k,j,i);  )

      EXPAND( Bx1 = Vc(BX1,k,j,i);  ,
              Bx2 = Vc(BX2,k,j,i);  ,
              Bx3 = Vc(BX3,k,j,i);  )

      // -- Compute inductive electric field

    #if DIMENSIONS == 3
      Ex1(k,j,i) = (vx3*Bx2 - vx2*Bx3);
      Ex2(k,j,i) = (vx1*Bx3 - vx3*Bx1);
    #endif
      Ex3(k,j,i) = (vx2*Bx1 - vx1*Bx2);

      // Compute Hall-induced centered EMF when needed
      if(haveHall) {
        // We store Jx1, Jx2, Jx3 in vx1, vx2 and vx3
        int ip1,jp1,kp1;

        ip1=i+1;
    #if DIMENSIONS >=2
        jp1 = j+1;
    #else
        jp1=j;
    #endif

    #if DIMENSIONS == 3
        kp1 = k+1;
    #else
        kp1 = k;
    #endif

        vx1 = AVERAGE_4D_YZ(J, IDIR, kp1, jp1, i);
        vx2 = AVERAGE_4D_XZ(J, JDIR, kp1, j, ip1);
        vx3 = AVERAGE_4D_XY(J, KDIR, k, jp1, ip1);

        if(haveHall == UserDefFunction)
          xH = xHallArr(k,j,i);
        else
          xH = xHConstant;

    #if DIMENSIONS == 3
        Ex1(k,j,i) += - xH*(vx3*Bx2 - vx2*Bx3);
        Ex2(k,j,i) += - xH*(vx1*Bx3 - vx3*Bx1);
    #endif
        Ex3(k,j,i) += - xH*(vx2*Bx1 - vx1*Bx2);
      }
    }
  );
  #endif // EMF_AVERAGE

  // 1. averaging scheme
  #if EMF_AVERAGE == UCT0
  idefix_for("CalcCornerEMF",
             data->beg[KDIR]-KOFFSET,data->end[KDIR]+KOFFSET,
             data->beg[JDIR]-JOFFSET,data->end[JDIR]+JOFFSET,
             data->beg[IDIR]-IOFFSET,data->end[IDIR]+IOFFSET,
    KOKKOS_LAMBDA (int k, int j, int i) {
    #if DIMENSIONS == 3
      exj(k,j,i) *= TWO_F;
      exk(k,j,i) *= TWO_F;
      eyi(k,j,i) *= TWO_F;
      eyk(k,j,i) *= TWO_F;

      exj(k,j,i) -= HALF_F*(Ex1(k,j-1,i) + Ex1(k,j,i));
      exk(k,j,i) -= HALF_F*(Ex1(k-1,j,i) + Ex1(k,j,i));

      eyi(k,j,i) -= HALF_F*(Ex2(k,j,i-1) + Ex2(k,j,i));
      eyk(k,j,i) -= HALF_F*(Ex2(k-1,j,i) + Ex2(k,j,i));
    #endif
      ezi(k,j,i) *= TWO_F;
      ezj(k,j,i) *= TWO_F;
      ezi(k,j,i) -= HALF_F*(Ex3(k,j,i-1) + Ex3(k,j,i));
      ezj(k,j,i) -= HALF_F*(Ex3(k,j-1,i) + Ex3(k,j,i));
    }
  );

  idefix_for("CalcCornerEMF",
             data->beg[KDIR],data->end[KDIR]+KOFFSET,
             data->beg[JDIR],data->end[JDIR]+JOFFSET,
             data->beg[IDIR],data->end[IDIR]+IOFFSET,
    KOKKOS_LAMBDA (int k, int j, int i) {
      // CT_EMF_ArithmeticAverage (emf, 0.25);
      real w = ONE_FOURTH_F;
    #if DIMENSIONS == 3
      ex(k,j,i) = w * (exj(k,j,i) + exj(k-1,j,i) + exk(k,j,i) + exk(k,j-1,i));
      ey(k,j,i) = w * (eyi(k,j,i) + eyi(k-1,j,i) + eyk(k,j,i) + eyk(k,j,i-1));
    #endif
    #if DIMENSIONS >= 2
      ez(k,j,i) = w * (ezi(k,j,i) + ezi(k,j-1,i) + ezj(k,j,i) + ezj(k,j,i-1));
    #else
      ez(k,j,i) = w * (TWO_F*ezi(k,j,i) + ezj(k,j,i) + ezj(k,j,i-1));
    #endif
    }
  );
  #endif // EMF_AVERAGE

  #if EMF_AVERAGE == UCT_CONTACT
  IdefixArray3D<int> svx = this->emf.svx;
  IdefixArray3D<int> svy = this->emf.svy;
  IdefixArray3D<int> svz = this->emf.svz;

  idefix_for("EMF_ArithmeticAverage",
             data->beg[KDIR],data->end[KDIR]+KOFFSET,
             data->beg[JDIR],data->end[JDIR]+JOFFSET,
             data->beg[IDIR],data->end[IDIR]+IOFFSET,
    KOKKOS_LAMBDA (int k, int j, int i) {
      // CT_EMF_ArithmeticAverage (emf, 1.0);
      real w = ONE_F;
    #if DIMENSIONS == 3
      ex(k,j,i) = w * (exj(k,j,i) + exj(k-1,j,i) + exk(k,j,i) + exk(k,j-1,i));
      ey(k,j,i) = w * (eyi(k,j,i) + eyi(k-1,j,i) + eyk(k,j,i) + eyk(k,j,i-1));
    #endif
    #if DIMENSIONS >= 2
      ez(k,j,i) = w * (ezi(k,j,i) + ezi(k,j-1,i) + ezj(k,j,i) + ezj(k,j,i-1));
    #else
      ez(k,j,i) = w * (TWO_F*ezi(k,j,i) + ezj(k,j,i) + ezj(k,j,i-1));
    #endif

      //CT_EMF_IntegrateToCorner (data, emf, grid);
      int iu, ju, ku;
      D_EXPAND( int sx;  ,
                int sy;  ,
                int sz;  )

      D_EXPAND( sx = svx(k,j,i);  ,
                sy = svy(k,j,i);  ,
                sz = svz(k,j,i);  )

      D_EXPAND( iu = sx > 0 ? i-1:i;  ,  // -- upwind index
                ju = sy > 0 ? j-1:j;  ,
                ku = sz > 0 ? k-1:k;  )

      // Span X - Faces:    dEz/dy, dEy/dz

      if (sx == 0) {
        ez(k,j,i) += HALF_F*(ezj(k,j,i-1) - Ex3(k,j-1,i-1) + ezj(k,j,i) - Ex3(k,j-1,i));
        ez(k,j,i) -= HALF_F*(Ex3(k,j,i-1) - ezj(k,j,i-1)   + Ex3(k,j,i) - ezj(k,j,i));
    #if DIMENSIONS == 3
        ey(k,j,i) += HALF_F*(eyk(k,j,i-1) - Ex2(k-1,j,i-1) + eyk(k,j,i) - Ex2(k-1,j,i));
        ey(k,j,i) -= HALF_F*(Ex2(k,j,i-1) - eyk(k,j,i-1)   + Ex2(k,j,i) - eyk(k,j,i));
    #endif
      } else {
        ez(k,j,i) += ezj(k,j,iu) - Ex3(k,j-1,iu);
        ez(k,j,i) -= Ex3(k,j,iu) - ezj(k,j,iu);
    #if DIMENSIONS == 3
        ey(k,j,i) += eyk(k,j,iu) - Ex2(k-1,j,iu);
        ey(k,j,i) -= Ex2(k,j,iu) - eyk(k,j,iu);
    #endif
      }

      // Span Y - Faces:    dEz/dx, dEx/dz

      if (sy == 0) {
        ez(k,j,i) += HALF_F*(ezi(k,j-1,i) - Ex3(k,j-1,i-1) + ezi(k,j,i) - Ex3(k,j,i-1));
        ez(k,j,i) -= HALF_F*(Ex3(k,j-1,i) - ezi(k,j-1,i)   + Ex3(k,j,i) - ezi(k,j,i));
    #if DIMENSIONS == 3
        ex(k,j,i) += HALF_F*(exk(k,j-1,i) - Ex1(k-1,j-1,i) + exk(k,j,i) - Ex1(k-1,j,i));
        ex(k,j,i) -= HALF_F*(Ex1(k,j-1,i) - exk(k,j-1,i)   + Ex1(k,j,i) - exk(k,j,i));
    #endif
      } else {
        ez(k,j,i) += ezi(k,ju,i) - Ex3(k,ju,i-1);
        ez(k,j,i) -= Ex3(k,ju,i) - ezi(k,ju,i);
    #if DIMENSIONS == 3
        ex(k,j,i) += exk(k,ju,i) - Ex1(k-1,ju,i);
        ex(k,j,i) -= Ex1(k,ju,i) - exk(k,ju,i);
    #endif
      }

      // Span Z - Faces:    dEx/dy, dEy/dx

    #if DIMENSIONS == 3
      if (sz == 0) {
        ex(k,j,i) += HALF_F*(exj(k-1,j,i) - Ex1(k-1,j-1,i) + exj(k,j,i) - Ex1(k,j-1,i));
        ex(k,j,i) -= HALF_F*(Ex1(k-1,j,i) - exj(k-1,j,i)   + Ex1(k,j,i) - exj(k,j,i));
        ey(k,j,i) += HALF_F*(eyi(k-1,j,i) - Ex2(k-1,j,i-1) + eyi(k,j,i) - Ex2(k,j,i-1));
        ey(k,j,i) -= HALF_F*(Ex2(k-1,j,i) - eyi(k-1,j,i)   + Ex2(k,j,i) - eyi(k,j,i));
      } else {
        ex(k,j,i) += exj(ku,j,i) - Ex1(ku,j-1,i);
        ex(k,j,i) -= Ex1(ku,j,i) - exj(ku,j,i);
        ey(k,j,i) += eyi(ku,j,i) - Ex2(ku,j,i-1);
        ey(k,j,i) -= Ex2(ku,j,i) - eyi(ku,j,i);
      }

      ex(k,j,i) *= ONE_FOURTH_F;
      ey(k,j,i) *= ONE_FOURTH_F;
    #endif
      ez(k,j,i) *= ONE_FOURTH_F;
    }
  );
  #endif // EMF_AVERAGE

  if(haveEmfBoundary)
    emfBoundaryFunc(*data, t);

#endif // MHD

  idfx::popRegion();
}
