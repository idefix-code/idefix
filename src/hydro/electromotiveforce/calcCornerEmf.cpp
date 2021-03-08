// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2021 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include "hydro.hpp"
#include "dataBlock.hpp"

// Compute Corner EMFs from the one stored in the Riemann step
void ElectroMotiveForce::CalcCornerEMF(real t) {
  idfx::pushRegion("ElectroMotiveForce::CalcCornerEMF");

  // Corned EMFs
  IdefixArray3D<real> ex = this->ex;
  IdefixArray3D<real> ey = this->ey;
  IdefixArray3D<real> ez = this->ez;

  // Face-centered EMFs
  IdefixArray3D<real> exj = this->exj;
  IdefixArray3D<real> exk = this->exk;
  IdefixArray3D<real> eyi = this->eyi;
  IdefixArray3D<real> eyk = this->eyk;
  IdefixArray3D<real> ezi = this->ezi;
  IdefixArray3D<real> ezj = this->ezj;

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
  IdefixArray4D<real> Vc = hydro->Vc;

  IdefixArray3D<real> Ex1 = this->Ex1;
  IdefixArray3D<real> Ex2 = this->Ex2;
  IdefixArray3D<real> Ex3 = this->Ex3;

  // Required by Hall effect
  IdefixArray4D<real> J = hydro->J;
  real xHConstant = hydro->xH;
  HydroModuleStatus haveHall = hydro->haveHall;
  IdefixArray3D<real> xHallArr = hydro->xHall;

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


  IdefixArray3D<int> svx = this->svx;
  IdefixArray3D<int> svy = this->svy;
  IdefixArray3D<int> svz = this->svz;


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
    });

  // We define atomic references to ex,ey, ez because the loop "EMF_Integrate to corner"
  // Is not incrementing only e(k,j,i) but also its neighbour. Hence race conditions
  // Could occur in this loop.

  IdefixAtomicArray3D<real> exAtomic = this->ex;
  IdefixAtomicArray3D<real> eyAtomic = this->ey;
  IdefixAtomicArray3D<real> ezAtomic = this->ez;
  idefix_for("EMF_Integrate to Corner",
             data->beg[KDIR]-KOFFSET,data->end[KDIR]+KOFFSET,
             data->beg[JDIR]-JOFFSET,data->end[JDIR]+JOFFSET,
             data->beg[IDIR]-IOFFSET,data->end[IDIR]+IOFFSET,
    KOKKOS_LAMBDA (int k, int j, int i) {
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
      // Reminder:
      //   - Ez(k,j,i) is centered on  (i    ,j    ,k)
      //   - svx(k,j,i) is centered on (i-1/2,j    ,k)
      //   - ezj(k,j,i) is centered on (i    ,j-1/2,k)
      //   - ezi(k,j,i) is centered on (i-1/2,j    ,k)
      //   - ez(k,j,i) is centered on  (i-1/2,j-1/2,k)

#define DEX_DYP(k,j,i)    (exj(k,j+1,i) - Ex1(k,j,i))
#define DEX_DZP(k,j,i)    (exk(k+1,j,i) - Ex1(k,j,i))

#define DEY_DXP(k,j,i)    (eyi(k,j,i+1) - Ex2(k,j,i))
#define DEY_DZP(k,j,i)    (eyk(k+1,j,i) - Ex2(k,j,i))

#define DEZ_DXP(k,j,i)    (ezi(k,j,i+1) - Ex3(k,j,i))
#define DEZ_DYP(k,j,i)    (ezj(k,j+1,i) - Ex3(k,j,i))

#define DEX_DYM(k,j,i)    (Ex1(k,j,i) - exj(k,j,i))
#define DEX_DZM(k,j,i)    (Ex1(k,j,i) - exk(k,j,i))

#define DEY_DXM(k,j,i)    (Ex2(k,j,i) - eyi(k,j,i))
#define DEY_DZM(k,j,i)    (Ex2(k,j,i) - eyk(k,j,i))

#define DEZ_DXM(k,j,i)    (Ex3(k,j,i) - ezi(k,j,i))
#define DEZ_DYM(k,j,i)    (Ex3(k,j,i) - ezj(k,j,i))


      if (sx == 0) {
        ezAtomic(k,j+1,i) += HALF_F*(DEZ_DYP(k,j,i-1) + DEZ_DYP(k,j,i));
        ezAtomic(k,j  ,i) -= HALF_F*(DEZ_DYM(k,j,i-1) + DEZ_DYM(k,j,i));

    #if DIMENSIONS == 3
        eyAtomic(k+1,j,i) += HALF_F*(DEY_DZP(k,j,i-1) + DEY_DZP(k,j,i));
        eyAtomic(k  ,j,i) -= HALF_F*(DEY_DZM(k,j,i-1) + DEY_DZM(k,j,i));
    #endif
      } else {
        ezAtomic(k,j+1,i) += DEZ_DYP(k,j,iu);
        ezAtomic(k,j,  i) -= DEZ_DYM(k,j,iu);
    #if DIMENSIONS == 3
        eyAtomic(k+1,j,i) += DEY_DZP(k,j,iu);
        eyAtomic(k  ,j,i) -= DEY_DZM(k,j,iu);
    #endif
      }

      // Span Y - Faces:    dEz/dx, dEx/dz

      if (sy == 0) {
        ezAtomic(k,j,i+1) += HALF_F*(DEZ_DXP(k,j-1,i) + DEZ_DXP(k,j,i));
        ezAtomic(k,j,i  ) -= HALF_F*(DEZ_DXM(k,j-1,i) + DEZ_DXM(k,j,i));
    #if DIMENSIONS == 3
        exAtomic(k+1,j,i) += HALF_F*(DEX_DZP(k,j-1,i) + DEX_DZP(k,j,i));
        exAtomic(k  ,j,i) -= HALF_F*(DEX_DZM(k,j-1,i) + DEX_DZM(k,j,i));
    #endif
      } else {
        ezAtomic(k,j,i+1) += DEZ_DXP(k,ju,i);
        ezAtomic(k,j,i  ) -= DEZ_DXM(k,ju,i);
    #if DIMENSIONS == 3
        exAtomic(k+1,j,i) += DEX_DZP(k,ju,i);
        exAtomic(k  ,j,i) -= DEX_DZM(k,ju,i);
    #endif
      }

      // Span Z - Faces:    dEx/dy, dEy/dx

    #if DIMENSIONS == 3
      if (sz == 0) {
        exAtomic(k,j+1,i  ) += HALF_F*(DEX_DYP(k-1,j,i) + DEX_DYP(k,j,i));
        exAtomic(k,j  ,i  ) -= HALF_F*(DEX_DYM(k-1,j,i) + DEX_DYM(k,j,i));
        eyAtomic(k,j  ,i+1) += HALF_F*(DEY_DXP(k-1,j,i) + DEY_DXP(k,j,i));
        eyAtomic(k,j  ,i  ) -= HALF_F*(DEY_DXM(k-1,j,i) + DEY_DXM(k,j,i));
      } else {
        exAtomic(k,j+1,i  ) += DEX_DYP(ku,j,i);
        exAtomic(k,j  ,i  ) -= DEX_DYM(ku,j,i);
        eyAtomic(k,j  ,i+1) += DEY_DXP(ku,j,i);
        eyAtomic(k,j  ,i  ) -= DEY_DXM(ku,j,i);
      }
    #endif // DIMENSIONS
    });

    idefix_for("EMF_Renormalize",
             data->beg[KDIR],data->end[KDIR]+KOFFSET,
             data->beg[JDIR],data->end[JDIR]+JOFFSET,
             data->beg[IDIR],data->end[IDIR]+IOFFSET,
    KOKKOS_LAMBDA (int k, int j, int i) {
    #if DIMENSIONS ==  3
      ex(k,j,i) *= ONE_FOURTH_F;
      ey(k,j,i) *= ONE_FOURTH_F;
    #endif
      ez(k,j,i) *= ONE_FOURTH_F;
    }
  );

  #undef DEX_DYP
  #undef DEX_DZP
  #undef DEY_DXP
  #undef DEY_DZP
  #undef DEZ_DXP
  #undef DEZ_DYP
  #undef DEX_DYM
  #undef DEX_DZM
  #undef DEY_DXM
  #undef DEY_DZM
  #undef DEZ_DXM
  #undef DEZ_DYM

  #endif // EMF_AVERAGE

  #if EMF_AVERAGE == UCT_HLL
    calcRiemannEmf();
  #endif

#endif // MHD

  idfx::popRegion();
}
