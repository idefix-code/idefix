// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************
#ifndef FLUID_CONSTRAINEDTRANSPORT_CALCCORNEREMF_HPP_
#define FLUID_CONSTRAINEDTRANSPORT_CALCCORNEREMF_HPP_

#include "fluid.hpp"
#include "dataBlock.hpp"

// Compute Corner EMFs from the one stored in the Riemann step
template<typename Phys>
void ConstrainedTransport<Phys>::CalcCornerEMF(real t) {
  idfx::pushRegion("ConstrainedTransport::CalcCornerEMF");

#if MHD == YES && DIMENSIONS >= 2

  if(averaging==arithmetic) {
    CalcArithmeticAverage();
  }
  if(averaging==uct_contact||averaging==uct0) {
    CalcCellCenteredEMF();
    if(averaging==uct0) {
      CalcUCT0Average();
    } else {
      CalcContactAverage();
    }
  }
  if(averaging==uct_hll || averaging==uct_hlld) {
    CalcRiemannAverage();
  }

#endif // MHD

  idfx::popRegion();
}



// Compute Corner EMFs from arithmetic averages
template<typename Phys>
void ConstrainedTransport<Phys>::CalcArithmeticAverage() {
  idfx::pushRegion("ConstrainedTransport::CalcCornerEMF");

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

  idefix_for("CalcArithmeticAverage",
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
    });
#endif // MHD

  idfx::popRegion();
}

template<typename Phys>
void ConstrainedTransport<Phys>::CalcCellCenteredEMF() {
  idfx::pushRegion("ConstrainedTransport::CalcCellCenteredEMF");
  IdefixArray4D<real> Vc = hydro->Vc;
    // cell-centered EMFs
  IdefixArray3D<real> Ex1 = this->Ex1;
  IdefixArray3D<real> Ex2 = this->Ex2;
  IdefixArray3D<real> Ex3 = this->Ex3;

#if MHD == YES && DIMENSIONS >= 2
  idefix_for("CalcCenterEMF",
            0,data->np_tot[KDIR],
            0,data->np_tot[JDIR],
            0,data->np_tot[IDIR],
    KOKKOS_LAMBDA (int k, int j, int i) {
      [[maybe_unused]] real vx1, vx2, vx3;
      [[maybe_unused]] real Bx1, Bx2, Bx3;

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
    }
  );

#endif // MHD

  idfx::popRegion();
}

template<typename Phys>
void ConstrainedTransport<Phys>::CalcUCT0Average() {
  idfx::pushRegion("ConstrainedTransport::CalcUCT0Average");
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

  // cell-centered EMFs
  IdefixArray3D<real> Ex1 = this->Ex1;
  IdefixArray3D<real> Ex2 = this->Ex2;
  IdefixArray3D<real> Ex3 = this->Ex3;

#if MHD == YES && DIMENSIONS >= 2
  idefix_for("CalcUCT0FaceCentered",
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

  idefix_for("CalcUCT0CornerEMF",
            data->beg[KDIR],data->end[KDIR]+KOFFSET,
            data->beg[JDIR],data->end[JDIR]+JOFFSET,
            data->beg[IDIR],data->end[IDIR]+IOFFSET,
    KOKKOS_LAMBDA (int k, int j, int i) {
      // CT_EMF_ArithmeticAverage (emf, 0.25);
      const real w = ONE_FOURTH_F;
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
  idfx::popRegion();
}

template<typename Phys>
void ConstrainedTransport<Phys>::CalcContactAverage() {
  idfx::pushRegion("ConstrainedTransport::CalcContactAverage");
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

  // cell-centered EMFs
  IdefixArray3D<real> Ex1 = this->Ex1;
  IdefixArray3D<real> Ex2 = this->Ex2;
  IdefixArray3D<real> Ex3 = this->Ex3;

  // sign of contact discontinuity
  IdefixArray3D<real> wsx = this->svx;
  IdefixArray3D<real> wsy = this->svy;
  IdefixArray3D<real> wsz = this->svz;

#if MHD == YES && DIMENSIONS >= 2

  idefix_for("EMF_Integrate_to_Corner",
            data->beg[KDIR],data->end[KDIR]+KOFFSET,
            data->beg[JDIR],data->end[JDIR]+JOFFSET,
            data->beg[IDIR],data->end[IDIR]+IOFFSET,
    KOKKOS_LAMBDA (int k, int j, int i) {
      real ez_l2 = (1-wsx(k,j-1,i)) * (ezj(k,j,i)   - Ex3(k,j-1,i)) +
                   (  wsx(k,j-1,i)) * (ezj(k,j,i-1) - Ex3(k,j-1,i-1));

      real ez_r2 = (1-wsx(k,j,i)) * (ezj(k,j,i)   - Ex3(k,j,i)) +
                   (  wsx(k,j,i)) * (ezj(k,j,i-1) - Ex3(k,j,i-1));

      real ez_l1 = (1-wsy(k,j,i-1)) * (ezi(k,j,i)   - Ex3(k,j,i-1)) +
                   (  wsy(k,j,i-1)) * (ezi(k,j-1,i) - Ex3(k,j-1,i-1));

      real ez_r1 = (1-wsy(k,j,i)) * (ezi(k,j,i)   - Ex3(k,j,i)) +
                   (  wsy(k,j,i)) * (ezi(k,j-1,i) - Ex3(k,j-1,i));

      ez(k,j,i) = ONE_FOURTH_F * (ez_l2 + ez_r2 + ez_l1 + ez_r1 +
                        #if DIMENSIONS >= 2
                          ezi(k,j,i) + ezi(k,j-1,i) + ezj(k,j,i) + ezj(k,j,i-1)
                        #else
                          TWO_F*ezi(k,j,i) + ezj(k,j,i) + ezj(k,j,i-1)
                        #endif
                        );

      #if DIMENSIONS == 3
        real ex_l3 = (1-wsy(k-1,j,i)) * (exk(k,j,i)   - Ex1(k-1,j,i)) +
                     (  wsy(k-1,j,i)) * (exk(k,j-1,i) - Ex1(k-1,j-1,i));

        real ex_r3 = (1-wsy(k,j,i)) * (exk(k,j,i)   - Ex1(k,j,i)) +
                     (  wsy(k,j,i)) * (exk(k,j-1,i) - Ex1(k,j-1,i));

        real ex_l2 = (1-wsz(k,j-1,i)) * (exj(k,j,i)   - Ex1(k,j-1,i)) +
                     (  wsz(k,j-1,i)) * (exj(k-1,j,i) - Ex1(k-1,j-1,i));

        real ex_r2 = (1-wsz(k,j,i)) * (exj(k,j,i)   - Ex1(k,j,i)) +
                     (  wsz(k,j,i)) * (exj(k-1,j,i) - Ex1(k-1,j,i));

        ex(k,j,i) = ONE_FOURTH_F * (ex_l3 + ex_r3 + ex_l2 + ex_r2 +
                            exj(k,j,i) + exj(k-1,j,i) + exk(k,j,i) + exk(k,j-1,i) );

        real ey_l3 = (1-wsx(k-1,j,i)) * (eyk(k,j,i)   - Ex2(k-1,j,i)) +
                     (  wsx(k-1,j,i)) * (eyk(k,j,i-1) - Ex2(k-1,j,i-1));

        real ey_r3 = (1-wsx(k,j,i)) * (eyk(k,j,i)   - Ex2(k,j,i)) +
                     (  wsx(k,j,i)) * (eyk(k,j,i-1) - Ex2(k,j,i-1));

        real ey_l1 = (1-wsz(k,j,i-1)) * (eyi(k,j,i)   - Ex2(k,j,i-1)) +
                     (  wsz(k,j,i-1)) * (eyi(k-1,j,i) - Ex2(k-1,j,i-1));

        real ey_r1 = (1-wsz(k,j,i)) * (eyi(k,j,i)   - Ex2(k,j,i)) +
                     (  wsz(k,j,i)) * (eyi(k-1,j,i) - Ex2(k-1,j,i));

        ey(k,j,i) = ONE_FOURTH_F * (ey_l3 + ey_r3 + ey_l1 + ey_r1 +
                            eyi(k,j,i) + eyi(k-1,j,i) + eyk(k,j,i) + eyk(k,j,i-1));


      #endif
    });

  /*
  idefix_for("EMF_Integrate to Corner",
            data->beg[KDIR]-KOFFSET,data->end[KDIR]+KOFFSET,
            data->beg[JDIR]-JOFFSET,data->end[JDIR]+JOFFSET,
            data->beg[IDIR]-IOFFSET,data->end[IDIR]+IOFFSET,
    KOKKOS_LAMBDA (int k, int j, int i) {
      [[maybe_unused]] int iu, ju, ku;
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
  */

#endif // MHD
  idfx::popRegion();
}
#endif // FLUID_CONSTRAINEDTRANSPORT_CALCCORNEREMF_HPP_
