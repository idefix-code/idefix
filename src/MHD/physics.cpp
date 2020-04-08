#include "../idefix.hpp"
#include "physics.hpp"

/********************************
 * Local Kokkos Inline function *
 * ******************************/

KOKKOS_INLINE_FUNCTION void K_Flux(real F[], real V[], real U[], real C2Iso, 
                                    const int VXn, const int VXt, const int VXb,
                                    const int BXn, const int BXt, const int BXb,
                                    const int MXn) {

    F[RHO] = U[MXn];
    EXPAND( F[MX1] = U[MX1]*V[VXn] - V[BXn]*V[BX1]; ,
            F[MX2] = U[MX2]*V[VXn] - V[BXn]*V[BX2]; ,
            F[MX3] = U[MX3]*V[VXn] - V[BXn]*V[BX3];)

    EXPAND(F[BXn] = ZERO_F;                             ,
           F[BXt] = V[VXn]*V[BXt] - V[BXn]*V[VXt];   ,
           F[BXb] = V[VXn]*V[BXb] - V[BXn]*V[VXb]; )

    real Bmag2 = EXPAND(V[BX1]*V[BX1] , + V[BX2]*V[BX2], + V[BX3]*V[BX3]);

    #if HAVE_ENERGY
    real ptot  = V[PRS] + HALF_F*Bmag2;
    #elif EOS == ISOTHERMAL
    real ptot  = C2Iso * V[RHO] + HALF_F*Bmag2;
    #else
    #error "K_Flux not defined for this EOS!"
    #endif

    #if HAVE_ENERGY
    F[ENG]   = (U[ENG] + ptot)*V[VXn] - V[BXn] * (EXPAND(V[VX1]*V[BX1] , + V[VX2]*V[BX2], + V[VX3]*V[BX3]));
    #endif

    // Add back pressure in the flux (not included in original PLUTO implementation)
    F[MXn]   += ptot;
} 

KOKKOS_INLINE_FUNCTION void K_ConsToPrim(real Vc[], real Uc[], real gamma_m1) {
    

    Vc[RHO] = Uc[RHO];

    EXPAND( Vc[VX1] = Uc[MX1]/Uc[RHO]; ,
            Vc[VX2] = Uc[MX2]/Uc[RHO]; ,
            Vc[VX3] = Uc[MX3]/Uc[RHO];)

    EXPAND( Vc[BX1] = Uc[BX1];  ,
            Vc[BX2] = Uc[BX2];  ,
            Vc[BX3] = Uc[BX3];  )
         

#if HAVE_ENERGY
    real kin,mag;
    kin = HALF_F / Uc[RHO] * (EXPAND(    Uc[MX1]*Uc[MX1] , 
                                    + Uc[MX2]*Uc[MX2] ,
                                    + Uc[MX3]*Uc[MX3] ));

    mag = HALF_F * (EXPAND(    Uc[BX1]*Uc[BX1] , 
                                    + Uc[BX2]*Uc[BX2] ,
                                    + Uc[BX3]*Uc[BX3]));


    Vc[PRS] = gamma_m1 * (Uc[ENG] - kin - mag);
    
    // Check pressure positivity
    if(Vc[PRS]<= ZERO_F) {
        Vc[PRS] = SMALL_PRESSURE_FIX;
        Uc[ENG] = Vc[PRS]/gamma_m1+kin+mag;
    }
#endif  // Have_energy
}

KOKKOS_INLINE_FUNCTION void K_PrimToCons(real Uc[], real Vc[], real gamma_m1) {

    Uc[RHO] = Vc[RHO];

    EXPAND( Uc[MX1] = Vc[VX1]*Vc[RHO]; ,
            Uc[MX2] = Vc[VX2]*Vc[RHO]; ,
            Uc[MX3] = Vc[VX3]*Vc[RHO];)


    EXPAND( Uc[BX1] = Vc[BX1];  ,
            Uc[BX2] = Vc[BX2];  ,
            Uc[BX3] = Vc[BX3];  )

#if HAVE_ENERGY

    Uc[ENG] = Vc[PRS] / gamma_m1 
                + HALF_F * Vc[RHO] * (EXPAND(  Vc[VX1]*Vc[VX1] , 
                                         + Vc[VX2]*Vc[VX2] ,
                                         + Vc[VX3]*Vc[VX3] ))
                + HALF_F * (EXPAND(    Uc[BX1]*Uc[BX1] , 
                                     + Uc[BX2]*Uc[BX2] ,
                                     + Uc[BX3]*Uc[BX3])); 
#endif  // Have_energy

}




Physics::Physics(Input &input, Setup &setup) {
    Kokkos::Profiling::pushRegion("Physics::Physics(DataBock)");

    this->gamma = 5.0/3.0;
    this->C2Iso = 1.0;

    this->mySetup=setup;

    Kokkos::Profiling::popRegion();
}

Physics::Physics() {

}


// Convect Conservative to Primitive variable
void Physics::ConvertConsToPrim(DataBlock & data) {

    Kokkos::Profiling::pushRegion("Physics::ConvertConsToPrim");

    IdefixArray4D<real> Vc = data.Vc;
    IdefixArray4D<real> Uc = data.Uc;
    real gamma_m1=this->gamma-ONE_F;

    idefix_for("ConsToPrim",0,data.np_tot[KDIR],0,data.np_tot[JDIR],0,data.np_tot[IDIR],
                        KOKKOS_LAMBDA (int k, int j, int i) 
            {
                real U[NVAR];
                real V[NVAR];
                for(int nv = 0 ; nv < NVAR; nv++) {
                    U[nv] = Uc(nv,k,j,i); 
                }

                K_ConsToPrim(V,U,gamma_m1);

                for(int nv = 0 ; nv<NVAR; nv++) {
                    Vc(nv,k,j,i) = V[nv];
                }
            });

    Kokkos::Profiling::popRegion();

}


// Convert Primitive to conservative variables
void Physics::ConvertPrimToCons(DataBlock & data) {

    Kokkos::Profiling::pushRegion("Physics::ConvertPrimToCons");

    IdefixArray4D<real> Vc = data.Vc;
    IdefixArray4D<real> Uc = data.Uc;
    real gamma_m1=this->gamma-ONE_F;

    idefix_for("ConvertPrimToCons",0,data.np_tot[KDIR],0,data.np_tot[JDIR],0,data.np_tot[IDIR],
                        KOKKOS_LAMBDA (int k, int j, int i) 
            {
                real U[NVAR];
                real V[NVAR];
                for(int nv = 0 ; nv < NVAR; nv++) {
                    V[nv] = Vc(nv,k,j,i); 
                }

                K_PrimToCons(U,V,gamma_m1);

                for(int nv = 0 ; nv<NVAR; nv++) {
                    Uc(nv,k,j,i) = U[nv];
                }
            });

    Kokkos::Profiling::popRegion();
}


// Build a left and right extrapolation of the primitive variables along direction dir

// These functions extrapolate the cell prim vars to the faces. Definitions are as followed
//
// |       cell i-1               interface i          cell i
// |-----------------------------------|------------------------------------||
//          Vc(i-1)           PrimL(i)  PrimR(i)       Vc(i)

void Physics::ExtrapolatePrimVar(DataBlock &data, int dir) {
    int ioffset,joffset,koffset;

    Kokkos::Profiling::pushRegion("Physics::ExtrapolatePrimVar");
    ioffset=joffset=koffset=0;
    // Determine the offset along which we do the extrapolation
    if(dir==IDIR) ioffset=1;
    if(dir==JDIR) joffset=1;
    if(dir==KDIR) koffset=1;

    IdefixArray4D<real> Vc = data.Vc;
    IdefixArray4D<real> PrimL = data.PrimL;
    IdefixArray4D<real> PrimR = data.PrimR;


#if ORDER == 1

    idefix_for("ExtrapolatePrimVar",0,NVAR,data.beg[KDIR],data.end[KDIR]+koffset,data.beg[JDIR],data.end[JDIR]+joffset,data.beg[IDIR],data.end[IDIR]+ioffset,
                        KOKKOS_LAMBDA (int n, int k, int j, int i) 
            {   
                
                PrimL(n,k,j,i) = Vc(n,k-koffset,j-joffset,i-ioffset);
                PrimR(n,k,j,i) = Vc(n,k,j,i);
            });

#elif ORDER == 2
    idefix_for("ExtrapolatePrimVar",0,NVAR,data.beg[KDIR]-koffset,data.end[KDIR]+koffset,data.beg[JDIR]-joffset,data.end[JDIR]+joffset,data.beg[IDIR]-ioffset,data.end[IDIR]+ioffset,
                        KOKKOS_LAMBDA (int n, int k, int j, int i) 
            {
                real dvm = Vc(n,k,j,i)-Vc(n,k-koffset,j-joffset,i-ioffset);
                real dvp = Vc(n,k+koffset,j+joffset,i+ioffset) - Vc(n,k,j,i);

                // Van Leer limiter
                real dv = (dvp*dvm > ZERO_F ? TWO_F*dvp*dvm/(dvp + dvm) : ZERO_F);

                PrimL(n,k+koffset,j+joffset,i+ioffset) = Vc(n,k,j,i) + HALF_F*dv;
                PrimR(n,k,j,i) = Vc(n,k,j,i) - HALF_F*dv;

            });
#else   
        #error ORDER should be 1 or 2
#endif



    Kokkos::Profiling::popRegion();
}

// Compute Riemann fluxes from states
void Physics::CalcRiemannFlux(DataBlock & data, int dir) {
    int ioffset,joffset,koffset;

    Kokkos::Profiling::pushRegion("Physics::CalcRiemannFlux");
    ioffset=joffset=koffset=0;
    // Determine the offset along which we do the extrapolation
    if(dir==IDIR) ioffset=1;
    if(dir==JDIR) joffset=1;
    if(dir==KDIR) koffset=1;

    IdefixArray4D<real> PrimL = data.PrimL;
    IdefixArray4D<real> PrimR = data.PrimR;
    IdefixArray4D<real> Flux = data.FluxRiemann;
    IdefixArray1D<real> dx = data.dx[dir];
    IdefixArray3D<real> invDt = data.InvDtHyp;

    // References to required emf components
    IdefixArray3D<real> Eb;
    IdefixArray3D<real> Et;


    real gamma_m1=this->gamma-ONE_F;
    real gamma=this->gamma;
    real C2Iso = this->C2Iso;

    // Define normal, tangent and bi-tanget indices

    int nDIR, tDIR, bDIR;
    int VXn, VXt, VXb;
    int BXn, BXt, BXb;
    int MXn, MXt, MXb;

    real st, sb;      // st and sb will be useful only when Hall is included
    switch(dir) {
        case(IDIR):
            EXPAND(VXn = MXn = VX1; 
                   BXn = BX1;        , 
                   VXt = MXt = VX2; 
                   BXt = BX2;        , 
                   VXb = MXb = VX3;
                   BXb = BX3;       )

            Et = data.emf.ezi;
            Eb = data.emf.eyi;

            st = -1.0;
            sb = +1.0;
            break;
        case(JDIR):
            EXPAND(VXn = MXn = VX2; 
                   BXn = BX2;        , 
                   VXt = MXt = VX1; 
                   BXt = BX1;        , 
                   VXb = MXb = VX3;
                   BXb = BX3;       )

            Et = data.emf.ezj;
            Eb = data.emf.exj;

            st = +1.0;
            sb = -1.0;
            break;
        case(KDIR):
            EXPAND(VXn = MXn = VX3; 
                   BXn = BX3;        , 
                   VXt = MXt = VX1; 
                   BXt = BX1;        , 
                   VXb = MXb = VX2;
                   BXb = BX2;       )

            Et = data.emf.eyk;
            Eb = data.emf.exk;

            st = -1.0;
            sb = +1.0;
            break;
        default:
            IDEFIX_ERROR("Wrong direction");
    }

    nDIR = VXn-VX1; tDIR = VXt-VX1; bDIR = VXb-VX1;


    idefix_for("CalcRiemannFlux",data.beg[KDIR],data.end[KDIR]+koffset,data.beg[JDIR],data.end[JDIR]+joffset,data.beg[IDIR],data.end[IDIR]+ioffset,
                        KOKKOS_LAMBDA (int k, int j, int i) 
            {

                // Primitive variables
                real vL[NVAR];
                real vR[NVAR];
                real vRL[NVAR];

                // Conservative variables
                real uL[NVAR];
                real uR[NVAR];

                // Flux (left and right)
                real fluxL[NVAR];
                real fluxR[NVAR];

                // Signal speeds
                real cRL, cmax;
                real gpr, Bt2, B2;

                // 1-- Store the primitive variables on the left, right, and averaged states
                for(int nv = 0 ; nv < NVAR; nv++) {
                    vL[nv] = PrimL(nv,k,j,i);
                    vR[nv] = PrimR(nv,k,j,i);
                    vRL[nv]=HALF_F*(vL[nv]+vR[nv]);
                }

                // 2-- Compute the conservative variables
                K_PrimToCons(uL, vL, gamma_m1);
                K_PrimToCons(uR, vR, gamma_m1);

                // 3-- Compute the left and right fluxes
                K_Flux(fluxL, vL, uL, C2Iso, VXn, VXt, VXb, BXn, BXt, BXb, MXn);
                K_Flux(fluxR, vR, uR, C2Iso, VXn, VXt, VXb, BXn, BXt, BXb, MXn);


                // 4-- Get the wave speed
                
            #if HAVE_ENERGY
                gpr=(gamma_m1+ONE_F)*vRL[PRS];
            #else
                gpr=C2Iso*vRL[RHO];
            #endif
                Bt2=EXPAND(ZERO_F    ,
                            + vRL[BXt]*vRL[BXt],
                            + vRL[BXb]*vRL[BXb]);

                B2=Bt2 + vRL[BXn]*vRL[BXn];

                cRL = gpr - B2;
                cRL = cRL + B2 + SQRT(cRL*cRL + FOUR_F*gpr*Bt2);
                cRL = SQRT(HALF_F * cRL/vRL[RHO]);

                cmax = FMAX(FABS(vRL[VXn]+cRL),FABS(vRL[VXn]-cRL));

                // 5-- Compute the flux from the left and right states
                for(int nv = 0 ; nv < NVAR; nv++) {
                    Flux(nv,k,j,i) = HALF_F*(fluxL[nv]+fluxR[nv] - cmax*(uR[nv]-uL[nv]));
                }

                //6-- Compute maximum dt for this sweep
                const int ig = ioffset*i + joffset*j + koffset*k;

                invDt(k,j,i) += cmax/dx(ig);

                // 7-- Store the flux in the emf components
                D_EXPAND(Et(k,j,i) = st*Flux(BXt,k,j,i); ,
                                                         ,
                         Eb(k,j,i) = st*Flux(BXb,k,j,i); )

            });

    Kokkos::Profiling::popRegion();

}

// Compute the right handside in direction dir from conservative equation, with timestep dt
void Physics::CalcRightHandSide(DataBlock &data, int dir, real dt) {

    Kokkos::Profiling::pushRegion("Physics::CalcRightHandSide");
    IdefixArray4D<real> Uc = data.Uc;
    IdefixArray1D<real> dx = data.dx[dir];
    IdefixArray4D<real> Flux = data.FluxRiemann;

    int ioffset,joffset,koffset;
    ioffset=joffset=koffset=0;
    // Determine the offset along which we do the extrapolation
    if(dir==IDIR) ioffset=1;
    if(dir==JDIR) joffset=1;
    if(dir==KDIR) koffset=1;

    idefix_for("CalcRightHandSide",0,NVAR,data.beg[KDIR],data.end[KDIR],data.beg[JDIR],data.end[JDIR],data.beg[IDIR],data.end[IDIR],
        KOKKOS_LAMBDA (int n, int k, int j, int i) {
            
            const int ig = ioffset*i + joffset*j + koffset*k;

            Uc(n,k,j,i) = Uc(n,k,j,i) - dt / dx(ig) * (Flux(n, k+koffset, j+joffset, i+ioffset) - Flux(n, k, j, i));

        });

    Kokkos::Profiling::popRegion();
}


// Set Boundary conditions
void Physics::SetBoundary(DataBlock &data, real t) {

    Kokkos::Profiling::pushRegion("Physics::SetBoundary");

    IdefixArray4D<real> Vc = data.Vc;
    int ibeg,iend,jbeg,jend,kbeg,kend;
    int ioffset,joffset,koffset;
    int ighost,jghost,kghost;

    ighost = data.nghost[IDIR];
    jghost = data.nghost[JDIR];
    kghost = data.nghost[KDIR];

    // X1 boundary conditions
    

    for(int dir=0 ; dir < DIMENSIONS ; dir++ ) {

        ioffset = (dir == IDIR) ? data.np_int[IDIR] : 0;
        joffset = (dir == JDIR) ? data.np_int[JDIR] : 0;
        koffset = (dir == KDIR) ? data.np_int[KDIR] : 0;


        // left boundary
        ibeg=0;
        iend= (dir == IDIR) ? ighost : data.np_tot[IDIR];
        jbeg=0;
        jend= (dir == JDIR) ? jghost : data.np_tot[JDIR];
        kbeg=0;
        kend= (dir == KDIR) ? kghost : data.np_tot[KDIR];

        switch(data.lbound[dir]) {
            case periodic:
                idefix_for("BoundaryBegPeriodic",0,NVAR,kbeg,kend,jbeg,jend,ibeg,iend,
                    KOKKOS_LAMBDA (int n, int k, int j, int i) {
                
                        Vc(n,k,j,i) = Vc(n,k+koffset,j+joffset,i+ioffset);
                    });
                break;
            case outflow:
                idefix_for("BoundaryBegOutflow",0,NVAR,kbeg,kend,jbeg,jend,ibeg,iend,
                    KOKKOS_LAMBDA (int n, int k, int j, int i) {
                        int iref= (dir==IDIR) ? ighost : i;
                        int jref= (dir==JDIR) ? jghost : j;
                        int kref= (dir==KDIR) ? kghost : k;

                        Vc(n,k,j,i) = Vc(n,kref,jref,iref);
                    });
                break;
            default:
                std::stringstream msg ("Boundary condition type is not yet implemented");
                IDEFIX_ERROR(msg);
        }

        // right boundary
        ibeg= (dir == IDIR) ? ioffset + ighost : 0;
        iend = data.np_tot[IDIR];
        jbeg= (dir == JDIR) ? joffset + jghost : 0;
        jend = data.np_tot[JDIR];
        kbeg= (dir == KDIR) ? koffset + kghost : 0;
        kend = data.np_tot[KDIR];

        switch(data.rbound[dir]) {
            case periodic:
                idefix_for("BoundaryEndPeriodic",0,NVAR,kbeg,kend,jbeg,jend,ibeg,iend,
                    KOKKOS_LAMBDA (int n, int k, int j, int i) {
                
                        Vc(n,k,j,i) = Vc(n,k-koffset,j-joffset,i-ioffset);
                    });
                break;
            case outflow:
                idefix_for("BoundaryEndOutflow",0,NVAR,kbeg,kend,jbeg,jend,ibeg,iend,
                    KOKKOS_LAMBDA (int n, int k, int j, int i) {
                        int iref= (dir==IDIR) ? ighost + ioffset - 1 : i;
                        int jref= (dir==JDIR) ? jghost + joffset - 1 : j;
                        int kref= (dir==KDIR) ? kghost + koffset - 1 : k;

                        Vc(n,k,j,i) = Vc(n,kref,jref,iref);
                    });
                break;
            default:
                std::stringstream msg("Boundary condition type is not yet implemented");
                IDEFIX_ERROR(msg);

        }


    }

    Kokkos::Profiling::popRegion();

}

real Physics::CheckDivB(DataBlock &data) {

    real divB;
    IdefixArray4D<real> Vc = data.Vc;
    IdefixArray1D<real> dx1 = data.dx[IDIR];
    IdefixArray1D<real> dx2 = data.dx[JDIR];
    IdefixArray1D<real> dx3 = data.dx[KDIR];

    Kokkos::parallel_reduce("CheckDivB",
                                Kokkos::MDRangePolicy<Kokkos::Rank<3, Kokkos::Iterate::Right, Kokkos::Iterate::Right>>
                                ({data.beg[KDIR]+1,data.beg[JDIR]+1,data.beg[IDIR]+1},{data.end[KDIR]+1, data.end[JDIR]+1, data.end[IDIR]+1}),
                                KOKKOS_LAMBDA (int k, int j, int i, real &divBmax) {
                real dB1,dB2,dB3;

                dB1=dB2=dB3=ZERO_F;

                D_EXPAND( dB1=(Vs(BX1s,k,j,i+1)-Vs(BX1s,k,j,i))/(dx1(i)); ,
                          dB2=(Vs(BX2s,k,j+1,i)-Vs(BX2s,k,j,i))/(dx2(j)); ,
                          dB3=(Vs(BX3s,k+1,j,i)-Vs(BX3s,k,j,i))/(dx3(k));  )
                

                divBmax=FMAX(FABS(D_EXPAND(dB1, +dB2, +dB3)),divBmax);

            }, Kokkos::Max<real>(divB) );

    return(divB);

}

void Physics::SetGamma(real newGamma) {
    this->gamma=newGamma;
}