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
    int iextend, jextend,kextend;
    int BXn;

    Kokkos::Profiling::pushRegion("Physics::ExtrapolatePrimVar");
    // Offset is in the direction of integration
    ioffset=joffset=koffset=0;

    // extension if perp to the direction of integration, as required by CT.
    iextend=jextend=kextend=0;

    // Determine the offset along which we do the extrapolation, as well as the perp extension
    if(dir==IDIR) {
        ioffset=1;
        BXn = BX1; 
        D_EXPAND(           ,
                   jextend = 1; ,
                   kextend = 1; )
    }
    if(dir==JDIR) { 
        joffset=1;
        BXn = BX2; 
        D_EXPAND( iextend = 1;   ,
                                 ,
                  kextend = 1;)
    }
    if(dir==KDIR) {
        koffset=1;
        BXn = BX3; 
        D_EXPAND( iextend = 1;  ,
                  jextend = 1;  ,
                   )
    }

    IdefixArray4D<real> Vc = data.Vc;
    IdefixArray4D<real> Vs = data.Vs;
    IdefixArray4D<real> PrimL = data.PrimL;
    IdefixArray4D<real> PrimR = data.PrimR;


#if ORDER == 1

    idefix_for("ExtrapolatePrimVar",0,NVAR,data.beg[KDIR]-kextend,data.end[KDIR]+koffset+kextend,data.beg[JDIR]-jextend,data.end[JDIR]+joffset+jextend,data.beg[IDIR]-iextend,data.end[IDIR]+ioffset+iextend,
                        KOKKOS_LAMBDA (int n, int k, int j, int i) 
            {   
                // If normal component, the use Staggered field 
                if(n==BXn) {
                    PrimL(n,k,j,i) = Vs(dir,k,j,i);
                    PrimR(n,k,j,i) = Vs(dir,k,j,i);
                }
                else {
                    PrimL(n,k,j,i) = Vc(n,k-koffset,j-joffset,i-ioffset);
                    PrimR(n,k,j,i) = Vc(n,k,j,i);
                }


            });

#elif ORDER == 2
    idefix_for("ExtrapolatePrimVar",0,NVAR,data.beg[KDIR]-koffset-kextend,data.end[KDIR]+koffset+kextend,data.beg[JDIR]-joffset-jextend,data.end[JDIR]+joffset+jextend,data.beg[IDIR]-ioffset-iextend,data.end[IDIR]+ioffset+iextend,
                        KOKKOS_LAMBDA (int n, int k, int j, int i) 
            {
                if(n==BXn) {
                    PrimL(n,k+koffset,j+joffset,i+ioffset) = Vs(dir,k+koffset,j+joffset,i+ioffset);
                    PrimR(n,k,j,i) = Vs(dir,k,j,i);
                }
                else {
                    real dvm = Vc(n,k,j,i)-Vc(n,k-koffset,j-joffset,i-ioffset);
                    real dvp = Vc(n,k+koffset,j+joffset,i+ioffset) - Vc(n,k,j,i);

                    // Van Leer limiter
                    real dv = (dvp*dvm > ZERO_F ? TWO_F*dvp*dvm/(dvp + dvm) : ZERO_F);

                    PrimL(n,k+koffset,j+joffset,i+ioffset) = Vc(n,k,j,i) + HALF_F*dv;
                    PrimR(n,k,j,i) = Vc(n,k,j,i) - HALF_F*dv;
                }

            });
#else   
        #error ORDER should be 1 or 2
#endif



    Kokkos::Profiling::popRegion();
}

// Compute Riemann fluxes from states
void Physics::CalcRiemannFlux(DataBlock & data, int dir) {
    int ioffset,joffset,koffset;
    int iextend, jextend,kextend;

    Kokkos::Profiling::pushRegion("Physics::CalcRiemannFlux");
    
    ioffset=joffset=koffset=0;
    // extension in perp to the direction of integration, as required by CT.
    iextend=jextend=kextend=0;

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
            ioffset = 1;
            D_EXPAND(           ,
                   jextend = 1; ,
                   kextend = 1; )

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
            joffset=1;
            D_EXPAND( iextend = 1;   ,
                                    ,
                    kextend = 1;)
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
            koffset=1;
            D_EXPAND( iextend = 1;               ,
                    jextend = 1;                ,
                    )
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


    idefix_for("CalcRiemannFlux",data.beg[KDIR]-kextend,data.end[KDIR]+koffset+kextend,data.beg[JDIR]-jextend,data.end[JDIR]+joffset+jextend,data.beg[IDIR]-iextend,data.end[IDIR]+ioffset+iextend,
                        KOKKOS_LAMBDA (int k, int j, int i) 
            {
                // Primitive variables 
                real v[NVAR];
                real u[NVAR];
                real flux[NVAR];
                real fluxRiemann[NVAR];

                // Store the average primitive variables
                for(int nv = 0 ; nv < NVAR; nv++) {
                    v[nv] = HALF_F*(PrimL(nv,k,j,i) + PrimR(nv,k,j,i));
                }


                // 4-- Get the wave speed
                // Signal speeds
                real cRL, cmax;
                real gpr, Bt2, B2;

                #if HAVE_ENERGY
                    gpr=(gamma_m1+ONE_F)*v[PRS];
                #else
                    gpr=C2Iso*v[RHO];
                #endif
                Bt2=EXPAND(ZERO_F    ,
                            + v[BXt]*v[BXt],
                            + v[BXb]*v[BXb]);

                B2=Bt2 + v[BXn]*v[BXn];

                cRL = gpr - B2;
                cRL = cRL + B2 + SQRT(cRL*cRL + FOUR_F*gpr*Bt2);
                cRL = SQRT(HALF_F * cRL/v[RHO]);

                cmax = FMAX(FABS(v[VXn]+cRL),FABS(v[VXn]-cRL));


                // Load the left state
                for(int nv = 0 ; nv < NVAR; nv++) {
                    v[nv] = PrimL(nv,k,j,i);
                }

                // 2-- Compute the conservative variables
                K_PrimToCons(u, v, gamma_m1);

                // 3-- Compute the left and right fluxes
                K_Flux(flux, v, u, C2Iso, VXn, VXt, VXb, BXn, BXt, BXb, MXn);
                

                // 5-- Compute the flux from the left and right states
                for(int nv = 0 ; nv < NVAR; nv++) {
                    fluxRiemann[nv] = flux[nv] + cmax*u[nv];
                }

                // Load the right state
                for(int nv = 0 ; nv < NVAR; nv++) {
                    v[nv] = PrimR(nv,k,j,i);
                }

                // 2-- Compute the conservative variables
                K_PrimToCons(u, v, gamma_m1);

                // 3-- Compute the left and right fluxes
                K_Flux(flux, v, u, C2Iso, VXn, VXt, VXb, BXn, BXt, BXb, MXn);
                
                // 5-- Compute the flux from the left and right states
                for(int nv = 0 ; nv < NVAR; nv++) {
                    Flux(nv,k,j,i) = HALF_F*(fluxRiemann[nv]+flux[nv] - cmax*u[nv]);
                }

                //6-- Compute maximum dt for this sweep
                const int ig = ioffset*i + joffset*j + koffset*k;

                invDt(k,j,i) = FMAX(cmax/dx(ig),invDt(k,j,i));

                // 7-- Store the flux in the emf components
                D_EXPAND(Et(k,j,i) = st*Flux(BXt,k,j,i); ,
                                                         ,
                         Eb(k,j,i) = sb*Flux(BXb,k,j,i); )

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

    idefix_for("CalcRightHandSide",data.beg[KDIR],data.end[KDIR],data.beg[JDIR],data.end[JDIR],data.beg[IDIR],data.end[IDIR],
        KOKKOS_LAMBDA (int k, int j, int i) {
            
            const int ig = ioffset*i + joffset*j + koffset*k;
            real dtdx=dt / dx(ig);

            for(int nv = 0 ; nv < NVAR ; nv++) {
                // Do not evolve the field components if they are computed by CT (i.e. if they are in Vs)
                
                D_EXPAND( if(nv == BX1) continue;   ,
                          if(nv == BX2) continue;   ,
                          if(nv == BX3) continue;  ) 
                

                Uc(nv,k,j,i) = Uc(nv,k,j,i) -  dtdx*(Flux(nv, k+koffset, j+joffset, i+ioffset) - Flux(nv, k, j, i));
            }

            

        });

    Kokkos::Profiling::popRegion();
}

// Compute Corner EMFs from the one stored in the Riemann step
void Physics::CalcCornerEMF(DataBlock &data, real t) {
        Kokkos::Profiling::pushRegion("Physics::CalcCornerEMF");

        // Corned EMFs
        IdefixArray3D<real> Ex = data.emf.ex;
        IdefixArray3D<real> Ey = data.emf.ey;
        IdefixArray3D<real> Ez = data.emf.ez;

        // Face-centered EMFs
        IdefixArray3D<real> exj = data.emf.exj;
        IdefixArray3D<real> exk = data.emf.exk;
        IdefixArray3D<real> eyi = data.emf.eyi;
        IdefixArray3D<real> eyk = data.emf.eyk;
        IdefixArray3D<real> ezi = data.emf.ezi;
        IdefixArray3D<real> ezj = data.emf.ezj;

        real w = ONE_FOURTH_F;


        idefix_for("CalcCornerEMF",data.beg[KDIR],data.end[KDIR]+KOFFSET,data.beg[JDIR],data.end[JDIR]+JOFFSET,data.beg[IDIR],data.end[IDIR]+IOFFSET,
                    KOKKOS_LAMBDA (int k, int j, int i) {
                        #if DIMENSIONS == 3
                        Ex(k,j,i) = w * (exj(k,j,i) + exj(k-1,j,i) + exk(k,j,i) + exk(k,j-1,i));
                        Ey(k,j,i) = w * (eyi(k,j,i) + eyi(k-1,j,i) + eyk(k,j,i) + eyk(k,j,i-1));
                        #endif
                        #if DIMENSIONS >= 2
                        Ez(k,j,i) = w * (ezi(k,j,i) + ezi(k,j-1,i) + ezj(k,j,i) + ezj(k,j,i-1));
                        #else
                        Ez(k,j,i) = w * (TWO_F*ezi(k,j,i) + ezj(k,j,i) + ezj(k,j,i-1));
                        #endif

                    });

        Kokkos::Profiling::popRegion();
}

// Evolve the magnetic field in Vs according to Constranied transport
void Physics::EvolveMagField(DataBlock &data, real t, real dt) {
    Kokkos::Profiling::pushRegion("Physics::EvolveMagField");

    // Corned EMFs
    IdefixArray3D<real> Ex1 = data.emf.ex;
    IdefixArray3D<real> Ex2 = data.emf.ey;
    IdefixArray3D<real> Ex3 = data.emf.ez;

    // Field
    IdefixArray4D<real> Vs = data.Vs;

    // Coordinates
    IdefixArray1D<real> x1=data.x[IDIR];
    IdefixArray1D<real> x2=data.x[JDIR];
    IdefixArray1D<real> x3=data.x[KDIR];

    IdefixArray1D<real> dx1=data.dx[IDIR];
    IdefixArray1D<real> dx2=data.dx[JDIR];
    IdefixArray1D<real> dx3=data.dx[KDIR];



    idefix_for("EvolvMagField",data.beg[KDIR],data.end[KDIR]+KOFFSET,data.beg[JDIR],data.end[JDIR]+JOFFSET,data.beg[IDIR],data.end[IDIR]+IOFFSET,
                    KOKKOS_LAMBDA (int k, int j, int i) {

                        Vs(BX1s,k,j,i) = Vs(BX1s,k,j,i) + ( D_EXPAND( ZERO_F                                     ,
                                                                      - dt/dx2(j) * (Ex3(k,j+1,i) - Ex3(k,j,i) )  ,
                                                                      + dt/dx3(k) * (Ex2(k+1,j,i) - Ex2(k,j,i) ) ));
                        #if DIMENSIONS >= 2
                        Vs(BX2s,k,j,i) = Vs(BX2s,k,j,i) + ( D_EXPAND(   dt/dx1(i) * (Ex3(k,j,i+1) - Ex3(k,j,i) )  ,
                                                                                                                  ,
                                                                      - dt/dx3(k) * (Ex1(k+1,j,i) - Ex1(k,j,i) ) ));
                        #endif
                        #if DIMENSIONS == 3
                        Vs(BX3s,k,j,i) = Vs(BX3s,k,j,i) + (  - dt/dx1(i) * (Ex2(k,j,i+1) - Ex2(k,j,i) )  
                                                             + dt/dx2(j) * (Ex1(k,j+1,i) - Ex1(k,j,i) ) );

                        #endif

                    });
    Kokkos::Profiling::popRegion();
}

void Physics::ReconstructVcField(DataBlock & data,  IdefixArray4D<real> &Vc) {
    Kokkos::Profiling::pushRegion("Physics::ReconstructVcField");
    IdefixArray4D<real> Vs=data.Vs;

    // Reconstruct cell average field when using CT
    idefix_for("ReconstructVcMagField",0,data.np_tot[KDIR],0,data.np_tot[JDIR],0,data.np_tot[IDIR],
                    KOKKOS_LAMBDA (int k, int j, int i) {
                        D_EXPAND(   Vc(BX1,k,j,i) = HALF_F * (Vs(BX1s,k,j,i) + Vs(BX1s,k,j,i+1)) ;  ,
                                    Vc(BX2,k,j,i) = HALF_F * (Vs(BX2s,k,j,i) + Vs(BX2s,k,j+1,i)) ;  ,
                                    Vc(BX3,k,j,i) = HALF_F * (Vs(BX3s,k,j,i) + Vs(BX3s,k+1,j,i)) ; )
                       
                    });
    Kokkos::Profiling::popRegion();
}



void Physics::ReconstructNormalField(DataBlock &data) {
    Kokkos::Profiling::pushRegion("Physics::ReconstructNormalField");

    // Reconstruct the field
    IdefixArray4D<real> Vc = data.Vc;
    IdefixArray4D<real> Vs = data.Vs;
    // Coordinates
    IdefixArray1D<real> x1=data.x[IDIR];
    IdefixArray1D<real> x2=data.x[JDIR];
    IdefixArray1D<real> x3=data.x[KDIR];
    IdefixArray1D<real> dx1=data.dx[IDIR];
    IdefixArray1D<real> dx2=data.dx[JDIR];
    IdefixArray1D<real> dx3=data.dx[KDIR];

    int nstart, nend;
    int nx1,nx2,nx3;

    // reconstruct BX1s
    nstart = data.nghost[IDIR]-1;
    nend = data.np_tot[IDIR] - data.nghost[IDIR]-1;

    nx1=data.np_tot[IDIR];
    nx2=data.np_tot[JDIR];
    nx3=data.np_tot[KDIR];


    idefix_for("ReconstructBX1s",0,nx3,0,nx2,
                    KOKKOS_LAMBDA (int k, int j) {
                        
                        for(int i = nstart ; i>=0 ; i-- ) {
                            Vs(BX1s,k,j,i) = Vs(BX1s,k,j,i+1) + dx1(i)*(  D_EXPAND(      ZERO_F                                       ,                    
                                                                                     +  (Vs(BX2s,k,j+1,i) - Vs(BX2s,k,j,i))/dx2(j)  , 
                                                                                     +  (Vs(BX3s,k+1,j,i) - Vs(BX3s,k,j,i))/dx3(k)) );
                        }

                        for(int i = nend ; i<nx1 ; i++ ) {
                            Vs(BX1s,k,j,i+1) = Vs(BX1s,k,j,i) -   dx1(i)*(  D_EXPAND(      ZERO_F                                     ,              
                                                                                     +  (Vs(BX2s,k,j+1,i) - Vs(BX2s,k,j,i))/dx2(j)  , 
                                                                                     +  (Vs(BX3s,k+1,j,i) - Vs(BX3s,k,j,i))/dx3(k)) );
                        }
                        

                    });

    #if DIMENSIONS >=2
    
    nstart = data.nghost[JDIR]-1;
    nend = data.np_tot[JDIR] - data.nghost[JDIR]-1;
    idefix_for("ReconstructBX2s",0,data.np_tot[KDIR],0,data.np_tot[IDIR],
                    KOKKOS_LAMBDA (int k, int i) {
                        for(int j = nstart ; j>=0 ; j-- ) {
                            Vs(BX2s,k,j,i) = Vs(BX2s,k,j+1,i) + dx2(j)*(  D_EXPAND(     (Vs(BX1s,k,j,i+1) - Vs(BX1s,k,j,i))/dx1(i)  ,                     
                                                                                                                                      , 
                                                                                     +  (Vs(BX3s,k+1,j,i) - Vs(BX3s,k,j,i))/dx3(k)) );
                        }
                        for(int j = nend ; j<nx2 ; j++ ) {
                            Vs(BX2s,k,j+1,i) = Vs(BX2s,k,j,i) -   dx2(j)*(  D_EXPAND(   (Vs(BX1s,k,j,i+1) - Vs(BX1s,k,j,i))/dx1(i)  ,                     
                                                                                                                                      , 
                                                                                     +  (Vs(BX3s,k+1,j,i) - Vs(BX3s,k,j,i))/dx3(k)) );
                        }
                        

                    });
    #endif

    #if DIMENSIONS == 3
    nstart = data.nghost[KDIR]-1;
    nend = data.np_tot[KDIR] - data.nghost[KDIR]-1;
    idefix_for("ReconstructBX3s",0,data.np_tot[JDIR],0,data.np_tot[IDIR],
                    KOKKOS_LAMBDA (int j, int i) {
                        for(int k = nstart ; k>=0 ; k-- ) {
                            Vs(BX3s,k,j,i) = Vs(BX3s,k+1,j,i) + dx3(k)*(     (Vs(BX1s,k,j,i+1) - Vs(BX1s,k,j,i+1))/dx1(i)                  
                                                                          +  (Vs(BX2s,k,j+1,i) - Vs(BX2s,k,j+1,i))/dx2(j) );
                        }
                        for(int k = nend ; k<nx3 ; k++ ) {
                            Vs(BX3s,k+1,j,i) = Vs(BX3s,k,j,i) -  dx3(k)*(     (Vs(BX1s,k,j,i+1) - Vs(BX1s,k,j,i+1))/dx1(i)                  
                                                                           +  (Vs(BX2s,k,j+1,i) - Vs(BX2s,k,j+1,i))/dx2(j) );
                        }
                        

                    });

    #endif       
    Kokkos::Profiling::popRegion();
}


// Set Boundary conditions
void Physics::SetBoundary(DataBlock &data, real t) {

    Kokkos::Profiling::pushRegion("Physics::SetBoundary");

    IdefixArray4D<real> Vc = data.Vc;
    IdefixArray4D<real> Vs = data.Vs;

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
                idefix_for("BoundaryBegPeriodicVs",0,DIMENSIONS,kbeg,kend,jbeg,jend,ibeg,iend,
                    KOKKOS_LAMBDA (int n, int k, int j, int i) {

                        // Don't touch the normal component !
                        if(n != dir) Vs(n,k,j,i) = Vs(n,k+koffset,j+joffset,i+ioffset);
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
                idefix_for("BoundaryBegOutflowVs",0,DIMENSIONS,kbeg,kend,jbeg,jend,ibeg,iend,
                    KOKKOS_LAMBDA (int n, int k, int j, int i) {
                        int iref= (dir==IDIR) ? ighost : i;
                        int jref= (dir==JDIR) ? jghost : j;
                        int kref= (dir==KDIR) ? kghost : k;

                        // Don't touch the normal component !
                        if(n != dir) Vs(n,k,j,i) = Vs(n,kref,jref,iref);
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
                idefix_for("BoundaryEndPeriodicVs",0,DIMENSIONS,kbeg,kend,jbeg,jend,ibeg,iend,
                    KOKKOS_LAMBDA (int n, int k, int j, int i) {
                        // Don't touch the normal component !
                        if(n != dir) Vs(n,k,j,i) = Vs(n,k-koffset,j-joffset,i-ioffset);                        
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
                idefix_for("BoundaryEndOutflowVs",0,DIMENSIONS,kbeg,kend,jbeg,jend,ibeg,iend,
                    KOKKOS_LAMBDA (int n, int k, int j, int i) {
                        int iref= (dir==IDIR) ? ighost + ioffset - 1 : i;
                        int jref= (dir==JDIR) ? jghost + joffset - 1 : j;
                        int kref= (dir==KDIR) ? kghost + koffset - 1 : k;

                        if(n != dir) Vs(n,k,j,i) = Vs(n,kref,jref,iref);
                    });

                break;
            default:
                std::stringstream msg("Boundary condition type is not yet implemented");
                IDEFIX_ERROR(msg);

        }
    }   // Loop on dimension ends

    // Reconstruct the normal field component when using CT
    ReconstructNormalField(data);
    
    // Remake the cell-centered field.
    ReconstructVcField(data, data.Vc);


    Kokkos::Profiling::popRegion();

}

real Physics::CheckDivB(DataBlock &data) {

    real divB;
    IdefixArray4D<real> Vs = data.Vs;
    IdefixArray1D<real> dx1 = data.dx[IDIR];
    IdefixArray1D<real> dx2 = data.dx[JDIR];
    IdefixArray1D<real> dx3 = data.dx[KDIR];

    Kokkos::parallel_reduce("CheckDivB",
                                Kokkos::MDRangePolicy<Kokkos::Rank<3, Kokkos::Iterate::Right, Kokkos::Iterate::Right>>
                                ({data.beg[KDIR],data.beg[JDIR],data.beg[IDIR]},{data.end[KDIR], data.end[JDIR], data.end[IDIR]}),
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
