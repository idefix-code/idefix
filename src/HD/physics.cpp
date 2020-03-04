#include "../idefix.hpp"
#include "physics.hpp"

/********************************
 * Local Kokkos Inline function *
 * ******************************/

KOKKOS_INLINE_FUNCTION void K_Flux(real F[], real V[], real U[], const int dir) {
    int VXn = VX1+dir;

    F[RHO] = U[VXn];

    EXPAND(F[MX1] = U[MX1]*V[VXn]; ,
           F[MX2] = U[MX2]*V[VXn]; ,
           F[MX3] = U[MX3]*V[VXn];)

#if HAVE_ENERGY
    F[ENG]     = (U[ENG] + V[PRS])*V[VXn];
#elif EOS == ISOTHERMAL
    // Should we do anything here?
#endif
} 

KOKKOS_INLINE_FUNCTION void K_ConsToPrim(real Vc[], real Uc[], real gamma_m1) {
    

    Vc[RHO] = Uc[RHO];

    EXPAND( Vc[VX1] = Uc[MX1]/Uc[RHO]; ,
            Vc[VX2] = Uc[MX2]/Uc[RHO]; ,
            Vc[VX3] = Uc[MX3]/Uc[RHO];)     

#if HAVE_ENERGY
    real kin;
    kin = 0.5 / Uc[RHO] * (EXPAND(    Uc[MX1]*Uc[MX1] , 
                                    + Uc[MX2]*Uc[MX2] ,
                                    + Uc[MX3]*Uc[MX3] ));

    Vc[PRS] = gamma_m1 * (Uc[ENG] - kin);
#endif  // Have_energy
}

KOKKOS_INLINE_FUNCTION void K_PrimToCons(real Uc[], real Vc[], real gamma_m1) {

    Uc[RHO] = Vc[RHO];

    EXPAND( Uc[MX1] = Vc[VX1]*Vc[RHO]; ,
            Uc[MX2] = Vc[VX2]*Vc[RHO]; ,
            Uc[MX3] = Vc[VX3]*Vc[RHO];)

#if HAVE_ENERGY

    Uc[ENG] = Vc[PRS] / gamma_m1 
                + 0.5 * Vc[RHO] * EXPAND(  Vc[VX1]*Vc[VX1] , 
                                         + Vc[VX2]*Vc[VX2] ,
                                         + Vc[VX3]*Vc[VX3] ); 
#endif  // Have_energy

}




Physics::Physics(Input &input) {
    Kokkos::Profiling::pushRegion("Physics::Physics(DataBock)");

    this->gamma = 5.0/3.0;
    this->C2Iso = 1.0;

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

void Physics::InitFlow(DataBlock & data) {
    // Create a host copy
    DataBlockHost d(data);

    for(int dir=0; dir<3; dir++) {
        std::cout << "np["<<dir<<"]="<< data.np_tot[dir]<<std::endl;
    }
    for(int k = 0; k < d.np_tot[KDIR] ; k++) {
        for(int j = 0; j < d.np_tot[JDIR] ; j++) {
            for(int i = 0; i < d.np_tot[IDIR] ; i++) {
                d.Vc(RHO,k,j,i) = ONE_F;
                EXPAND(\
                d.Vc(VX1,k,j,i) = cos(2.0*M_PI*d.x[JDIR](j)); ,\
                d.Vc(VX2,k,j,i) = randm(); ,\
                d.Vc(VX3,k,j,i) = ZERO_F; )
#if HAVE_ENERGY 
                d.Vc(PRS,k,j,i) = ONE_F;
#endif
            }
        }
    }
    
    d.SyncToDevice();

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

    real gamma_m1=this->gamma-ONE_F;
    real C2Iso = this->C2Iso;

    idefix_for("CalcRiemannFlux",data.beg[KDIR],data.end[KDIR]+koffset,data.beg[JDIR],data.end[JDIR]+joffset,data.beg[IDIR],data.end[IDIR]+ioffset,
                        KOKKOS_LAMBDA (int k, int j, int i) 
            {
                int VXn = VX1+dir;
                int MXn = VXn;
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
                real c2RL, cmax;

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
                K_Flux(fluxL, vL, uL, dir);
                K_Flux(fluxR, vR, uR, dir);

                // 4-- Get the wave speed
            #if HAVE_ENERGY
                c2RL = (gamma_m1+ONE_F)*(vRL[PRS]/vRL[RHO]);
            #else
                c2RL = C2Iso;
            #endif
                cmax = FMAX(c2RL+vRL[VXn],c2RL-vRL[VXn]);

                // 5-- Compute the flux from the left and right states
                for(int nv = 0 ; nv < NVAR; nv++) {
                    Flux(nv,k,j,i) = HALF_F*(fluxL[nv]+fluxR[nv] - cmax*(uR[nv]-uL[nv]));
                }

                //6-- Compute maximum dt for this sweep
                const int ig = ioffset*i + joffset*j + koffset*k;

                invDt(k,j,i) = FMAX(invDt(k,j,i),cmax/dx(ig));

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
void Physics::SetBoundary(DataBlock &data) {

    Kokkos::Profiling::pushRegion("Physics::SetBoundary");

    IdefixArray4D<real> Vc = data.Vc;
    // X1 boundary conditions
    int offset = data.np_int[IDIR];
    int nghost = data.nghost[IDIR];

    idefix_for("Boundary_X1",0,NVAR,data.beg[KDIR],data.end[KDIR],data.beg[JDIR],data.end[JDIR],0,nghost,
        KOKKOS_LAMBDA (int n, int k, int j, int i) {
            
            Vc(n,k,j,i) = Vc(n,k,j,i+offset);
            Vc(n,k,j,i+offset+nghost) = Vc(n,k,j,i+nghost); 

        });

    if(DIMENSIONS>=2) {
        // X2 boundary conditions
        int offset = data.np_int[JDIR];
        int nghost = data.nghost[JDIR];

        idefix_for("Boundary_X2",0,NVAR,data.beg[KDIR],data.end[KDIR],0,nghost,0,data.np_tot[IDIR],
            KOKKOS_LAMBDA (int n, int k, int j, int i) {
                
                Vc(n,k,j,i) = Vc(n,k,j+offset,i);
                Vc(n,k,j+offset+nghost,i) = Vc(n,k,j+nghost,i); 

            });
    }

    if(DIMENSIONS==3) {
        // X3 boundary conditions
        int offset = data.np_int[KDIR];
        int nghost = data.nghost[KDIR];

        idefix_for("Boundary_X3",0,NVAR,0,nghost,0,data.np_tot[JDIR],0,data.np_tot[IDIR],
            KOKKOS_LAMBDA (int n, int k, int j, int i) {
                
                Vc(n,k,j,i) = Vc(n,k+offset,j,i);
                Vc(n,k+offset+nghost,j,i) = Vc(n,k+nghost,j,i); 

            });
    }

    Kokkos::Profiling::popRegion();

}

/*********************************************/
/**
Customized random number generator
Allow one to have consistant random numbers
generators on different architectures.
**/
/*********************************************/
real Physics::randm(void) {
	const int a	=	16807;
	const int m =	2147483647;
	static int in0 = 13763;
	int q;

	/* find random number  */
	q= (int) fmod((real) a * in0, m);
	in0=q;

	return((real)q/(real)m);
}
