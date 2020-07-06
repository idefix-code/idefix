#include "../idefix.hpp"
#include "hydro.hpp"

#if MHD == YES
#include "solversMHD.hpp"
#else
#include "solversHD.hpp"
#endif


Hydro::Hydro(Input &input, Grid &grid) {
    idfx::pushRegion("Hydro::Hydro(DataBock)");

    this->gamma = 5.0/3.0;
    this->C2Iso = 1.0;
    
    // read Solver from input file
    std::string solverString = input.GetString("Hydro","Solver",0);
    
    if (solverString.compare("tvdlf") == 0)     mySolver = TVDLF;
    else if (solverString.compare("hll") == 0)  mySolver = HLL;
    #if MHD == YES
        else if (solverString.compare("hlld") == 0) mySolver = HLLD;
    #else
        else if (solverString.compare("hllc") == 0) mySolver = HLLC;
    #endif
    else if (solverString.compare("roe") == 0)  mySolver = ROE;
    else {
        std::stringstream msg;
    #if MHD == YES
        msg << "Unknown MHD solver type " << solverString;
    #else
        msg << "Unknown HD solver type " << solverString;
    #endif
        IDEFIX_ERROR(msg);
    }

    // Source terms
    this->haveSourceTerms = false;

    // Check whether we have rotation
    int rotation = input.CheckEntry("Hydro","Rotation");

    if(rotation>=0 ) {
        this->haveSourceTerms = true;
        this->haveRotation = true;
        if(rotation != 3) IDEFIX_ERROR("Rotation needs a 3 components vector in idefix.ini");
        this->OmegaX1 = input.GetReal("Hydro","Rotation",0);
        this->OmegaX2 = input.GetReal("Hydro","Rotation",1);
        this->OmegaX3 = input.GetReal("Hydro","Rotation",2);

        idfx::cout << "Rotation enabled with Omega=(" << this->OmegaX1 << ", " << this->OmegaX2 << ", " << this->OmegaX3 << ")" << std::endl;
    }
    else {
        this->haveRotation = false;
    }

    // Check whether we have shearing box
    int shearingbox = input.CheckEntry("Hydro","ShearingBox");

    if(shearingbox>=0 ) {
        this->haveShearingBox = true;
        this->haveSourceTerms = true;
        if(shearingbox != 1) IDEFIX_ERROR("Shearing box needs a scalar value for the shear rate in idefix.ini");
        this->sbS = input.GetReal("Hydro","ShearingBox",0);

        // Get box size
        this->sbLx = grid.xend[IDIR] - grid.xbeg[IDIR];

        idfx::cout << "xbeg=" << grid.xbeg[IDIR] << "xend=" << grid.xend[IDIR] <<std::endl;
        idfx::cout << "ShearingBox enabled with Shear rate=" << this->sbS <<  "and Lx=" << sbLx << std::endl;
    }
    else {
        this->haveShearingBox = false;
    }

    idfx::popRegion();
}

Hydro::Hydro() {

}

real Hydro::GetGamma() {
    return(this->gamma);
}

real Hydro::GetC2iso() {
    return(this->C2Iso);
}


// Convect Conservative to Primitive variable
void Hydro::ConvertConsToPrim(DataBlock & data) {

    idfx::pushRegion("Hydro::ConvertConsToPrim");

    IdefixArray4D<real> Vc = data.Vc;
    IdefixArray4D<real> Uc = data.Uc;
    real gamma_m1=this->gamma-ONE_F;

    idefix_for("ConsToPrim",0,data.np_tot[KDIR],0,data.np_tot[JDIR],0,data.np_tot[IDIR],
                        KOKKOS_LAMBDA (int k, int j, int i) 
            {
                real U_RHO;
                EXPAND(
                    real U_MX1; ,
                    real U_MX2; ,
                    real U_MX3; )
#if MHD == YES
                EXPAND(
                    real U_BX1; ,
                    real U_BX2; ,
                    real U_BX3; )
#endif
                
                real V_RHO;
                EXPAND(
                    real V_VX1; ,
                    real V_VX2; ,
                    real V_VX3; )
#if MHD == YES
                EXPAND(
                    real V_BX1; ,
                    real V_BX2; ,
                    real V_BX3; )
#endif
                
#if HAVE_ENERGY
                real U_ENG = Uc(ENG,k,j,i);
                real V_PRS;
#endif

                U_RHO = Uc(RHO,k,j,i);
                EXPAND (
                    U_MX1 = Uc(MX1,k,j,i); ,
                    U_MX2 = Uc(MX2,k,j,i); ,
                    U_MX3 = Uc(MX3,k,j,i);
                )
#if MHD == YES
                EXPAND (
                    U_BX1 = Uc(BX1,k,j,i); ,
                    U_BX2 = Uc(BX2,k,j,i); ,
                    U_BX3 = Uc(BX3,k,j,i);
                )
                
                K_ConsToPrim(V_RHO, ARG_EXPAND(V_VX1, V_VX2, V_VX3, V_BX1, V_BX2, V_BX3, V_PRS),
                             U_RHO, ARG_EXPAND(U_MX1, U_MX2, U_MX3, U_BX1, U_BX2, U_BX3, U_ENG),
                             gamma_m1);
#else
                K_ConsToPrim(V_RHO, ARG_EXPAND(V_VX1, V_VX2, V_VX3, V_PRS),
                             U_RHO, ARG_EXPAND(U_MX1, U_MX2, U_MX3, U_ENG),
                             gamma_m1);
#endif

                Vc(RHO,k,j,i) = V_RHO;
                EXPAND (
                    Vc(VX1,k,j,i) = V_VX1; ,
                    Vc(VX2,k,j,i) = V_VX2; ,
                    Vc(VX3,k,j,i) = V_VX3;
                )
#if MHD == YES
                EXPAND (
                    Vc(BX1,k,j,i) = V_BX1; ,
                    Vc(BX2,k,j,i) = V_BX2; ,
                    Vc(BX3,k,j,i) = V_BX3;
                )
#endif
#if HAVE_ENERGY
                Vc(PRS,k,j,i) = V_PRS;
#endif

            });

    idfx::popRegion();

}


// Convert Primitive to conservative variables
void Hydro::ConvertPrimToCons(DataBlock & data) {

    idfx::pushRegion("Hydro::ConvertPrimToCons");

    IdefixArray4D<real> Vc = data.Vc;
    IdefixArray4D<real> Uc = data.Uc;
    real gamma_m1=this->gamma-ONE_F;

    idefix_for("ConvertPrimToCons",0,data.np_tot[KDIR],0,data.np_tot[JDIR],0,data.np_tot[IDIR],
                        KOKKOS_LAMBDA (int k, int j, int i) 
            {
                real U_RHO;
                EXPAND(
                    real U_MX1; ,
                    real U_MX2; ,
                    real U_MX3; )
#if MHD == YES
                EXPAND(
                    real U_BX1; ,
                    real U_BX2; ,
                    real U_BX3; )
#endif
                
                real V_RHO;
                EXPAND(
                    real V_VX1; ,
                    real V_VX2; ,
                    real V_VX3; )
#if MHD == YES
                EXPAND(
                    real V_BX1; ,
                    real V_BX2; ,
                    real V_BX3; )
#endif

#if HAVE_ENERGY
                real U_ENG;
                real V_PRS = Vc(PRS,k,j,i);
#endif
                
                V_RHO = Vc(RHO,k,j,i);
                EXPAND (
                    V_VX1 = Vc(VX1,k,j,i); ,
                    V_VX2 = Vc(VX2,k,j,i); ,
                    V_VX3 = Vc(VX3,k,j,i);
                )

#if MHD == YES
                EXPAND (
                    V_BX1 = Vc(BX1,k,j,i); ,
                    V_BX2 = Vc(BX2,k,j,i); ,
                    V_BX3 = Vc(BX3,k,j,i);
                )
                
                K_PrimToCons(U_RHO, ARG_EXPAND(U_MX1, U_MX2, U_MX3, U_BX1, U_BX2, U_BX3, U_ENG),
                             V_RHO, ARG_EXPAND(V_VX1, V_VX2, V_VX3, V_BX1, V_BX2, V_BX3, V_PRS),
                             gamma_m1);
#else
                K_PrimToCons(U_RHO, ARG_EXPAND(U_MX1, U_MX2, U_MX3, U_ENG),
                             V_RHO, ARG_EXPAND(V_VX1, V_VX2, V_VX3, V_PRS),
                             gamma_m1);
#endif
                

                Uc(RHO,k,j,i) = U_RHO;
                EXPAND (
                    Uc(MX1,k,j,i) = U_MX1; ,
                    Uc(MX2,k,j,i) = U_MX2; ,
                    Uc(MX3,k,j,i) = U_MX3;
                )
#if MHD == YES
                EXPAND (
                    Uc(BX1,k,j,i) = U_BX1; ,
                    Uc(BX2,k,j,i) = U_BX2; ,
                    Uc(BX3,k,j,i) = U_BX3;
                )
#endif
#if HAVE_ENERGY
                Uc(ENG,k,j,i) = U_ENG;
#endif
                
            });

    idfx::popRegion();
}


// Build a left and right extrapolation of the primitive variables along direction dir

// These functions extrapolate the cell prim vars to the faces. Definitions are as followed
//
// |       cell i-1               interface i          cell i
// |-----------------------------------|------------------------------------||
//          Vc(i-1)           PrimL(i)  PrimR(i)       Vc(i)

void Hydro::ExtrapolatePrimVar(DataBlock &data, int dir) {
    int ioffset,joffset,koffset;
    int iextend, jextend,kextend;
    int BXn;

    idfx::pushRegion("Hydro::ExtrapolatePrimVar");
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



    idfx::popRegion();
}

// Compute Riemann fluxes from states
void Hydro::CalcRiemannFlux(DataBlock & data, int dir) {

    idfx::pushRegion("Hydro::CalcRiemannFlux");
    
    switch (mySolver) {
    #if MHD == YES
        case TVDLF: TvdlfMHD(data, dir, this->gamma, this->C2Iso);
            break;
        case HLL:   HllMHD(data, dir, this->gamma, this->C2Iso);
            break;
        case HLLD:  HlldMHD(data, dir, this->gamma, this->C2Iso);
            break;
        case ROE:   RoeMHD(data, dir, this->gamma, this->C2Iso);
            break;
    #else
        case TVDLF: TvdlfHD(data, dir, this->gamma, this->C2Iso);
            break;
        case HLL:   HllHD(data, dir, this->gamma, this->C2Iso);
            break;
        case HLLC:  HllcHD(data, dir, this->gamma, this->C2Iso);
            break;
        case ROE:   RoeHD(data, dir, this->gamma, this->C2Iso);
            break;
    #endif
        default: // do nothing
            break;
    }

    idfx::popRegion();

}

// Compute the right handside in direction dir from conservative equation, with timestep dt
void Hydro::CalcRightHandSide(DataBlock &data, int dir, real dt) {

    idfx::pushRegion("Hydro::CalcRightHandSide");
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

    idfx::popRegion();
}

// Add source terms
void Hydro::AddSourceTerms(DataBlock &data, real dt) {

    idfx::pushRegion("Hydro::AddSourceTerms");
    IdefixArray4D<real> Uc = data.Uc;
    IdefixArray4D<real> Vc = data.Vc;

    real OmegaX1 = this->OmegaX1;
    real OmegaX2 = this->OmegaX2;
    real OmegaX3 = this->OmegaX3;
    bool haveRotation = this->haveRotation;
    
    idefix_for("AddSourceTerms",data.beg[KDIR],data.end[KDIR],data.beg[JDIR],data.end[JDIR],data.beg[IDIR],data.end[IDIR],
        KOKKOS_LAMBDA (int k, int j, int i) {
            if(haveRotation) {
                #if COMPONENTS == 3
                Uc(MX1,k,j,i) += TWO_F * dt * Vc(RHO,k,j,i) * (OmegaX3 * Vc(VX2,k,j,i) - OmegaX2 * Vc(VX3,k,j,i));
                Uc(MX2,k,j,i) += TWO_F * dt * Vc(RHO,k,j,i) * (OmegaX1 * Vc(VX3,k,j,i) - OmegaX3 * Vc(VX1,k,j,i));
                Uc(MX3,k,j,i) += TWO_F * dt * Vc(RHO,k,j,i) * (OmegaX2 * Vc(VX1,k,j,i) - OmegaX1 * Vc(VX2,k,j,i));
                #endif
                #if COMPONENTS == 2
                Uc(MX1,k,j,i) += TWO_F * dt * Vc(RHO,k,j,i) * (   OmegaX3 * Vc(VX2,k,j,i) );
                Uc(MX2,k,j,i) += TWO_F * dt * Vc(RHO,k,j,i) * ( - OmegaX3 * Vc(VX1,k,j,i) );
                #endif

            }
        });
    
    idfx::popRegion();
    
}

// Compute Corner EMFs from the one stored in the Riemann step
void Hydro::CalcCornerEMF(DataBlock &data, real t) {
        idfx::pushRegion("Hydro::CalcCornerEMF");

    // Corned EMFs
    IdefixArray3D<real> ex = data.emf.ex;
    IdefixArray3D<real> ey = data.emf.ey;
    IdefixArray3D<real> ez = data.emf.ez;

    // Face-centered EMFs
    IdefixArray3D<real> exj = data.emf.exj;
    IdefixArray3D<real> exk = data.emf.exk;
    IdefixArray3D<real> eyi = data.emf.eyi;
    IdefixArray3D<real> eyk = data.emf.eyk;
    IdefixArray3D<real> ezi = data.emf.ezi;
    IdefixArray3D<real> ezj = data.emf.ezj;

#if EMF_AVERAGE == ARITHMETIC
    idefix_for("CalcCornerEMF",
                data.beg[KDIR],data.end[KDIR]+KOFFSET,
                data.beg[JDIR],data.end[JDIR]+JOFFSET,
                data.beg[IDIR],data.end[IDIR]+IOFFSET,
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
#endif
    
#if EMF_AVERAGE == UCT_CONTACT || EMF_AVERAGE == UCT0

    // 0. Compute cell-centered emf.
    IdefixArray4D<real> Vc = data.Vc;
    
    IdefixArray3D<real> Ex1 = data.emf.Ex1;
    IdefixArray3D<real> Ex2 = data.emf.Ex2;
    IdefixArray3D<real> Ex3 = data.emf.Ex3;

    idefix_for("CalcCenterEMF",
                0,data.np_tot[KDIR],
                0,data.np_tot[JDIR],
                0,data.np_tot[IDIR],
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
                });
#endif

    // 1. averaging scheme
#if EMF_AVERAGE == UCT0
    idefix_for("CalcCornerEMF",
                data.beg[KDIR]-KOFFSET,data.end[KDIR]+KOFFSET,
                data.beg[JDIR]-JOFFSET,data.end[JDIR]+JOFFSET,
                data.beg[IDIR]-IOFFSET,data.end[IDIR]+IOFFSET,
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
                });
    
    idefix_for("CalcCornerEMF",
                data.beg[KDIR],data.end[KDIR]+KOFFSET,
                data.beg[JDIR],data.end[JDIR]+JOFFSET,
                data.beg[IDIR],data.end[IDIR]+IOFFSET,
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
#endif
    
#if EMF_AVERAGE == UCT_CONTACT
    IdefixArray3D<int> svx = data.emf.svx;
    IdefixArray3D<int> svy = data.emf.svy;
    IdefixArray3D<int> svz = data.emf.svz;
    
    idefix_for("EMF_ArithmeticAverage",
                data.beg[KDIR],data.end[KDIR]+KOFFSET,
                data.beg[JDIR],data.end[JDIR]+JOFFSET,
                data.beg[IDIR],data.end[IDIR]+IOFFSET,
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
    
    idefix_for("CalcCornerEMF",
                data.beg[KDIR],data.end[KDIR]+KOFFSET,
                data.beg[JDIR],data.end[JDIR]+JOFFSET,
                data.beg[IDIR],data.end[IDIR]+IOFFSET,
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

                    if (sx == 0) {
                        ez(k,j,i) += HALF_F*(ezj(k,j,i-1)   - Ex3(k,j,i-1) + ezj(k,j,i)   - Ex3(k,j,i));
                        ez(k,j,i) -= HALF_F*(Ex3(k,j+1,i-1) - ezj(k,j,i-1) + Ex3(k,j+1,i) - ezj(k,j,i));
    #if DIMENSIONS == 3
                        ey(k,j,i) += HALF_F*(eyk(k,j,i-1)   - Ex2(k,j,i-1) + eyk(k,j,i)   - Ex2(k,j,i));
                        ey(k,j,i) -= HALF_F*(Ex2(k+1,j,i-1) - eyk(k,j,i-1) + Ex2(k+1,j,i) - eyk(k,j,i));
    #endif
                    }
                    else {
                        ez(k,j,i) += ezj(k,j,iu)   - Ex3(k,j,iu);
                        ez(k,j,i) -= Ex3(k,j+1,iu) - ezj(k,j,iu);
    #if DIMENSIONS == 3
                        ey(k,j,i) += eyk(k,j,iu)   - Ex2(k,j,iu);
                        ey(k,j,i) -= Ex2(k+1,j,iu) - eyk(k,j,iu);
    #endif
                    }

                    // Span Y - Faces:    dEz/dx, dEx/dz

                    if (sy == 0) {
                        ez(k,j,i) += HALF_F*(ezi(k,j-1,i)   - Ex3(k,j-1,i) + ezi(k,j,i)   - Ex3(k,j,i));
                        ez(k,j,i) -= HALF_F*(Ex3(k,j-1,i+1) - ezi(k,j-1,i) + Ex3(k,j,i+1) - ezi(k,j,i));
    #if DIMENSIONS == 3
                        ex(k,j,i) += HALF_F*(exk(k,j-1,i)   - Ex1(k,j-1,i) + exk(k,j,i)   - Ex1(k,j,i));
                        ex(k,j,i) -= HALF_F*(Ex1(k+1,j-1,i) - exk(k,j-1,i) + Ex1(k+1,j,i) - exk(k,j,i));
    #endif
                    }
                    else {
                        ez(k,j,i) += ezi(k,ju,i)   - Ex3(k,ju,i);
                        ez(k,j,i) -= Ex3(k,ju,i+1) - ezi(k,ju,i);
    #if DIMENSIONS == 3
                        ex(k,j,i) += exk(k,ju,i)   - Ex1(k,ju,i);
                        ex(k,j,i) -= Ex1(k+1,ju,i) - exk(k,ju,i);
    #endif
                    }

                    // Span Z - Faces:    dEx/dy, dEy/dx

    #if DIMENSIONS == 3
                    if (sz == 0) {
                        ex(k,j,i) += HALF_F*(exj(k-1,j,i)   - Ex1(k-1,j,i) + exj(k,j,i)   - Ex1(k,j,i));
                        ex(k,j,i) -= HALF_F*(Ex1(k-1,j+1,i) - exj(k-1,j,i) + Ex1(k,j+1,i) - exj(k,j,i));
                        ey(k,j,i) += HALF_F*(eyi(k-1,j,i)   - Ex2(k-1,j,i) + eyi(k,j,i)   - Ex2(k,j,i));
                        ey(k,j,i) -= HALF_F*(Ex2(k-1,j,i+1) - eyi(k-1,j,i) + Ex2(k,j,i+1) - eyi(k,j,i));
                    }
                    else {
                        ex(k,j,i) += exj(ku,j,i)   - Ex1(ku,j,i);
                        ex(k,j,i) -= Ex1(ku,j+1,i) - exj(ku,j,i);
                        ey(k,j,i) += eyi(ku,j,i)   - Ex2(ku,j,i);
                        ey(k,j,i) -= Ex2(ku,j,i+1) - eyi(ku,j,i);
                    }

                    ex(k,j,i) *= ONE_FOURTH_F;
                    ey(k,j,i) *= ONE_FOURTH_F;
    #endif
                    ez(k,j,i) *= ONE_FOURTH_F;

                });
#endif

    idfx::popRegion();
}

// Evolve the magnetic field in Vs according to Constranied transport
void Hydro::EvolveMagField(DataBlock &data, real t, real dt) {
    idfx::pushRegion("Hydro::EvolveMagField");

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
    idfx::popRegion();
}

void Hydro::ReconstructVcField(DataBlock & data,  IdefixArray4D<real> &Vc) {
    idfx::pushRegion("Hydro::ReconstructVcField");
    IdefixArray4D<real> Vs=data.Vs;

    // Reconstruct cell average field when using CT
    idefix_for("ReconstructVcMagField",0,data.np_tot[KDIR],0,data.np_tot[JDIR],0,data.np_tot[IDIR],
                    KOKKOS_LAMBDA (int k, int j, int i) {
                        D_EXPAND(   Vc(BX1,k,j,i) = HALF_F * (Vs(BX1s,k,j,i) + Vs(BX1s,k,j,i+1)) ;  ,
                                    Vc(BX2,k,j,i) = HALF_F * (Vs(BX2s,k,j,i) + Vs(BX2s,k,j+1,i)) ;  ,
                                    Vc(BX3,k,j,i) = HALF_F * (Vs(BX3s,k,j,i) + Vs(BX3s,k+1,j,i)) ; )
                       
                    });
    idfx::popRegion();
}



void Hydro::ReconstructNormalField(DataBlock &data) {
    idfx::pushRegion("Hydro::ReconstructNormalField");

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
    nstart = data.beg[IDIR]-1;
    nend = data.end[IDIR];

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
    
    nstart = data.beg[JDIR]-1;
    nend = data.end[JDIR];
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
    
    nstart = data.beg[KDIR]-1;
    nend = data.end[KDIR];

    idefix_for("ReconstructBX3s",0,data.np_tot[JDIR],0,data.np_tot[IDIR],
                    KOKKOS_LAMBDA (int j, int i) {
                        for(int k = nstart ; k>=0 ; k-- ) {
                            Vs(BX3s,k,j,i) = Vs(BX3s,k+1,j,i) + dx3(k)*(     (Vs(BX1s,k,j,i+1) - Vs(BX1s,k,j,i))/dx1(i)                  
                                                                          +  (Vs(BX2s,k,j+1,i) - Vs(BX2s,k,j,i))/dx2(j) );
                        }
                        for(int k = nend ; k<nx3 ; k++ ) {
                            Vs(BX3s,k+1,j,i) = Vs(BX3s,k,j,i) -  dx3(k)*(     (Vs(BX1s,k,j,i+1) - Vs(BX1s,k,j,i))/dx1(i)                  
                                                                           +  (Vs(BX2s,k,j+1,i) - Vs(BX2s,k,j,i))/dx2(j) );
                        }
                        

                    });

    #endif  
    
    idfx::popRegion();
}


// Set Boundary conditions
void Hydro::SetBoundary(DataBlock &data, real t) {

    idfx::pushRegion("Hydro::SetBoundary");

    IdefixArray4D<real> Vc = data.Vc;
    IdefixArray4D<real> Vs = data.Vs;

    int ibeg,iend,jbeg,jend,kbeg,kend;
    int ioffset,joffset,koffset;
    int ighost,jghost,kghost;

    ighost = data.nghost[IDIR];
    jghost = data.nghost[JDIR];
    kghost = data.nghost[KDIR];

    real sbLx = this->sbLx;
    real sbS = this->sbS;

    // X1 boundary conditions
    

    for(int dir=0 ; dir < DIMENSIONS ; dir++ ) {
        // MPI Exchange data when needed
        if(data.mygrid->nproc[dir]>1) {
            switch(dir) {
                case 0:
                    data.ExchangeX1();
                    break;
                case 1:
                    data.ExchangeX2();
                    break;
                case 2:
                    data.ExchangeX3();
                    break;
            }
        }

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
            case internal:
                // internal is used for MPI-enforced boundary conditions. Nothing to be done here.
                break;
            case periodic:
                if(data.mygrid->nproc[dir] > 1) break; // Periodicity already enforced by MPI calls
                idefix_for("BoundaryBegPeriodic",0,NVAR,kbeg,kend,jbeg,jend,ibeg,iend,
                    KOKKOS_LAMBDA (int n, int k, int j, int i) {
                
                        Vc(n,k,j,i) = Vc(n,k+koffset,j+joffset,i+ioffset);
                    });
                #if MHD == YES
                idefix_for("BoundaryBegPeriodicVs",0,DIMENSIONS,kbeg,kend,jbeg,jend,ibeg,iend,
                    KOKKOS_LAMBDA (int n, int k, int j, int i) {

                        // Don't touch the normal component !
                        if(n != dir) Vs(n,k,j,i) = Vs(n,k+koffset,j+joffset,i+ioffset);
                    });
                #endif
                break;
            case outflow:
                idefix_for("BoundaryBegOutflow",0,NVAR,kbeg,kend,jbeg,jend,ibeg,iend,
                    KOKKOS_LAMBDA (int n, int k, int j, int i) {
                        int iref= (dir==IDIR) ? ighost : i;
                        int jref= (dir==JDIR) ? jghost : j;
                        int kref= (dir==KDIR) ? kghost : k;

                        if(n==VX1+dir) Vc(n,k,j,i) = ZERO_F;
                        else Vc(n,k,j,i) = Vc(n,kref,jref,iref);
                    });
                #if MHD == YES
                idefix_for("BoundaryBegOutflowVs",0,DIMENSIONS,kbeg,kend,jbeg,jend,ibeg,iend,
                    KOKKOS_LAMBDA (int n, int k, int j, int i) {
                        int iref= (dir==IDIR) ? ighost : i;
                        int jref= (dir==JDIR) ? jghost : j;
                        int kref= (dir==KDIR) ? kghost : k;

                        // Don't touch the normal component !
                        if(n != dir) Vs(n,k,j,i) = Vs(n,kref,jref,iref);
                    });
                #endif
                break;
            case shearingbox:
                if(data.mygrid->nproc[dir] > 1) {
                    // if shearing box enabled, the MPI call has already enforced strict periodicicty, so we just need to enforce the offset
                    real voffset=-sbLx*sbS;

                    idefix_for("BoundaryBegShearingBox",kbeg,kend,jbeg,jend,ibeg,iend,
                        KOKKOS_LAMBDA (int k, int j, int i) {
                            Vc(VX2,k,j,i) = Vc(VX2,k,j,i) + voffset;
                        });
                }
                else {
                    idefix_for("BoundaryBegShearingBox",0,NVAR,kbeg,kend,jbeg,jend,ibeg,iend,
                        KOKKOS_LAMBDA (int n, int k, int j, int i) {
                            real voffset= (n == VX2) ? - sbLx * sbS : ZERO_F;
                            Vc(n,k,j,i) = Vc(n,k+koffset,j+joffset,i+ioffset) + voffset;
                        });
                    #if MHD == YES
                    idefix_for("BoundaryBegShearingBoxVs",0,DIMENSIONS,kbeg,kend,jbeg,jend,ibeg,iend,
                        KOKKOS_LAMBDA (int n, int k, int j, int i) {

                            // Don't touch the normal component !
                            if(n != dir) Vs(n,k,j,i) = Vs(n,k+koffset,j+joffset,i+ioffset);
                        });
                    #endif
                }
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
            case internal:
                // internal is used for MPI-enforced boundary conditions. Nothing to be done here.
                break;
            case periodic:
                if(data.mygrid->nproc[dir] > 1) break; // Periodicity already enforced by MPI calls
                idefix_for("BoundaryEndPeriodic",0,NVAR,kbeg,kend,jbeg,jend,ibeg,iend,
                    KOKKOS_LAMBDA (int n, int k, int j, int i) {
                
                        Vc(n,k,j,i) = Vc(n,k-koffset,j-joffset,i-ioffset);
                    });
                #if MHD == YES
                idefix_for("BoundaryEndPeriodicVs",0,DIMENSIONS,kbeg,kend,jbeg,jend,ibeg,iend,
                    KOKKOS_LAMBDA (int n, int k, int j, int i) {
                        // Don't touch the normal component !
                        if(n != dir) Vs(n,k,j,i) = Vs(n,k-koffset,j-joffset,i-ioffset);                        
                    });
                #endif
                break;
            case outflow:
                idefix_for("BoundaryEndOutflow",0,NVAR,kbeg,kend,jbeg,jend,ibeg,iend,
                    KOKKOS_LAMBDA (int n, int k, int j, int i) {
                        int iref= (dir==IDIR) ? ighost + ioffset - 1 : i;
                        int jref= (dir==JDIR) ? jghost + joffset - 1 : j;
                        int kref= (dir==KDIR) ? kghost + koffset - 1 : k;

                        if(n==VX1+dir) Vc(n,k,j,i) = ZERO_F;
                        else Vc(n,k,j,i) = Vc(n,kref,jref,iref);
                    });
                #if MHD == YES
                idefix_for("BoundaryEndOutflowVs",0,DIMENSIONS,kbeg,kend,jbeg,jend,ibeg,iend,
                    KOKKOS_LAMBDA (int n, int k, int j, int i) {
                        int iref= (dir==IDIR) ? ighost + ioffset - 1 : i;
                        int jref= (dir==JDIR) ? jghost + joffset - 1 : j;
                        int kref= (dir==KDIR) ? kghost + koffset - 1 : k;

                        if(n != dir) Vs(n,k,j,i) = Vs(n,kref,jref,iref);
                    });
                #endif
                break;
            case shearingbox:
                if(data.mygrid->nproc[dir] > 1) {
                    // if shearing box enabled, the MPI call has already enforced strict periodicicty, so we just need to enforce the offset
                    real voffset=sbLx*sbS;

                    idefix_for("BoundaryEndShearingBox",kbeg,kend,jbeg,jend,ibeg,iend,
                        KOKKOS_LAMBDA (int k, int j, int i) {
                            Vc(VX2,k,j,i) = Vc(VX2,k,j,i) + voffset;
                        });
                }
                else {
                    idefix_for("BoundaryEndShearingBox",0,NVAR,kbeg,kend,jbeg,jend,ibeg,iend,
                        KOKKOS_LAMBDA (int n, int k, int j, int i) {
                            real voffset= (n == VX2) ? + sbLx * sbS : ZERO_F;

                            Vc(n,k,j,i) = Vc(n,k-koffset,j-joffset,i-ioffset) + voffset;
                        });
                    #if MHD == YES
                    idefix_for("BoundaryEndShearingBoxVs",0,DIMENSIONS,kbeg,kend,jbeg,jend,ibeg,iend,
                        KOKKOS_LAMBDA (int n, int k, int j, int i) {
                            // Don't touch the normal component !
                            if(n != dir) Vs(n,k,j,i) = Vs(n,k-koffset,j-joffset,i-ioffset);                        
                        });
                    #endif
                }
                break;
            default:
                std::stringstream msg("Boundary condition type is not yet implemented");
                IDEFIX_ERROR(msg);

        }
    }   // Loop on dimension ends

    #if MHD == YES
    // Reconstruct the normal field component when using CT
    ReconstructNormalField(data);

    // Remake the cell-centered field.
    ReconstructVcField(data, data.Vc);
    #endif

    idfx::popRegion();

}



real Hydro::CheckDivB(DataBlock &data) {
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

    #ifdef WITH_MPI
    if(idfx::psize>1) {
        MPI_Allreduce(MPI_IN_PLACE, &divB, 1, realMPI, MPI_MAX, MPI_COMM_WORLD);
    }
    #endif
    return(divB);
}

/*
real Hydro::CheckDivB(DataBlock &data) {

    real divB=0;
    IdefixArray4D<real> Vs = data.Vs;
    IdefixArray1D<real> dx1 = data.dx[IDIR];
    IdefixArray1D<real> dx2 = data.dx[JDIR];
    IdefixArray1D<real> dx3 = data.dx[KDIR];

    int iref,jref,kref;
    
    for(int k = data.beg[KDIR] ; k < data.end[KDIR] ; k++) {
        for(int j = data.beg[JDIR] ; j < data.end[JDIR] ; j++) {
            for(int i = data.beg[IDIR] ; i < data.end[IDIR] ; i++) {
                real dB1,dB2,dB3;

                dB1=dB2=dB3=ZERO_F;

                D_EXPAND( dB1=(Vs(BX1s,k,j,i+1)-Vs(BX1s,k,j,i))/(dx1(i)); ,
                          dB2=(Vs(BX2s,k,j+1,i)-Vs(BX2s,k,j,i))/(dx2(j)); ,
                          dB3=(Vs(BX3s,k+1,j,i)-Vs(BX3s,k,j,i))/(dx3(k));  )
                
                if(FABS(D_EXPAND(dB1, +dB2, +dB3)) > divB) {
                    iref=i;
                    jref=j;
                    kref=k;
                    divB=FABS(D_EXPAND(dB1, +dB2, +dB3));
                }
            }
        }
    }
    //idfx::cout << "divB=" << divB << "(i,j,k)=(" << iref << "," << jref << "," << kref << ")" << std::endl;
    return(divB);

}

*/


void Hydro::SetGamma(real newGamma) {
    this->gamma=newGamma;
}
