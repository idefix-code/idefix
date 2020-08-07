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
    
    // No userdefBoundary by default
    this->haveUserDefBoundary = false;

    // Source terms (always activated when non-cartesian geometry because of curvature source terms)
    #if GEOMETRY == CARTESIAN
    this->haveSourceTerms = false;
    #else
    this->haveSourceTerms = true;
    #endif

    // Check whether we have rotation
    int rotation = input.CheckEntry("Hydro","Rotation");

    if(rotation>=0 ) {
        this->haveSourceTerms = true;
        this->haveRotation = true;
        if(rotation != 3) IDEFIX_ERROR("Rotation needs a 3 components vector in idefix.ini");
        this->OmegaX1 = input.GetReal("Hydro","Rotation",0);
        this->OmegaX2 = input.GetReal("Hydro","Rotation",1);
        this->OmegaX3 = input.GetReal("Hydro","Rotation",2);

        idfx::cout << "Hydro: Rotation enabled with Omega=(" << this->OmegaX1 << ", " << this->OmegaX2 << ", " << this->OmegaX3 << ")" << std::endl;
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

        idfx::cout << "Hydro: ShearingBox enabled with Shear rate=" << this->sbS <<  "and Lx=" << sbLx << std::endl;
    }
    else {
        this->haveShearingBox = false;
    }

    // Gravitational potential
    this->haveGravPotential = false;
    this->gravPotentialFunc = nullptr;
    int gravPotential = input.CheckEntry("Hydro","GravPotential");
    if(gravPotential>=0) {
        std::string potentialString = input.GetString("Hydro","GravPotential",0);
        if(potentialString.compare("userdef") == 0) {
            this->haveGravPotential = true;        
            idfx::cout << "Hydro:: Enabling user-defined gravitational potential" << std::endl;
        }
        else {
            IDEFIX_ERROR("Unknown type of gravitational potential in idefix.ini. Only userdef is implemented");
        }
    }

    idfx::popRegion();
}

void Hydro::EnrollUserDefBoundary(UserDefBoundaryFunc myFunc) {
    this->userDefBoundaryFunc = myFunc;
    this->haveUserDefBoundary = true;
    idfx::cout << "Hydro: User-defined boundary condition has been enrolled" << std::endl;
}

void Hydro::EnrollGravPotential(GravPotentialFunc myFunc) {
    if(!this->haveGravPotential) IDEFIX_ERROR("In order to enroll your gravitational potential, you need to enable it first in the .ini file.");
    this->gravPotentialFunc = myFunc;
    idfx::cout << "Hydro: User-defined gravitational potential has been enrolled" << std::endl;
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
void Hydro::CalcRightHandSide(DataBlock &data, int dir, real t, real dt) {

    idfx::pushRegion("Hydro::CalcRightHandSide");
    IdefixArray4D<real> Uc = data.Uc;
    IdefixArray4D<real> Vc = data.Vc;
    IdefixArray4D<real> Flux = data.FluxRiemann;
    IdefixArray3D<real> A = data.A[dir];
    IdefixArray3D<real> dV = data.dV;
    IdefixArray1D<real> x1m = data.xl[IDIR];
    IdefixArray1D<real> x1 = data.x[IDIR];
    IdefixArray1D<real> sm = data.sm;
    IdefixArray1D<real> rt = data.rt;
    IdefixArray1D<real> dmu = data.dmu;
    IdefixArray1D<real> s = data.s;
    IdefixArray1D<real> dx = data.dx[dir];
    IdefixArray1D<real> dx2 = data.dx[JDIR];
    
    // Gravitational potential
    IdefixArray3D<real> phiP = data.phiP;
    bool needPotential = this->haveGravPotential;

    if(needPotential) {
        IdefixArray1D<real> x1,x2,x3;

        if(dir==IDIR) x1 = data.xl[IDIR];
        else          x1 = data.x[IDIR];
        if(dir==JDIR) x2 = data.xl[JDIR];
        else          x2 = data.x[JDIR];
        if(dir==KDIR) x3 = data.xl[KDIR];
        else          x3 = data.x[KDIR];
        
        if(this->gravPotentialFunc == nullptr) IDEFIX_ERROR("Gravitational potential is enabled, but no user-defined potential has been enrolled.");

        gravPotentialFunc(data, t, x1, x2, x3, phiP);
    }

    

    int ioffset,joffset,koffset;
    ioffset=joffset=koffset=0;
    // Determine the offset along which we do the extrapolation
    if(dir==IDIR) ioffset=1;
    if(dir==JDIR) joffset=1;
    if(dir==KDIR) koffset=1;
    


    idefix_for("CalcTotalFlux",data.beg[KDIR],data.end[KDIR]+koffset,data.beg[JDIR],data.end[JDIR]+joffset,data.beg[IDIR],data.end[IDIR]+ioffset,
            KOKKOS_LAMBDA (int k, int j, int i) {
                
                // TODO: Should add gravitational potential here and Fargo source terms when needed
                #if HAVE_ENERGY
                  if(needPotential) Flux(ENG, k, j, i) += Flux(RHO, k, j, i) * phiP(k,j,i) ;  // Potential at the cell face
                #endif


                for(int nv = 0 ; nv < NVAR ; nv++) {
                    Flux(nv,k,j,i) = Flux(nv,k,j,i) * A(k,j,i);
                }

                
                // Curvature terms
                #if    (GEOMETRY == POLAR       && COMPONENTS >= 2) \
                        || (GEOMETRY == CYLINDRICAL && COMPONENTS == 3) 
                    if(dir==IDIR) {
                        Flux(iMPHI,k,j,i) = Flux(iMPHI,k,j,i) * FABS(x1m(i));   // Conserve angular momentum, hence flux is R*Bphi
                        #if MHD == YES
                        Flux(iBPHI,k,j,i) = Flux(iBPHI,k,j,i) / A(k,j,i);   // No area for this one
                        #endif // MHD
                    }
                #endif // GEOMETRY==POLAR OR CYLINDRICAL

                #if GEOMETRY == SPHERICAL
                    if(dir==IDIR) {
                        #if COMPONENTS == 3
                            Flux(iMPHI,k,j,i) = Flux(iMPHI,k,j,i) * FABS(x1m(i));
                        #endif // COMPONENTS == 3
                        #if MHD == YES
                            EXPAND(                                            ,
                                Flux(iBTH,k,j,i)  = Flux(iBTH,k,j,i) * x1m(i) / A(k,j,i);  ,
                                Flux(iBPHI,k,j,i) = Flux(iBPHI,k,j,i) * x1m(i) / A(k,j,i); )
                        #endif // MHD
                    }
                       
                    if(dir==JDIR) {
                        #if COMPONENTS == 3  
                            Flux(iMPHI,k,j,i) = Flux(iMPHI,k,j,i) * FABS(sm(j));
                            #if PHYSICS == MHD
                                Flux(iBPHI,k,j,i) = Flux(iBPHI,k,j,i)  / A(k,j,i);
                            #endif // MHD
                        #endif // COMPONENTS = 3
                    }
                #endif // GEOMETRY == SPHERICAL


            });

    idefix_for("CalcRightHandSide",data.beg[KDIR],data.end[KDIR],data.beg[JDIR],data.end[JDIR],data.beg[IDIR],data.end[IDIR],
        KOKKOS_LAMBDA (int k, int j, int i) {
            
            real dtdV=dt / dV(k,j,i);
            real rhs[NVAR];

            for(int nv = 0 ; nv < NVAR ; nv++) {
                rhs[nv] = -  dtdV*(Flux(nv, k+koffset, j+joffset, i+ioffset) - Flux(nv, k, j, i));
            }

            #if GEOMETRY != CARTESIAN
                if(dir==IDIR) {
                    #ifdef iMPHI
                        rhs[iMPHI] = rhs[iMPHI] / x1(i); 
                    #endif
                    #if (GEOMETRY == POLAR || GEOMETRY == CYLINDRICAL) &&  (defined iBPHI)
                        rhs[iBPHI] = - dt / dx(i) * (Flux(iBPHI, k, j, i+1) - Flux(iBPHI, k, j, i) );

                    #elif (GEOMETRY == SPHERICAL) && (PHYSICS == MHD)
                        real q = dt / (x1(i)*dx(i));
                        EXPAND(                                                                     ,
                                rhs[iBTH]  = -q * ((Flux(iBTH, k, j, i+1)  - Flux(iBTH, k, j, i) ));  ,
                                rhs[iBPHI] = -q * ((Flux(iBPHI, k, j, i+1) - Flux(iBPHI, k, j, i) )); )
                    #endif
                }
                if(dir==JDIR) {
                    #if (GEOMETRY == SPHERICAL) && (COMPONENTS == 3)
                        rhs[iMPHI] /= FABS(s(j));
                        #if PHYSICS == MHD
                            rhs[iBPHI] = -dt / (rt(i)*dx(j)) * (Flux(iBPHI, k, j+1, i) - Flux(iBPHI, k, j, i));
                        #endif // MHD
                    #endif // GEOMETRY
                }
            // Nothing for KDIR

            #endif // GEOMETRY != CARTESIAN

            // Potential terms
            if(needPotential) {
                const int ig = ioffset*i + joffset*j + koffset*k;
                real dl = dx(ig);
                #if GEOMETRY == POLAR
                    if(dir==JDIR) dl = dl*x1(i);
                #elif GEOMETRY == SPHERICAL
                    if(dir==JDIR) dl = dl*rt(i);
                    if(dir==KDIR) dl = dl*rt(i)*dmu(j)/dx2(j);
                #endif
                rhs[1+dir] -= dt/dl * Vc(RHO,k,j,i) * (phiP(k+koffset,j+joffset,i+ioffset) - phiP(k,j,i));      // Gravitational force in direction i
                #if HAVE_ENERGY
                    rhs[ENG] -=  HALF_F * (phiP(k+koffset,j+joffset,i+ioffset) + phiP(k,j,i)) * rhs[RHO];        // We conserve total energy without potential
                #endif
            }            
            // Evolve the field components
            for(int nv = 0 ; nv < NVAR ; nv++) {
                // Do not evolve the field components if they are computed by CT (i.e. if they are in Vs)
                
                D_EXPAND( if(nv == BX1) continue;   ,
                          if(nv == BX2) continue;   ,
                          if(nv == BX3) continue;  ) 
                

                Uc(nv,k,j,i) = Uc(nv,k,j,i) + rhs[nv];
            }
            

        });

    idfx::popRegion();
}

// Add source terms
void Hydro::AddSourceTerms(DataBlock &data, real t, real dt) {

    idfx::pushRegion("Hydro::AddSourceTerms");
    IdefixArray4D<real> Uc = data.Uc;
    IdefixArray4D<real> Vc = data.Vc;
    IdefixArray1D<real> x1 = data.x[IDIR]; 
    IdefixArray1D<real> x2 = data.x[JDIR]; 
    #if GEOMETRY == SPHERICAL
    IdefixArray1D<real> s  = data.s;
    IdefixArray1D<real> rt = data.rt;
    #endif

    real OmegaX1 = this->OmegaX1;
    real OmegaX2 = this->OmegaX2;
    real OmegaX3 = this->OmegaX3;
    bool haveRotation = this->haveRotation;
    
    idefix_for("AddSourceTerms",data.beg[KDIR],data.end[KDIR],data.beg[JDIR],data.end[JDIR],data.beg[IDIR],data.end[IDIR],
        KOKKOS_LAMBDA (int k, int j, int i) {

            #if GEOMETRY == CARTESIAN
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
            #elif GEOMETRY == CYLINDRICAL
                #if COMPONENTS == 3
                    real vphi,Sm;
                    vphi = Vc(iVPHI,k,j,i);
                    if(haveRotation) vphi += OmegaX3*x1(i);
                    Sm = Vc(RHO,k,j,i) * vphi*vphi; // Centrifugal     
                    Sm += Vc(PRS,k,j,i);            // Presure (because pressure is included in the flux, additional source terms arise)
                    #if MHD==YES
                        Sm -=  Vc(iBPHI,k,j,i)*Vc(iBPHI,k,j,i); // Hoop stress
                        Sm += HALF_F*(EXPAND(Vc(BX1,k,j,i)*Vc(BX1,k,j,i) , +Vc(BX2,k,j,i)*Vc(BX2,k,j,i), +Vc(BX3,k,j,i)*Vc(BX3,k,j,i))); // Magnetic pressure
                    #endif // MHD
                    Uc(MX1,k,j,i) += dt * Sm / x1(i);
                #endif // COMPONENTS
            #elif GEOMETRY == POLAR
                real vphi,Sm;
                vphi = Vc(iVPHI,k,j,i);
                if(haveRotation) vphi += OmegaX3*x1(i);
                Sm = Vc(RHO,k,j,i) * vphi*vphi;     // Centrifugal
                Sm += Vc(PRS,k,j,i);               // Pressure (because we're including pressure in the flux, we need that to get the radial pressure gradient)       
                #if MHD==YES
                    Sm -=  Vc(iBPHI,k,j,i)*Vc(iBPHI,k,j,i); // Hoop stress
                    // Magnetic pressus
                    Sm += HALF_F*(EXPAND(Vc(BX1,k,j,i)*Vc(BX1,k,j,i) , +Vc(BX2,k,j,i)*Vc(BX2,k,j,i), +Vc(BX3,k,j,i)*Vc(BX3,k,j,i)));
                #endif // MHD
                Uc(MX1,k,j,i) += dt * Sm / x1(i);

            #elif GEOMETRY == SPHERICAL
                real vphi,Sm,ct; 
                vphi = SELECT(ZERO_F, ZERO_F, Vc(iVPHI,k,j,i));
                if(haveRotation) vphi += OmegaX3*x1(i)*s(j);
                Sm = Vc(RHO,k,j,i) * (EXPAND( ZERO_F, + Vc(VX2,k,j,i)*Vc(VX2,k,j,i), + vphi*vphi)); // Centrifugal
                Sm += 2.0*Vc(PRS,k,j,i);    // Pressure curvature
                #if MHD == YES
                    Sm -= EXPAND( ZERO_F, + Vc(iBTH,k,j,i)*Vc(iBTH,k,j,i), + Vc(iBPHI,k,j,i)*Vc(iBPHI,k,j,i)); // Hoop stress
                    Sm += EXPAND(Vc(BX1,k,j,i)*Vc(BX1,k,j,i) , +Vc(BX2,k,j,i)*Vc(BX2,k,j,i), +Vc(BX3,k,j,i)*Vc(BX3,k,j,i)); // 2* mag pressure curvature
                #endif
                Uc(MX1,k,j,i) += dt*Sm/x1(i);
                
                ct = 1.0/TAN(x2(j));
                Sm = Vc(RHO,k,j,i) * (EXPAND( ZERO_F, - Vc(iVTH,k,j,i)*Vc(iVR,k,j,i), + ct*vphi*vphi)); // Centrifugal
                Sm += ct * Vc(PRS,k,j,i);       // Pressure curvature
                #if MHD == YES
                    Sm += EXPAND( ZERO_F, + Vc(iBTH,k,j,i)*Vc(iBR,k,j,i), - ct*Vc(iBPHI,k,j,i)*Vc(iBPHI,k,j,i)); // Hoop stress
                    Sm += HALF_F*ct*EXPAND(Vc(BX1,k,j,i)*Vc(BX1,k,j,i) , +Vc(BX2,k,j,i)*Vc(BX2,k,j,i), +Vc(BX3,k,j,i)*Vc(BX3,k,j,i)); // Magnetic pressure
                #endif
                Uc(MX2,k,j,i) += dt*Sm / rt(i);
            #endif

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
                    }
                    else {
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
                    }
                    else {
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
                    }
                    else {
                        ex(k,j,i) += exj(ku,j,i) - Ex1(ku,j-1,i);
                        ex(k,j,i) -= Ex1(ku,j,i) - exj(ku,j,i);
                        ey(k,j,i) += eyi(ku,j,i) - Ex2(ku,j,i-1);
                        ey(k,j,i) -= Ex2(ku,j,i) - eyi(ku,j,i);
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

    IdefixArray1D<real> x1p=data.xr[IDIR];
    IdefixArray1D<real> x2p=data.xr[JDIR];
    IdefixArray1D<real> x3p=data.xr[KDIR];

    IdefixArray1D<real> x1m=data.xl[IDIR];
    IdefixArray1D<real> x2m=data.xl[JDIR];
    IdefixArray1D<real> x3m=data.xl[KDIR];

    IdefixArray1D<real> dx1=data.dx[IDIR];
    IdefixArray1D<real> dx2=data.dx[JDIR];
    IdefixArray1D<real> dx3=data.dx[KDIR];



    idefix_for("EvolvMagField",data.beg[KDIR],data.end[KDIR]+KOFFSET,data.beg[JDIR],data.end[JDIR]+JOFFSET,data.beg[IDIR],data.end[IDIR]+IOFFSET,
                    KOKKOS_LAMBDA (int k, int j, int i) {
                        
                        real rhsx1, rhsx2, rhsx3;

                        #if GEOMETRY == CARTESIAN
                            rhsx1 = D_EXPAND( ZERO_F                                     ,
                                            - dt/dx2(j) * (Ex3(k,j+1,i) - Ex3(k,j,i) )  ,
                                            + dt/dx3(k) * (Ex2(k+1,j,i) - Ex2(k,j,i) ) );
                                            
                            #if DIMENSIONS >= 2
                            rhsx2 =  D_EXPAND(   dt/dx1(i) * (Ex3(k,j,i+1) - Ex3(k,j,i) )  ,
                                                                                            ,
                                                - dt/dx3(k) * (Ex1(k+1,j,i) - Ex1(k,j,i) ) );
                            #endif
                            #if DIMENSIONS == 3
                            rhsx3 =    - dt/dx1(i) * (Ex2(k,j,i+1) - Ex2(k,j,i) )  
                                       + dt/dx2(j) * (Ex1(k,j+1,i) - Ex1(k,j,i) ) );
                            #endif

                        #elif GEOMETRY == CYLINDRICAL
                            rhsx1 = - dt/dx2(j) * (Ex3(k,j+1,i) - Ex3(k,j,i) );
                            #if DIMENSIONS >= 2
                            rhsx2 = dt * (FABS(x1p(i)) * Ex3(k,j,i+1) - FABS(x1m(i)) * Ex3(k,j,i)) / FABS(x1(i)*dx1(i));
                            #endif

                        #elif GEOMETRY == POLAR

                            rhsx1 = D_EXPAND( ZERO_F                                     ,
                                            - dt/(FABS(x1m(i)) * dx2(j)) * (Ex3(k,j+1,i) - Ex3(k,j,i) )  ,
                                            + dt/dx3(k) * (Ex2(k+1,j,i) - Ex2(k,j,i) ) );

                            #if DIMENSIONS >= 2
                            rhsx2 =  D_EXPAND(   dt/dx1(i) * (Ex3(k,j,i+1) - Ex3(k,j,i) )  ,
                                                                                            ,
                                                - dt/dx3(k) * (Ex1(k+1,j,i) - Ex1(k,j,i) ) );
                            #endif
                            #if DIMENSIONS == 3
                            rhsx3 =    dt/(FABS(x1(i))) * (
                                            -  (x1m(i+1)*Ex2(k,j,i+1) - x1m(i)*Ex2(k,j,i) ) / dx1(i) 
                                            +  (Ex1(k,j+1,i) - Ex1(k,j,i) ) / dx2(j) );
                            #endif

                        #elif GEOMETRY == SPHERICAL
                            real dV2 = FABS(cos(x2m(j)) - cos(x2p(j)));
                            real Ax2p = FABS(sin(x2p(j)));
                            real Ax2m = FABS(sin(x2m(j)));

                            rhsx1 = D_EXPAND( ZERO_F                                     ,
                                            - dt/(x1p(i)*dV2) * ( Ax2p*Ex3(k,j+1,i) - Ax2m*Ex3(k,j,i) )  ,
                                            + dt*dx2(j)/(x1m(i)*dV2*dx3(k))*(Ex2(k+1,j,i) - Ex2(k,j,i) ) );
                                            
                            #if DIMENSIONS >= 2
                            rhsx2 =  D_EXPAND(   dt/(x1(i)*dx1(i)) * (x1m(i+1)*Ex3(k,j,i+1) - x1m(i)*Ex3(k,j,i) )  ,
                                                                                            ,
                                                - dt/(x1(i)*Ax2m*dx3(k)) * (Ex1(k+1,j,i) - Ex1(k,j,i) ) );
                            #endif
                            #if DIMENSIONS == 3
                            rhsx3 =    - dt/(x1(i)*dx1(i)) * (x1m(i+1)*Ex2(k,j,i+1) - x1m(i)*Ex2(k,j,i) )  
                                       + dt/(x1(i)*dx2(j)) * (Ex1(k,j+1,i) - Ex1(k,j,i) ) );                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     * (Ex2(k+1,j,i) - Ex2(k,j,i) ) );
                            #endif
                        #endif



                        Vs(BX1s,k,j,i) = Vs(BX1s,k,j,i) + rhsx1;

                        #if DIMENSIONS >= 2
                        Vs(BX2s,k,j,i) = Vs(BX2s,k,j,i) + rhsx2;
                        #endif
                        #if DIMENSIONS == 3
                        Vs(BX3s,k,j,i) = Vs(BX3s,k,j,i) + rhsx3;
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



void Hydro::ReconstructNormalField(DataBlock &data, int dir) {
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

    IdefixArray3D<real> Ax1=data.A[IDIR];
    IdefixArray3D<real> Ax2=data.A[JDIR];
    IdefixArray3D<real> Ax3=data.A[KDIR];

    int nstart, nend;
    int nx1,nx2,nx3;

    // reconstruct BX1s
    nstart = data.beg[IDIR]-1;
    nend = data.end[IDIR];

    nx1=data.np_tot[IDIR];
    nx2=data.np_tot[JDIR];
    nx3=data.np_tot[KDIR];

    if(dir==IDIR) {
        idefix_for("ReconstructBX1s",0,nx3,0,nx2,
                        KOKKOS_LAMBDA (int k, int j) {


                            for(int i = nstart ; i>=0 ; i-- ) {
                                Vs(BX1s,k,j,i) = 1/ Ax1(k,j,i) * (   Ax1(k,j,i+1)*Vs(BX1s,k,j,i+1)  +   (D_EXPAND( ZERO_F                                       ,                    
                                                                                                +  Ax2(k,j+1,i) * Vs(BX2s,k,j+1,i) - Ax2(k,j,i) * Vs(BX2s,k,j,i)  , 
                                                                                                +  Ax3(k+1,j,i) * Vs(BX3s,k+1,j,i) - Ax3(k,j,i) * Vs(BX3s,k,j,i) ))) ;
                            }

                            for(int i = nend ; i<nx1 ; i++ ) {
                                Vs(BX1s,k,j,i+1) = 1/ Ax1(k,j,i+1) * (   Ax1(k,j,i)*Vs(BX1s,k,j,i)  -   (D_EXPAND(      ZERO_F                                       ,                    
                                                                                                +  Ax2(k,j+1,i) * Vs(BX2s,k,j+1,i) - Ax2(k,j,i) * Vs(BX2s,k,j,i)  , 
                                                                                                +  Ax3(k+1,j,i) * Vs(BX3s,k+1,j,i) - Ax3(k,j,i) * Vs(BX3s,k,j,i) ))) ;
                            }
                            

                        });
    }

    #if DIMENSIONS >=2
    
    if(dir==JDIR) {
        nstart = data.beg[JDIR]-1;
        nend = data.end[JDIR];
        idefix_for("ReconstructBX2s",0,data.np_tot[KDIR],0,data.np_tot[IDIR],
                        KOKKOS_LAMBDA (int k, int i) {
                            for(int j = nstart ; j>=0 ; j-- ) {
                                Vs(BX2s,k,j,i) = 1/ Ax2(k,j,i) * (   Ax2(k,j+1,i)*Vs(BX2s,k,j+1,i)  +   (D_EXPAND( Ax1(k,j,i+1) * Vs(BX1s,k,j,i+1) - Ax1(k,j,i) * Vs(BX1s,k,j,i)  , 
                                                                                                                                                                                , 
                                                                                                                +  Ax3(k+1,j,i) * Vs(BX3s,k+1,j,i) - Ax3(k,j,i) * Vs(BX3s,k,j,i) ))) ;
                            }
                            for(int j = nend ; j<nx2 ; j++ ) {
                                Vs(BX2s,k,j+1,i) = 1/ Ax2(k,j+1,i) * (   Ax2(k,j,i)*Vs(BX2s,k,j,i)  -   (D_EXPAND( Ax1(k,j,i+1) * Vs(BX1s,k,j,i+1) - Ax1(k,j,i) * Vs(BX1s,k,j,i)  , 
                                                                                                                                                                                , 
                                                                                                                +  Ax3(k+1,j,i) * Vs(BX3s,k+1,j,i) - Ax3(k,j,i) * Vs(BX3s,k,j,i) ))) ;
                            }
                            

                        });
        #endif
    }
    #if DIMENSIONS == 3
    
    if(dir==KDIR) {
        nstart = data.beg[KDIR]-1;
        nend = data.end[KDIR];

        idefix_for("ReconstructBX3s",0,data.np_tot[JDIR],0,data.np_tot[IDIR],
                        KOKKOS_LAMBDA (int j, int i) {
                            for(int k = nstart ; k>=0 ; k-- ) {
                                Vs(BX3s,k,j,i) = 1/ Ax3(k,j,i) * (   Ax3(k+1,j,i)*Vs(BX3s,k+1,j,i)  +   ( Ax1(k,j,i+1) * Vs(BX1s,k,j,i+1) - Ax1(k,j,i) * Vs(BX1s,k,j,i)
                                                                                                        +  Ax2(k,j+1,i) * Vs(BX2s,k,j+1,i) - Ax2(k,j,i) * Vs(BX2s,k,j,i) )) ;
                            }
                            for(int k = nend ; k<nx3 ; k++ ) {
                                Vs(BX3s,k+1,j,i) = 1/ Ax3(k+1,j,i) * (   Ax3(k,j,i)*Vs(BX3s,k,j,i)  -   ( Ax1(k,j,i+1) * Vs(BX1s,k,j,i+1) - Ax1(k,j,i) * Vs(BX1s,k,j,i)
                                                                                                        +  Ax2(k,j+1,i) * Vs(BX2s,k,j+1,i) - Ax2(k,j,i) * Vs(BX2s,k,j,i) )) ;
                            }
                            

                        });
    }

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
            case userdef:
                if(this->haveUserDefBoundary) this->userDefBoundaryFunc(data, dir, left, t);
                else IDEFIX_ERROR("No function has been enrolled to define your own boundary conditions");
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
            case userdef:
                if(this->haveUserDefBoundary) this->userDefBoundaryFunc(data, dir, right, t);
                else IDEFIX_ERROR("No function has been enrolled to define your own boundary conditions");
                break;
            default:
                std::stringstream msg("Boundary condition type is not yet implemented");
                IDEFIX_ERROR(msg);

        }
        #if MHD == YES
        // Reconstruct the normal field component when using CT
        ReconstructNormalField(data,dir);
        #endif

    }   // Loop on dimension ends

    #if MHD == YES
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
    IdefixArray3D<real> Ax1 = data.A[IDIR];
    IdefixArray3D<real> Ax2 = data.A[JDIR];
    IdefixArray3D<real> Ax3 = data.A[KDIR];
    IdefixArray3D<real> dV = data.dV;
    

    Kokkos::parallel_reduce("CheckDivB",
                                Kokkos::MDRangePolicy<Kokkos::Rank<3, Kokkos::Iterate::Right, Kokkos::Iterate::Right>>
                                ({data.beg[KDIR],data.beg[JDIR],data.beg[IDIR]},{data.end[KDIR], data.end[JDIR], data.end[IDIR]}),
                                KOKKOS_LAMBDA (int k, int j, int i, real &divBmax) {
                real dB1,dB2,dB3;

                dB1=dB2=dB3=ZERO_F;

                D_EXPAND( dB1=(Ax1(k,j,i+1)*Vs(BX1s,k,j,i+1)-Ax1(k,j,i)*Vs(BX1s,k,j,i)); ,
                          dB2=(Ax2(k,j+1,i)*Vs(BX2s,k,j+1,i)-Ax2(k,j,i)*Vs(BX2s,k,j,i)); ,
                          dB3=(Ax3(k+1,j,i)*Vs(BX3s,k+1,j,i)-Ax3(k,j,i)*Vs(BX3s,k,j,i));  )
                
                divBmax=FMAX(FABS(D_EXPAND(dB1, +dB2, +dB3))/dV(k,j,i),divBmax);

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
