#include "../idefix.hpp"
#include "physics.hpp"

#include "solvers.hpp"

Physics::Physics(Input &input, Setup &setup) {
    Kokkos::Profiling::pushRegion("Physics::Physics(DataBock)");

    this->gamma = 5.0/3.0;
    this->C2Iso = 1.0;

    this->mySetup=setup;
    
    // read Solver from input file
    std::string solverString = input.GetString("Solver","Solver",0);
    
    if (solverString.compare("tvdlf") == 0)     mySolver = TVDLF;
    else if (solverString.compare("hll") == 0)  mySolver = HLL;
    else if (solverString.compare("hllc") == 0) mySolver = HLLC;
    else if (solverString.compare("roe") == 0)  mySolver = ROE;
    else {
        std::stringstream msg;
        msg << "Unknown HD solver type " << solverString;
        IDEFIX_ERROR(msg);
    }

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

    Kokkos::Profiling::pushRegion("Physics::CalcRiemannFlux");

    switch (mySolver) {
        case TVDLF: Tvdlf(data, dir, this->gamma, this->C2Iso);
            break;
        case HLL:   Hll(data, dir, this->gamma, this->C2Iso);
            break;
        case HLLC:  Hllc(data, dir, this->gamma, this->C2Iso);
            break;
        case ROE:   Roe(data, dir, this->gamma, this->C2Iso);
            break;
        default: // do nothing
            break;
    }
    
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


