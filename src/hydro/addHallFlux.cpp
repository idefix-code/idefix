#include "../idefix.hpp"

void Hydro::AddHallFlux(DataBlock &data, int dir, const real t) {

    idfx::pushRegion("Hydro::AddHallFlux");
    
    IdefixArray4D<real> PrimL = data.PrimL;
    IdefixArray4D<real> PrimR = data.PrimR;
    IdefixArray4D<real> Flux = data.FluxRiemann;
    IdefixArray4D<real> J = data.J;
    IdefixArray3D<real> cMax = data.cMax;
    IdefixArray3D<real> xHallArr = data.xHall;
    IdefixArray1D<real> dx = data.dx[dir];
    IdefixArray1D<real> x1 = data.x[IDIR];
    IdefixArray1D<real> rt = data.rt;
    IdefixArray1D<real> dmu = data.dmu;
    ParabolicType haveHall = this->haveHall;
    
#if EMF_AVERAGE != ARITHMETIC
    IDEFIX_ERROR("the Hall effect module is demonstrated stable only when using EMF_AVERAGE=ARITHMETIC");
#endif

#if HAVE_ENERGY
    real gamma_m1 = this->gamma-1;
#endif
    // References to required emf components
    IdefixArray3D<real> Eb;
    IdefixArray3D<real> Et;

    // Only used when using constant Hall effect diffusivity
    real xHConstant = this->xH;

    int ioffset,joffset,koffset;
    int iextend, jextend,kextend;

    ioffset=joffset=koffset=0;
    // extension in perp to the direction of integration, as required by CT.
    iextend=jextend=kextend=0;


    if(haveHall == UserDefFunction && dir == IDIR) {
        if(hallDiffusivityFunc) hallDiffusivityFunc(data, t, xHallArr);
        else IDEFIX_ERROR("No user-defined Hall diffusivity function has been enrolled");
    }


    // Define normal, tangent and bi-tanget indices
    int nDir, tDir, bDir;
    real st,sb;



    switch(dir) {
        case(IDIR):
            ioffset = 1;
            D_EXPAND(               ,
                   jextend = 1;     ,
                   kextend = 1; )

            nDir = IDIR;
            tDir = JDIR;
            bDir = KDIR;

            Et = data.emf.ezi;
            Eb = data.emf.eyi;

            D_EXPAND( st = -1.0;  ,
                                  ,
                      sb = +1.0;  )
            break;
        case(JDIR):
            joffset=1;
            D_EXPAND( iextend = 1;  ,
                                    ,
                    kextend = 1;)

            nDir = JDIR;
            tDir = IDIR;
            bDir = KDIR;

            Et = data.emf.ezj;
            Eb = data.emf.exj;

            D_EXPAND( st = +1.0;  ,
                                  ,
                      sb = -1.0;  )
            break;
        case(KDIR):
            koffset=1;
            D_EXPAND( iextend = 1;  ,
                      jextend = 1;    ,
                    )
            
            nDir = KDIR;
            tDir = IDIR;
            bDir = JDIR;

            Et = data.emf.eyk;
            Eb = data.emf.exk;

            st = -1.0;  
            sb = +1.0;

            break;
        default:
            IDEFIX_ERROR("Wrong direction");
    }


    idefix_for("CalcHallFlux",data.beg[KDIR]-kextend,data.end[KDIR]+koffset+kextend,data.beg[JDIR]-jextend,data.end[JDIR]+joffset+jextend,data.beg[IDIR]-iextend,data.end[IDIR]+ioffset+iextend,
                    KOKKOS_LAMBDA (int k, int j, int i) 
        {
            int ip1, jp1, kp1;  // Offset indieces

            real Jx1, Jx2, Jx3;
            real Bx1, Bx2, Bx3;

            real vH[3];
            real Eh[3];
            real FluxL_BXt, FluxL_BXb;
            real FluxR_BXt, FluxR_BXb;

            real Bmag2L,Bmag2R;

            #if HAVE_ENERGY
                real FluxL_ENG, FluxR_ENG;
                real vHB;
            #endif

            real xH;

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

           
            xH = xHConstant;

             // Compute current on the right face
            if(dir == IDIR) {
                Jx1 = AVERAGE_4D_XYZ(J, IDIR, kp1, jp1, i);
                Jx2 = AVERAGE_4D_Z(J, JDIR, kp1, j, i);
                Jx3 = AVERAGE_4D_Y(J, KDIR, k, jp1, i);   

                if(haveHall == UserDefFunction) xH = AVERAGE_3D_X(xHallArr,k,j,i);
            }
            if(dir == JDIR) {
                Jx1 = AVERAGE_4D_Z(J, IDIR, kp1, j, i);
                Jx2 = AVERAGE_4D_XYZ(J, JDIR, kp1, j, ip1);
                Jx3 = AVERAGE_4D_X(J, KDIR, k, j, ip1);     

                if(haveHall == UserDefFunction) xH = AVERAGE_3D_Y(xHallArr,k,j,i);             
            }
            if(dir == KDIR) {
                Jx1 = AVERAGE_4D_Y(J, IDIR, k, jp1, i);
                Jx2 = AVERAGE_4D_X(J, JDIR, k, j, ip1);
                Jx3 = AVERAGE_4D_XYZ(J, KDIR, k, jp1, ip1);

                if(haveHall == UserDefFunction) xH = AVERAGE_3D_Z(xHallArr,k,j,i);  
            }
            
            // Compute drift speed
            vH[IDIR] = -Jx1*xH;
            vH[JDIR] = -Jx2*xH;
            vH[KDIR] = -Jx3*xH;

            // Up to now, everything is assumed to be homogeneous on the cell face.
            // We now compute the fluxes on both sides of the cells, and use an HLL solver to get the riemann flux
            
            // Left Side
            EXPAND(                 ,
                    FluxL_BXt = vH[nDir] * PrimL(tDir+BX1,k,j,i) - PrimL(nDir+BX1,k,j,i) * vH[tDir]; ,
                    FluxL_BXb = vH[nDir] * PrimL(bDir+BX1,k,j,i) - PrimL(nDir+BX1,k,j,i) * vH[bDir]; )

            Bmag2L = EXPAND(PrimL(BX1,k,j,i) * PrimL(BX1,k,j,i), 
                         + PrimL(BX2,k,j,i) * PrimL(BX2,k,j,i), 
                         + PrimL(BX3,k,j,i) * PrimL(BX3,k,j,i));

            #if HAVE_ENERGY
                    vHB = EXPAND( vH[IDIR] * PrimL(BX1,k,j,i)    ,
                                + vH[JDIR] * PrimL(BX2,k,j,i)    ,
                                + vH[KDIR] * PrimL(BX3,k,j,i) );
                    
                    FluxL_ENG = vH[bDir] * Bmag2L - vHB * PrimL(nDir+BX1,k,j,i);
            #endif

            // Right Side
            EXPAND(                 ,
                    FluxR_BXt = vH[nDir] * PrimR(tDir+BX1,k,j,i) - PrimR(nDir+BX1,k,j,i) * vH[tDir]; ,
                    FluxR_BXb = vH[nDir] * PrimR(bDir+BX1,k,j,i) - PrimR(nDir+BX1,k,j,i) * vH[bDir]; )

            Bmag2R = EXPAND(PrimR(BX1,k,j,i) * PrimR(BX1,k,j,i), 
                         + PrimR(BX2,k,j,i) * PrimR(BX2,k,j,i), 
                         + PrimR(BX3,k,j,i) * PrimR(BX3,k,j,i));  

            #if HAVE_ENERGY
                    vHB = EXPAND( vH[IDIR] * PrimR(BX1,k,j,i)    ,
                                + vH[JDIR] * PrimR(BX2,k,j,i)    ,
                                + vH[KDIR] * PrimR(BX3,k,j,i) );
                    
                    FluxR_ENG = vH[bDir] * Bmag2R - vHB * PrimR(nDir+BX1,k,j,i);
            #endif

            // Compute whistler wave speed
            // Get the physical grid lengthcscale
            const int ig = ioffset*i + joffset*j + koffset*k;
            real dl = dx(ig);
            #if GEOMETRY == POLAR
                if(dir==JDIR) dl = dl*x1(i);
            #elif GEOMETRY == SPHERICAL
                if(dir==JDIR) dl = dl*rt(i);
                if(dir==KDIR) dl = dl*rt(i)*dmu(j)/dx2(j);
            #endif

            real cwR = FABS(xH) * sqrt(Bmag2R) / (dl);
            real cwL = FABS(xH) * sqrt(Bmag2L) / (dl);            

            real SR = FMAX(cwR,cwL);   

            real SL = -SR;
            real Sdiff = SR-SL;
            if(Sdiff < 1e-10) Sdiff = 1e-10;

            // Compute left/right Conservative variable used for the HLL flux
            EXPAND(                                         ,
                    real uR_BXt = PrimR(BX1+tDir,k,j,i);
                    real uL_BXt = PrimL(BX1+tDir,k,j,i);    ,
                    real uR_BXb = PrimR(BX1+bDir,k,j,i);
                    real uL_BXb = PrimL(BX1+bDir,k,j,i);  )

            #if HAVE_ENERGY    
                real uL_ENG = PrimL(PRS,k,j,i) / gamma_m1 
                    + HALF_F * PrimL(RHO,k,j,i) * (EXPAND(  PrimL(VX1,k,j,i) * PrimL(VX1,k,j,i) ,
                                                          + PrimL(VX2,k,j,i) * PrimL(VX2,k,j,i) ,
                                                          + PrimL(VX3,k,j,i) * PrimL(VX3,k,j,i) ))
                    + HALF_F * (EXPAND(  PrimL(BX1,k,j,i) * PrimL(BX1,k,j,i) ,
                                       + PrimL(BX2,k,j,i) * PrimL(BX2,k,j,i) ,
                                       + PrimL(BX3,k,j,i) * PrimL(BX3,k,j,i) ));
                
                real uR_ENG = PrimR(PRS,k,j,i) / gamma_m1 
                    + HALF_F * PrimR(RHO,k,j,i) * (EXPAND(  PrimR(VX1,k,j,i) * PrimR(VX1,k,j,i) ,
                                                          + PrimR(VX2,k,j,i) * PrimR(VX2,k,j,i) ,
                                                          + PrimR(VX3,k,j,i) * PrimR(VX3,k,j,i) ))
                    + HALF_F * (EXPAND(  PrimR(BX1,k,j,i) * PrimR(BX1,k,j,i) ,
                                       + PrimR(BX2,k,j,i) * PrimR(BX2,k,j,i) ,
                                       + PrimR(BX3,k,j,i) * PrimR(BX3,k,j,i) ));
            #endif

            // Compute HLL Flux and add it to Riemann flux, knowing in advance we're in the star zone
            EXPAND(                                                                                                   ,
                     real Fluxt = ( SL*SR*(uR_BXt - uL_BXt) + SR * FluxL_BXt - SL * FluxR_BXt ) /Sdiff;
                     Flux(BX1+tDir,k,j,i) += Fluxt;                                                                   ,
                     real Fluxb = ( SL*SR*(uR_BXb - uL_BXb) + SR * FluxL_BXb - SL * FluxR_BXb ) /Sdiff;
                     Flux(BX1+bDir,k,j,i) +=  Fluxb;  )
            
            #if HAVE_ENERGY
                Flux(ENG,k,j,i) += (SL*SR*(uR_ENG - uL_ENG) + SR * FluxL_ENG - SL * FluxR_ENG) / Sdiff;
            #endif

            // Add whistler wave speed to the max signal speed
            cMax(k,j,i) += TWO_F*TWO_F*SR; // TWO beacuse the CFL for Hall is like a second order diffusion operator, for which there's a factor 2

            // Store the new flux in the emf components
            D_EXPAND(Et(k,j,i) += st*Fluxt; ,
                                                    ,
                    Eb(k,j,i) += sb*Fluxb; )

        });


    idfx::popRegion();

}
