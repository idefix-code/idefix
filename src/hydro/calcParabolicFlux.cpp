#include "../idefix.hpp"
#include "hydro.hpp"

// Compute parabolic fluxes
void Hydro::CalcParabolicFlux(DataBlock &data, int dir) {
    idfx::pushRegion("Hydro::CalcParabolicFlux");

    int ioffset,joffset,koffset;

    IdefixArray4D<real> Flux = data.FluxRiemann;
    IdefixArray4D<real> Vc = data.Vc;
    IdefixArray4D<real> Vs = data.Vs;
    IdefixArray3D<real> dMax = data.dMax;
    IdefixArray4D<real> J = data.J;
    real etaConstant = this->etaO;

    ioffset=joffset=koffset=0;

    switch(dir) {
        case(IDIR):
            ioffset = 1;
            break;
        case(JDIR):
            joffset = 1;
            break;
        case(KDIR):
            koffset=1;
            break;
        default:
            IDEFIX_ERROR("Wrong direction");
    }

    // Note the flux follows the same sign convention as the hyperbolic flux
    // HEnce signs are reversed compared to the parabolic fluxes found in Pluto 4.3
    idefix_for("CalcParabolicFlux",data.beg[KDIR],data.end[KDIR]+koffset,data.beg[JDIR],data.end[JDIR]+joffset,data.beg[IDIR],data.end[IDIR]+ioffset,
                KOKKOS_LAMBDA (int k, int j, int i) 
            {
                real Jx1, Jx2, Jx3;
                real Bx1, Bx2, Bx3;

                int ip1, jp1, kp1;  // Offset indieces

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

                
                Jx1=Jx2=Jx3=ZERO_F;

                real eta = etaConstant;
                dMax(k,j,i) = ZERO_F;

                if(dir==IDIR) {
                    EXPAND(                                            ,
                            Jx3 = AVERAGE_4D_Y(J, KDIR, k, jp1, i);    ,
                            Jx2 = AVERAGE_4D_Z(J, JDIR, kp1, j, i);
                            Jx1 = AVERAGE_4D_XYZ(J, IDIR, kp1, jp1, i);  )

                    EXPAND( Bx1 = Vs(BX1s,k,j,i);                             ,
                            Bx2 = HALF_F*( Vc(BX2,k,j,i) + Vc(BX2,k,j,ip1)); ,
                            Bx3 = HALF_F*( Vc(BX3,k,j,i) + Vc(BX3,k,j,ip1)); )

                    EXPAND(                                           ,
                            Flux(BX2,k,j,i) += - eta * Jx3;               ,
                            Flux(BX3,k,j,i) +=   eta * Jx2; )
                    
                    #if HAVE_ENERGY
                            Flux(ENG,k,j,i) += EXPAND( ZERO_F               ,
                                                       - Bx2 * eta * Jx3    ,
                                                       + Bx3 * eta * Jx2);
                    #endif

                    dMax(k,j,i) += eta;

                }
                if(dir==JDIR) {
                    EXPAND( Jx3 = AVERAGE_4D_X(J, KDIR, k, j, ip1);   , 
                                                                      ,
                            Jx1 = AVERAGE_4D_Z(J, IDIR, kp1, j, i);
                            Jx2 = AVERAGE_4D_XYZ(J, JDIR, kp1, j, ip1);  )

                    EXPAND( Bx1 = HALF_F*( Vc(BX1,k,j,i) + Vc(BX1,k,jp1,i));    ,
                            Bx2 = Vs(BX2s,k,j,i);                               ,
                            Bx3 = HALF_F*( Vc(BX3,k,j,i) + Vc(BX3,k,jp1,i)); )

                    EXPAND( Flux(BX1,k,j,i) += eta * Jx3;   ,  
                                                            ,
                            Flux(BX3,k,j,i) += - eta * Jx1; )
                    
                    #if HAVE_ENERGY
                            Flux(ENG,k,j,i) += EXPAND(  Bx1 * eta * Jx3     , 
                                                                            ,
                                                       - Bx3 * eta * Jx1); 
                    #endif
                    dMax(k,j,i) += eta;

                }
                if(dir==KDIR) {
                    Jx1 = AVERAGE_4D_X(J, IDIR, k, j, ip1);
                    Jx2 = AVERAGE_4D_Y(J, JDIR, k, jp1, i);
                    Jx3 = AVERAGE_4D_XYZ(J, KDIR, k, jp1, ip1);

                    Bx1 = HALF_F*( Vc(BX1,k,j,i) + Vc(BX1,kp1,j,i)); 
                    Bx2 = HALF_F*( Vc(BX2,k,j,i) + Vc(BX2,kp1,j,i)); 
                    Bx3 = Vs(BX3s,k,j,i);

                    Flux(BX1,k,j,i) += -eta * Jx2;
                    Flux(BX2,k,j,i) += eta * Jx1;

                    
                    #if HAVE_ENERGY
                        Flux(ENG,k,j,i) += - Bx1 * eta * Jx2 + Bx2 * eta * Jx1;  
                    #endif
                    dMax(k,j,i) += eta;
                }

            });

}
