#include "../idefix.hpp"
#include "hydro.hpp"

#if MHD == YES
#include "solversMHD.hpp"
#else
#include "solversHD.hpp"
#endif

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
