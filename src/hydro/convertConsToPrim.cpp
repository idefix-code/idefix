#include "../idefix.hpp"
#include "hydro.hpp"

#if MHD == YES
#include "solversMHD.hpp"
#else
#include "solversHD.hpp"
#endif

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