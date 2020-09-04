#include "../idefix.hpp"
#include "hydro.hpp"


// Compute Corner EMFs from nonideal MHD
void Hydro::CalcNonidealEMF(DataBlock &data, real t) {
        idfx::pushRegion("Hydro::CalcNonidealEMF");

    // Corned EMFs
    IdefixArray3D<real> ex = data.emf.ex;
    IdefixArray3D<real> ey = data.emf.ey;
    IdefixArray3D<real> ez = data.emf.ez;
    IdefixArray4D<real> J = data.J;

    real etaConstant = this->etaO;


#if MHD == YES
    idefix_for("CalcNIEMF",
                data.beg[KDIR],data.end[KDIR]+KOFFSET,
                data.beg[JDIR],data.end[JDIR]+JOFFSET,
                data.beg[IDIR],data.end[IDIR]+IOFFSET,
                KOKKOS_LAMBDA (int k, int j, int i) {

                    // CT_EMF_ArithmeticAverage (emf, 0.25);
                    real eta = etaConstant;
    #if DIMENSIONS == 3
                    ex(k,j,i) += eta * J(IDIR, k, j, i);
                    ey(k,j,i) += eta * J(JDIR, k, j, i);
    #endif
                    ez(k,j,i) += eta * J(KDIR, k, j, i);

                });
#endif

    idfx::popRegion();
}