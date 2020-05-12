#ifndef HYDRO_HPP
#define HYDRO_HPP
#include "../idefix.hpp"



#define     SMALL_PRESSURE_FIX      (1.0e-5)

// Solver type
#if MHD == YES
enum Solver {TVDLF=1, HLL, HLLD, ROE};
#else
enum Solver {TVDLF=1, HLL, HLLC, ROE};
#endif

class Hydro {
public:
    Hydro();
    Hydro(Input &, Setup &);
    void ConvertConsToPrim(DataBlock &);
    void ConvertPrimToCons(DataBlock &);
    void ExtrapolatePrimVar(DataBlock &, int);
    void CalcRiemannFlux(DataBlock &, int);
    void CalcRightHandSide(DataBlock &, int, real );
    void ReconstructVcField(DataBlock &, IdefixArray4D<real> &);
    void ReconstructNormalField(DataBlock &);
    void EvolveMagField(DataBlock &, real, real);
    void CalcCornerEMF(DataBlock &, real );
    void SetBoundary(DataBlock &, real);
    void SetGamma(real);
    real CheckDivB(DataBlock &);

private:

    real C2Iso;
    real gamma;

    Setup mySetup;
    Solver mySolver;

};



#endif
