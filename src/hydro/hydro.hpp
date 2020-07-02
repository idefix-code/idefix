#ifndef HYDRO_HPP
#define HYDRO_HPP
#include "../idefix.hpp"



#define     SMALL_PRESSURE_FIX      (1.0e-5)
#define     eps_UCT_CONTACT          1.e-6

// Solver type
#if MHD == YES
enum Solver {TVDLF=1, HLL, HLLD, ROE};
#else
enum Solver {TVDLF=1, HLL, HLLC, ROE};
#endif

class Hydro {
public:
    Hydro();
    Hydro(Input &, Grid &);
    void ConvertConsToPrim(DataBlock &);
    void ConvertPrimToCons(DataBlock &);
    void ExtrapolatePrimVar(DataBlock &, int);
    void CalcRiemannFlux(DataBlock &, int);
    void CalcRightHandSide(DataBlock &, int, real );
    void AddSourceTerms(DataBlock &, real );
    void ReconstructVcField(DataBlock &, IdefixArray4D<real> &);
    void ReconstructNormalField(DataBlock &);
    void EvolveMagField(DataBlock &, real, real);
    void CalcCornerEMF(DataBlock &, real );
    void SetBoundary(DataBlock &, real);
    void SetGamma(real);
    real GetGamma();
    real GetC2iso();
    real CheckDivB(DataBlock &);

    // Source terms
    bool haveSourceTerms;
private:

    real C2Iso;
    real gamma;

    Solver mySolver;

    // Rotation vector
    bool haveRotation;
    real OmegaX1, OmegaX2, OmegaX3;

    bool haveShearingBox;
    // Shear rate for shearing box problems
    real sbS;
    // Box width for shearing box problems
    real sbLx;
};



#endif
