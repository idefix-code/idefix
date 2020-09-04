#ifndef HYDRO_HPP
#define HYDRO_HPP
#include "../idefix.hpp"

#define     SMALL_PRESSURE_FIX      (1.0e-5)
#define     eps_UCT_CONTACT         (1.0e-6)

// Solver type
#if MHD == YES
enum Solver {TVDLF=1, HLL, HLLD, ROE};
#else
enum Solver {TVDLF=1, HLL, HLLC, ROE};
#endif

/*---- EMFs -----*/
#define ARITHMETIC   1
#define UCT0         2
#define UCT_CONTACT  3

#ifndef EMF_AVERAGE
    #define EMF_AVERAGE     UCT_CONTACT
#endif

// How are set parabolic terms
enum ParabolicType {Disabled, Constant, UserDefFunction };

using UserDefBoundaryFunc = void (*) (DataBlock &, int dir, BoundarySide side, const real t);
using GravPotentialFunc = void (*) (DataBlock &, const real t, IdefixArray1D<real>&, IdefixArray1D<real>&, IdefixArray1D<real>&, IdefixArray3D<real> &);
using SrcTermFunc = void (*) (DataBlock &, const real t, const real dt);

class Hydro {
public:
    Hydro();
    Hydro(Input &, Grid &);
    void ConvertConsToPrim(DataBlock &);
    void ConvertPrimToCons(DataBlock &);
    void ExtrapolatePrimVar(DataBlock &, int);
    void CalcRiemannFlux(DataBlock &, int);
    void CalcParabolicFlux(DataBlock &, int);
    void CalcRightHandSide(DataBlock &, int, real, real );
    void CalcCurrent(DataBlock &);
    void AddSourceTerms(DataBlock &, real, real );
    void ReconstructVcField(DataBlock &, IdefixArray4D<real> &);
    void ReconstructNormalField(DataBlock &, int);
    void EvolveMagField(DataBlock &, real, real);
    void CalcCornerEMF(DataBlock &, real );
    void CalcNonidealEMF(DataBlock &, real );
    void SetBoundary(DataBlock &, real);
    void SetGamma(real);
    real GetGamma();
    real GetC2iso();
    real CheckDivB(DataBlock &);

    // Source terms
    bool haveSourceTerms;

    // Parabolic terms
    bool haveParabolicTerms;
    
    // Current
    bool needCurrent;

    // Nonideal MHD effects coefficients
    ParabolicType haveResistivity, haveAmbipolar, haveHall;

    // Enroll user-defined boundary conditions
    void EnrollUserDefBoundary(UserDefBoundaryFunc);

    // Enroll user-defined gravitational potential
    void EnrollGravPotential(GravPotentialFunc);

    // Add some user source terms
    void EnrollUserSourceTerm(SrcTermFunc);

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

    // User defined Boundary conditions
    UserDefBoundaryFunc userDefBoundaryFunc;
    bool haveUserDefBoundary;

    // User defined gravitational potential
    GravPotentialFunc gravPotentialFunc;
    bool haveGravPotential;

    // User defined source term
    SrcTermFunc userSourceTerm;
    bool haveUserSourceTerm;

    real etaO, xH, xA;  // Ohmic resistivity, Hall, ambipolar

};



#endif
