#ifndef PHYSICS_HPP
#define PHYSICS_HPP

#include "../idefix.hpp"

// Solver type
enum Solver {TVDLF=1, HLL, HLLC, ROE};

class Physics {
public:
    Physics();
    Physics(Input &, Setup &);
    void ConvertConsToPrim(DataBlock &);
    void ConvertPrimToCons(DataBlock &);
    void ExtrapolatePrimVar(DataBlock &, int);
    void CalcRiemannFlux(DataBlock &, int);
    void CalcRightHandSide(DataBlock &, int, real );
    void SetBoundary(DataBlock &, real);

private:

    real C2Iso;
    real gamma;

    Setup mySetup;
    Solver mySolver;

};


#endif
