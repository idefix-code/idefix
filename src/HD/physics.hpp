#ifndef PHYSICS_HPP
#define PHYSICS_HPP
#include "../idefix.hpp"

class Physics {
public:
    Physics();
    Physics(DataBlock &);
    void ConvertConsToPrim(DataBlock &);
    void ConvertPrimToCons(DataBlock &);
    void ExtrapolatePrimVar(DataBlock &, int);
    void CalcRiemannFlux(DataBlock &, IdefixArray3D<real> &, int);
    void CalcRightHandSide(DataBlock &, int, real );
    void InitFlow(DataBlock &);
    void SetBoundary(DataBlock &);

private:
    IdefixArray4D<real> PrimL;     // Left State (primitive variables)
    IdefixArray4D<real> PrimR;     // Right State (primitive variables)
    IdefixArray4D<real> FluxRiemann;      // Riemann flux at each interface

    real C2Iso;
    real gamma;

};



#endif