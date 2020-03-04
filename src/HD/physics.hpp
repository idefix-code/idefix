#ifndef PHYSICS_HPP
#define PHYSICS_HPP
#include "../idefix.hpp"

class Physics {
public:
    Physics();
    Physics(Input &);
    void ConvertConsToPrim(DataBlock &);
    void ConvertPrimToCons(DataBlock &);
    void ExtrapolatePrimVar(DataBlock &, int);
    void CalcRiemannFlux(DataBlock &, int);
    void CalcRightHandSide(DataBlock &, int, real );
    void InitFlow(DataBlock &);
    void SetBoundary(DataBlock &);

private:

    real C2Iso;
    real gamma;

    real randm(void);

};



#endif