#ifndef SETUP_HPP
#define SETUP_HPP
#include "idefix.hpp"

class Setup {
public:
    Setup();
    Setup(Input &, Grid &, DataBlock &);

    void InitFlow(DataBlock &);
    void MakeAnalysis(DataBlock&, real);
    void SetUserdefBoundary(DataBlock&, int, BoundarySide, real );
    void ComputeUserStep(DataBlock&, real, real);

};



#endif