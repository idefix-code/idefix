#ifndef TIMEINTEGRATOR_HPP
#define TIMEINTEGRATOR_HPP
#include "idefix.hpp"


class TimeIntegrator {
public:
    real getDt();      // Get current times step of time integrator
    void setDt(real );      // Set specific dt for time integrator  
    real getT();       // Get current time of time integrator


    // Constructor from input and given datablock
    TimeIntegrator(Input &, DataBlock &);

    void Stage();
    // Do one integration cycle
    void Cycle();

private:
    int nstages;
    // Weights of time integrator
    real w0[2];
    real wc[2];


    real dt;    // Current timestep
    real t;     // Current time
    
    IdefixArray4D<real> V0;     // Temporary Flow structure for multi-step integrators
    IdefixArray3D<real> InvDtHyp;  // Inverse Dt at each point for hyperbolic terms
    IdefixArray3D<real> InvDtPar;   // Inverse Dt at each point for parabolic terms
    DataBlock data;
};

#endif