#ifndef TIMEINTEGRATOR_HPP
#define TIMEINTEGRATOR_HPP
#include "idefix.hpp"


class TimeIntegrator {
public:
    real getDt();      // Get current times step of time integrator
    long int getNcycles();  // Get current number of cycles
    void setDt(real );      // Set specific dt for time integrator  
    real getT();       // Get current time of time integrator


    // Constructor from input and given datablock
    TimeIntegrator(Input &, Hydro &, Setup &);

    void Stage(DataBlock &);
    // Do one integration cycle
    void Cycle(DataBlock &);
    void ReinitInvDt(DataBlock & );
private:
    int nstages;
    // Weights of time integrator
    real w0[2];
    real wc[2];


    real dt;    // Current timestep
    real t;     // Current time
    real cfl;   // CFL number
    long int ncycles;   // # of cycles
    double lastLog;     // # time for the last log
    Kokkos::Timer timer;    // Internal timer of the integrator
    
    DataBlock data;
    Hydro hydro;
    Setup mySetup;

    
};

#endif