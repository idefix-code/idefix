// ********************************************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ********************************************************************************************************

#ifndef TIMEINTEGRATOR_HPP
#define TIMEINTEGRATOR_HPP
#include "idefix.hpp"



class TimeIntegrator {
    friend class OutputDump;

public:
    real getDt();      // Get current times step of time integrator
    int64_t getNcycles();  // Get current number of cycles
    void setDt(real );      // Set specific dt for time integrator  
    real getT();       // Get current time of time integrator


    // Constructor from input and given datablock
    TimeIntegrator(Input &, Hydro &);

    void Stage(DataBlock &);
    // Do one integration cycle
    void Cycle(DataBlock &);
    void ReinitInvDt(DataBlock & );

    // Return the hydro object used by present TimeIntegrator
    Hydro &GetHydro();

    

private:
    int nstages;
    // Weights of time integrator
    real w0[2];
    real wc[2];


    real dt;    // Current timestep
    real t;     // Current time
    real cfl;   // CFL number
    int64_t ncycles;   // # of cycles
    double lastLog;     // # time for the last log
    double lastMpiLog;  // # time for the last MPI log
    long int cyclePeriod;   // # of cycles between two logs
    Kokkos::Timer timer;    // Internal timer of the integrator
    
    Hydro *hydro;
};

#endif