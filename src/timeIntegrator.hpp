// ***********************************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************************

#ifndef TIMEINTEGRATOR_HPP_
#define TIMEINTEGRATOR_HPP_

#include "idefix.hpp"



class TimeIntegrator {
  friend class OutputDump;

 public:
  int64_t getNcycles();   // Get current number of cycles

  // Constructor from input and given datablock
  TimeIntegrator(Input &, DataBlock &);

  // Do one integration cycle
  void Cycle(DataBlock &);


 private:
  int nstages;
  // Weights of time integrator
  real w0[2];
  real wc[2];

  real cfl;   // CFL number
  real cflMaxVar; // Max CFL variation number
  int64_t ncycles;        // # of cycles
  double lastLog;         // # time for the last log
  double lastMpiLog;      // # time for the last MPI log
  int64_t cyclePeriod;    // # of cycles between two logs
  Kokkos::Timer timer;    // Internal timer of the integrator
};

#endif // TIMEINTEGRATOR_HPP_
