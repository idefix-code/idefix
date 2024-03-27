// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef TIMEINTEGRATOR_HPP_
#define TIMEINTEGRATOR_HPP_

#include "idefix.hpp"
#include "dataBlock.hpp"
#include "rkl.hpp"



class TimeIntegrator {
 public:
  int64_t GetNCycles();   // Get current number of cycles

  // Constructor from input and given datablock
  TimeIntegrator(Input &, DataBlock &);

  // Do one integration cycle
  void Cycle(DataBlock &);

  // check whether we have reached the maximum runtime
  bool CheckForMaxRuntime();

  void ShowLog(DataBlock &);    //<  Display progress log
  void ShowConfig();            //< Show configuration of time integrator

  bool isSilent{false};   // Whether the integration should proceed silently

 private:
  double ComputeBalance(); // Compute the compute balance between MPI processes

  // Whether we have RKL
  bool haveRKL{false};

  int nstages;
  // Weights of time integrator
  real w0[2];
  real wc[2];

  int checkNanPeriodicity{1};

  bool haveFixedDt = false;
  real fixedDt;

  real cfl;   // CFL number
  real cflMaxVar; // Max CFL variation number
  real maxdivB{0};   // Maximum allowed divB
  int64_t ncycles;        // # of cycles

  double computeLastLog;  // Timer for actual computeTime

  double lastLog;         // time for the last log (s)
  double lastMpiLog;      // time for the last MPI log (s)
  double lastSGLog;      // time for the last SelfGravity log (s)
  double maxRuntime;      // Maximum runtime requested (disabled when negative)
  int64_t cyclePeriod;    // # of cycles between two logs
  Kokkos::Timer timer;    // Internal timer of the integrator
};

#endif // TIMEINTEGRATOR_HPP_
