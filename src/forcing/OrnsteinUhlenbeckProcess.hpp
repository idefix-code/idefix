// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef FORCING_ORNSTEINUHLENBECKPROCESS_HPP_
#define FORCING_ORNSTEINUHLENBECKPROCESS_HPP_

#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "idefix.hpp"
#include "random.hpp"

class OrnsteinUhlenbeckProcesses {
 private:
  std::vector<std::vector<std::string>> names;
  int nSeries;

 public:
  // OU parameters per (mode, component). They are replicated on every rank
  // and stored in the restart dump (see agent/adr/0002).
  IdefixArray2D<real> means;
  IdefixArray2D<real> tcorrs;
  IdefixArray2D<real> epsilons;

  // Complete OU state, replicated on every rank and stored in the restart
  // dump as "ouValues" (complex) and "ouRngState" (uint64_t, via rng.state()).
  IdefixArray2D<Kokkos::complex<real>> ouValues;
  idfx::RandomXorShift64 rng;
  IdefixHostArray2D<Kokkos::complex<real>> ouValuesHost;
  IdefixArray2D<real> normalValuesReal;
  IdefixArray2D<real> normalValuesImag;
  IdefixHostArray2D<real> normalValuesRealHost;
  IdefixHostArray2D<real> normalValuesImagHost;

  OrnsteinUhlenbeckProcesses();  // Default (empty) constructor
  void InitProcesses(std::string, int, int, std::vector<std::vector<std::string>>);
  void SetProcesses(IdefixArray2D<real>, IdefixArray2D<real>, IdefixArray2D<real>);
  void UpdateProcessesValues(real);

  std::string ouFilename;
  std::string normalFilename;
  void ResetProcessesValues();
  void WriteProcessesValues(real);
  void ResetNormalValues();
  void WriteNormalValues(real);
  int precision;

  std::ofstream file;
};

#endif  // FORCING_ORNSTEINUHLENBECKPROCESS_HPP_
