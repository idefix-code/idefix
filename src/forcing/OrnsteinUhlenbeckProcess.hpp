// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef FORCING_ORNSTEINUHLENBECKPROCESS_HPP_
#define FORCING_ORNSTEINUHLENBECKPROCESS_HPP_

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <Kokkos_Random.hpp>

#include "idefix.hpp"

class OrnsteinUhlenbeckProcesses {
 private:
  std::vector<std::vector<std::string>> names;
  int nSeries;

  IdefixArray2D<real> means;
  IdefixArray2D<real> tcorrs;
  IdefixArray2D<real> epsilons;

  Kokkos::Random_XorShift64_Pool<> random_pool;

 public:
  IdefixArray2D<Kokkos::complex<real>> ouValues;
  IdefixHostArray2D<Kokkos::complex<real>> ouValuesHost;
  IdefixArray2D<real> normalValuesReal;
  IdefixArray2D<real> normalValuesImag;
  IdefixHostArray2D<real> normalValuesRealHost;
  IdefixHostArray2D<real> normalValuesImagHost;

  OrnsteinUhlenbeckProcesses();  // Default (empty) constructor
  void InitProcesses(std::string, int, int, std::vector<std::vector<std::string>>);
  void SetProcesses(IdefixArray2D<real>, IdefixArray2D<real>, IdefixArray2D<real>);
  void UpdateProcessesValues(real);
  void AdvanceProcessesValues();

  std::string ouFilename;
  std::string normalFilename;
  std::string timestepFilename;
  void ResetProcessesValues();
  void WriteProcessesValues(real);
  void ResetNormalValues();
  void WriteNormalValues(real);
  void ResetTimestep();
  void WriteTimestep(real, real);
  int precision;

  std::ofstream file;
};

#endif  // FORCING_ORNSTEINUHLENBECKPROCESS_HPP_
