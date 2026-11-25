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

// Self-managed xorshift64 random draws (see agent/adr/0002).
// One 64-bit state per (mode, component) is stored in rngState and advanced
// with the exact same step as Kokkos::Random_XorShift64::urand64(), and the
// normal draws use the same Box–Muller transform as
// Kokkos::Random_XorShift64::normal(). The distribution is therefore unchanged,
// while the sequence becomes deterministic on every backend and thread count
// and identical on all MPI ranks (unlike Kokkos::Random_XorShift64_Pool, whose
// no-argument get_state() selects a state through a race on GPU).

/// Advance `state` with the xorshift64 step and return a double in [0,1),
/// exactly like Kokkos::Random_XorShift64::drand() (urand64()/MAX_URAND64)
KOKKOS_INLINE_FUNCTION
double XorShift64Drand(uint64_t &state) {
  state ^= state >> 12;
  state ^= state << 25;
  state ^= state >> 27;
  return static_cast<double>((state * 2685821657736338717ULL) - 1)
         / static_cast<double>(std::numeric_limits<uint64_t>::max());
}

/// Standard normal draw through the Box–Muller transform, exactly like
/// Kokkos::Random_XorShift64::normal()
KOKKOS_INLINE_FUNCTION
double XorShift64Normal(uint64_t &state) {
  constexpr double two_pi = 2.0 * Kokkos::numbers::pi_v<double>;
  const double u = XorShift64Drand(state);
  const double v = XorShift64Drand(state);
  const double r = Kokkos::sqrt(-2.0 * Kokkos::log(u));
  const double theta = v * two_pi;
  return r * Kokkos::cos(theta);
}

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
  // dump as "ouValues" (complex) and "ouRngState" (uint64_t).
  IdefixArray2D<Kokkos::complex<real>> ouValues;
  IdefixArray2D<uint64_t> rngState;
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
