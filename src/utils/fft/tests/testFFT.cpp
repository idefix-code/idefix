// ***********************************************************************************
// Idefix MHD astrophysical code
//
// Source file src/utils/fft/tests/testFFT.cpp
//
// Last modified : 11/2026
//
// Copyright(C) by :
// - Geoffroy Lesur <geoffroy.lesur@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2021 - 2026)
// and other code contributors
//
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

// headers
#include <array>
#include <cmath>
#include "gtest-idfx.hpp"
#include "fft.hpp"
#include "transpose.hpp"



// Serial-only tests: under WITH_MPI, FFT::R2C/C2R route through the distributed
// (MPI) path and expect local-sized arrays — that path is covered by the tests
// in testFFTMPI.cpp (unittests_mpi executable).
#ifndef WITH_MPI

namespace {
constexpr real kPi = M_PI;
}  // namespace

// R2C followed by C2R must recover the original array (roundtrip identity).
TEST(FFT, R2C_C2R_RoundTrip) {
  const int n = 16;
  std::array<int,3> nr = {n, n, n};
  std::array<int,3> nf = {n, n, n/2+1};

  FFT fft(nr, nf);

  IdefixArray3D<real> in("in", n, n, n);
  IdefixArray3D<complex> freq("freq", n, n, n/2+1);
  IdefixArray3D<real> out("out", n, n, n);

  // Fill on host and copy to device
  auto inHost = Kokkos::create_mirror_view(in);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      for (int k = 0; k < n; k++)
        inHost(i,j,k) = std::cos(2.0*kPi*i/n) * std::sin(2.0*kPi*j/n) + 0.5*k;
  Kokkos::deep_copy(in, inHost);

  fft.R2C(in, freq);
  fft.C2R(freq, out);

  ASSERT_THAT_3D_ARRAYS_ARGS(in, testing::DoubleNear, out, 1e-10);
}

// The R2C of cos(2*pi*i/n) (constant in j,k) must peak at i_freq = 1 and
// i_freq = n-1 with amplitude n^3/2, and vanish everywhere else.
TEST(FFT, R2C_AnalyticalCosine) {
  const int n = 16;
  std::array<int,3> nr = {n, n, n};
  std::array<int,3> nf = {n, n, n/2+1};

  FFT fft(nr, nf);

  IdefixArray3D<real> in("in", n, n, n);
  IdefixArray3D<complex> freq("freq", n, n, n/2+1);

  auto inHost = Kokkos::create_mirror_view(in);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      for (int k = 0; k < n; k++)
        inHost(i,j,k) = std::cos(2.0*kPi*i/n);
  Kokkos::deep_copy(in, inHost);

  fft.R2C(in, freq);

  auto freqHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), freq);

  const real expected = static_cast<real>(n*n*n)/2.0;
  // Peak at (1, 0, 0)
  ASSERT_NEAR(freqHost(1,0,0).real(), expected, 1e-8);
  ASSERT_NEAR(freqHost(1,0,0).imag(), 0.0, 1e-8);
  // Peak at (n-1, 0, 0)
  ASSERT_NEAR(freqHost(n-1,0,0).real(), expected, 1e-8);
  ASSERT_NEAR(freqHost(n-1,0,0).imag(), 0.0, 1e-8);
  // A few elements that must be zero
  ASSERT_NEAR(freqHost(0,0,0).real(), 0.0, 1e-8);
  ASSERT_NEAR(freqHost(2,0,0).real(), 0.0, 1e-8);
  ASSERT_NEAR(freqHost(1,1,0).real(), 0.0, 1e-8);
  ASSERT_NEAR(freqHost(1,0,1).real(), 0.0, 1e-8);
}

#endif  // WITH_MPI

// TransposeLocal swaps the last two indices: out(i,k,j) = in(i,j,k).
// Works identically with or without MPI (purely local index permutation).
TEST(FFT, TransposeLocal) {
  const int n = 8;
  IdefixArray3D<real> in("in", n, n, n);
  IdefixArray3D<real> out("out", n, n, n);
  IdefixArray3D<real> expected("expected", n, n, n);

  // in(i,j,k) = i + 100*j + 10000*k
  // out(i,j,k) = in(i,k,j) = i + 100*k + 10000*j
  idefix_for("fill", 0, n, 0, n, 0, n, KOKKOS_LAMBDA(const int i, const int j, const int k) {
    in(i,j,k) = i + 100*j + 10000*k;
    expected(i,j,k) = i + 100*k + 10000*j;
  });

  FFT fft;
  fft.TransposeLocal(in, out);

  ASSERT_THAT_3D_ARRAYS(out, testing::Eq, expected);
}
