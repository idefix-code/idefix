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

namespace {
constexpr real kPi = 3.14159265358979323846;
}  // namespace

// Serial-only tests: under WITH_MPI, FFT::R2C/C2R route through the distributed
// (MPI) path and expect local-sized arrays — that path is covered by
// FFT.MPI_Distributed_R2C_C2R below.
#ifndef WITH_MPI

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

#ifdef WITH_MPI

// MPI: Transpose<real> — same logic as Transpose<T>::Test() but with gtest assertions.
// in(i,j,k) = k + 10*j + 100*(i + rank*ni)  →  out(j,i,k) = k + 10*(j + rank*nj) + 100*i
TEST(FFT, TransposeReal_MPI) {
  const int p = idfx::psize;
  const int rank = idfx::prank;
  const int ni = 3;
  const int nj = 4;
  const int nk = 9;

  Transpose<real> transpose({ni, nj*p, nk});
  IdefixArray3D<real> in("in", ni, nj*p, nk);
  IdefixArray3D<real> out("out", nj, ni*p, nk);
  IdefixArray3D<real> expected("expected", nj, ni*p, nk);

  idefix_for("fill in", 0, ni, 0, nj*p, 0, nk,
    KOKKOS_LAMBDA(const int i, const int j, const int k) {
      in(i,j,k) = k + 10*j + 100*(i + rank*ni);
    });
  idefix_for("fill expected", 0, nj, 0, ni*p, 0, nk,
    KOKKOS_LAMBDA(const int j, const int i, const int k) {
      expected(j,i,k) = k + 10*(j + rank*nj) + 100*i;
    });

  transpose.Apply(in, out);

  ASSERT_THAT_3D_ARRAYS(out, testing::Eq, expected);
}

// MPI: Transpose<complex> — same pattern, with both real and imaginary parts carrying
// the test signal. The per-element error is computed on device as a real array, then
// checked with ASSERT_THAT_3D_ARRAY.
TEST(FFT, TransposeComplex_MPI) {
  const int p = idfx::psize;
  const int rank = idfx::prank;
  const int ni = 3;
  const int nj = 4;
  const int nk = 9;

  Transpose<complex> transpose({ni, nj*p, nk});
  IdefixArray3D<complex> in("in", ni, nj*p, nk);
  IdefixArray3D<complex> out("out", nj, ni*p, nk);
  IdefixArray3D<real> err("err", nj, ni*p, nk);

  idefix_for("fill in", 0, ni, 0, nj*p, 0, nk,
    KOKKOS_LAMBDA(const int i, const int j, const int k) {
      in(i,j,k) = complex(k + 10*j + 100*(i + rank*ni),
                         2*k + 20*j + 200*(i + rank*ni));
    });

  transpose.Apply(in, out);

  idefix_for("check", 0, nj, 0, ni*p, 0, nk,
    KOKKOS_LAMBDA(const int j, const int i, const int k) {
      complex expected(k + 10*(j + rank*nj) + 100*i,
                       2*k + 20*(j + rank*nj) + 200*i);
      complex d = out(j,i,k) - expected;
      err(j,i,k) = std::sqrt(d.real()*d.real() + d.imag()*d.imag());
    });

  ASSERT_THAT_3D_ARRAY(err, testing::Le, 1e-12);
}

// MPI: distributed R2C/C2R must agree with the serial FFT on the full global array.
// Reimplements FFT::TestMPI() with gtest assertions: every rank fills the same global
// array deterministically, extracts its local slice, runs R2C_MPI/C2R_MPI, and the
// per-element error against the serial reference is checked on device.
TEST(FFT, MPI_Distributed_R2C_C2R) {
  const int n = 16;
  const int p = idfx::psize;
  std::array<int,3> nr_glob = {n*p, n, n};
  std::array<int,3> nf_glob = {n*p, n, n/2+1};

  FFT fft(nr_glob, nf_glob);

  // Global arrays — filled identically on all ranks (deterministic, no scatter needed)
  IdefixArray3D<real> globalReal("global real", nr_glob[0], nr_glob[1], nr_glob[2]);
  IdefixArray3D<complex> globalComplex("global complex", nf_glob[0], nf_glob[1], nf_glob[2]);

  // hoist into a plain int: std::array::operator[] is not callable from a device lambda on nvcc
  const int nx = nr_glob[0];
  idefix_for("fill global", 0, nx, 0, nr_glob[1], 0, nr_glob[2],
    KOKKOS_LAMBDA(const int i, const int j, const int k) {
      globalReal(i,j,k) = std::cos(2.0*kPi*static_cast<real>(i)/static_cast<real>(nx))
                        + 0.5*std::sin(2.0*kPi*static_cast<real>(j)/static_cast<real>(n))
                        + 0.25*static_cast<real>(k);
    });

  // Local arrays (shape [n, n, n] / [n, n, n/2+1] on every rank)
  IdefixArray3D<real> localReal("local real", n, n, n);
  IdefixArray3D<complex> localComplex("local complex", n, n, n/2+1);

  // Extract the local slice of the global array
  const int offset = n * idfx::prank;
  idefix_for("extract local", 0, n, 0, n, 0, n,
    KOKKOS_LAMBDA(const int i, const int j, const int k) {
      localReal(i,j,k) = globalReal(i + offset, j, k);
    });

  // Serial reference FFT on the full global array
  KokkosFFT::rfftn(Kokkos::DefaultExecutionSpace(), globalReal, globalComplex);

  // Parallel R2C (transpose=true → output shape [n, n, n/2+1])
  fft.R2C_MPI(localReal, localComplex, true);

  // Per-element error vs serial reference, computed on device as a real array
  IdefixArray3D<real> errC("err C", n, n, n/2+1);
  idefix_for("check R2C", 0, n, 0, n, 0, n/2+1,
    KOKKOS_LAMBDA(const int i, const int j, const int k) {
      complex d = localComplex(i,j,k) - globalComplex(i + offset, j, k);
      errC(i,j,k) = std::sqrt(d.real()*d.real() + d.imag()*d.imag());
    });
  ASSERT_THAT_3D_ARRAY(errC, testing::Le, 1e-8);

  // Parallel C2R back to real
  fft.C2R_MPI(localComplex, localReal, true);

  // Per-element error vs the original global array
  IdefixArray3D<real> errR("err R", n, n, n);
  idefix_for("check C2R", 0, n, 0, n, 0, n,
    KOKKOS_LAMBDA(const int i, const int j, const int k) {
      errR(i,j,k) = std::abs(localReal(i,j,k) - globalReal(i + offset, j, k));
    });
  ASSERT_THAT_3D_ARRAY(errR, testing::Le, 1e-8);
}

#endif  // WITH_MPI
