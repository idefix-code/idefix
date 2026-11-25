// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include <cstdint>
#include <string>
#include <vector>

#include "OrnsteinUhlenbeckProcess.hpp"

OrnsteinUhlenbeckProcesses::OrnsteinUhlenbeckProcesses() {  // Default (empty) constructor
}

void OrnsteinUhlenbeckProcesses::InitProcesses(std::string folder, int seed, int nSeries,
                                                 std::vector<std::vector<std::string>> modeNames) {
  this->names = modeNames;
  this->nSeries = nSeries;
  this->ouValues = IdefixArray2D<Kokkos::complex<real>>("ouValues", nSeries, COMPONENTS);
  this->normalValuesReal = IdefixArray2D<real>("normalValuesReal", nSeries, COMPONENTS);
  this->normalValuesImag = IdefixArray2D<real>("normalValuesImag", nSeries, COMPONENTS);
  this->means = IdefixArray2D<real>("ouMeans", nSeries, COMPONENTS);
  this->tcorrs = IdefixArray2D<real>("ouTcorrs", nSeries, COMPONENTS);
  this->epsilons = IdefixArray2D<real>("ouEpsilons", nSeries, COMPONENTS);

  // The OU amplitudes start at the (zero) mean. On a fresh start,
  // SetProcesses() resets them to the actual mean; on a restart, the dump
  // overwrites them with the exact values at the time of the dump.
  Kokkos::deep_copy(this->ouValues, 0);

  // Deterministic initialisation of one xorshift64 state per (mode, component)
  // from the configured seed (see agent/adr/0002 and idfx::RandomXorShift64).
  // The generator stores a flat 1D state; the (l, dir) pair is mapped to the
  // effective index l*COMPONENTS+dir (the LayoutRight offset of (l, dir)).
  this->rng = idfx::RandomXorShift64(nSeries*COMPONENTS);
  this->rng.InitFromSeed(seed);

  this->ouFilename = folder + "/ou_prank" + std::to_string(idfx::prank)
                     + "_seed" + std::to_string(seed) + ".dat";
  this->normalFilename = folder + "/normal_prank" + std::to_string(idfx::prank)
                         + "_seed" + std::to_string(seed) + ".dat";
  this->precision = 10;
  this->ouValuesHost =
      IdefixHostArray2D<Kokkos::complex<real>>("ouValuesHost", nSeries, COMPONENTS);
  this->normalValuesRealHost =
      IdefixHostArray2D<real>("normalValuesRealHost", nSeries, COMPONENTS);
  this->normalValuesImagHost =
      IdefixHostArray2D<real>("normalValuesImagHost", nSeries, COMPONENTS);
}

void OrnsteinUhlenbeckProcesses::SetProcesses(IdefixArray2D<real> mean,
                                                IdefixArray2D<real> tcorr,
                                                IdefixArray2D<real> epsilon) {
  IdefixArray2D<real> means = this->means;
  IdefixArray2D<real> tcorrs = this->tcorrs;
  IdefixArray2D<real> epsilons = this->epsilons;
  IdefixArray2D<Kokkos::complex<real>> ouValues = this->ouValues;
  int nSeries = this->nSeries;
  idefix_for("SetProcesses", 0, nSeries, 0, COMPONENTS,
              KOKKOS_LAMBDA (int l, int dir) {
        means(l, dir) = mean(l, dir);
        tcorrs(l, dir) = tcorr(l, dir);
        // so that the total amplitude is independent of the number of modes used to force
        epsilons(l, dir) = epsilon(l, dir)/nSeries;
        ouValues(l, dir) = mean(l, dir);
  });
}

void OrnsteinUhlenbeckProcesses::UpdateProcessesValues(real dt) {
  IdefixArray2D<real> means = this->means;
  IdefixArray2D<real> tcorrs = this->tcorrs;
  IdefixArray2D<real> epsilons = this->epsilons;
  IdefixArray2D<Kokkos::complex<real>> ouValues = this->ouValues;
  idfx::RandomXorShift64 rng = this->rng;
  IdefixArray2D<real> normalValuesReal = this->normalValuesReal;
  IdefixArray2D<real> normalValuesImag = this->normalValuesImag;
  idefix_for("UpdateProcesses", 0, nSeries, 0, COMPONENTS,
              KOKKOS_LAMBDA (int l, int dir) {
      // Each (mode, component) draws from its own xorshift64 state: the
      // sequence is deterministic on every backend and identical on all
      // ranks, and it is part of the restart dump (see agent/adr/0002).
      // The by-value copy of rng shares the same underlying state array.
      // Effective 1D index of the (l, dir) stream in the flat state array.
      const int idx = l*COMPONENTS + dir;
      real normalReal = static_cast<real>(rng.Normal(idx));
      real normalImag = static_cast<real>(rng.Normal(idx));
      normalValuesReal(l, dir) = normalReal;
      normalValuesImag(l, dir) = normalImag;
      real expTerm = std::exp(-dt/tcorrs(l, dir));
      real dou = std::sqrt(epsilons(l, dir)/sqrt(2.)/tcorrs(l, dir)*(1. - expTerm*expTerm));
      real realPart = ouValues(l, dir).real();
      real normalisedMean = means(l, dir)/sqrt(2.);
      real newRealPart = normalisedMean + (realPart-normalisedMean)*expTerm + dou*normalReal;

      real imagPart = ouValues(l, dir).imag();
      real newImagPart = normalisedMean + (imagPart-normalisedMean)*expTerm + dou*normalImag;
      Kokkos::complex<real> newCompValue(newRealPart, newImagPart);
      ouValues(l, dir) = newCompValue;
  });
}

void OrnsteinUhlenbeckProcesses::ResetProcessesValues() {
  if (idfx::prank==0) {
    file.open(ouFilename, std::ios::trunc);
    int col_width = 3*precision + 10;
    file << "t";
    for (int l=0; l<nSeries; l++) {
      for (int dir=0; dir<COMPONENTS; dir++) {
        std::string current_name = names[l][dir];
        file << std::setw(col_width) << current_name;
      }
    }
    file << std::endl;
    file.close();
  }
}

void OrnsteinUhlenbeckProcesses::WriteProcessesValues(real time) {
  if (idfx::prank==0) {
    int col_width = 3*precision + 10;
    file.open(ouFilename, std::ios::app);
    file.precision(precision);
    this->file << std::scientific << time;
    Kokkos::deep_copy(ouValuesHost, ouValues);
    for (int l=0; l<nSeries; l++) {
      for (int dir=0; dir<COMPONENTS; dir++) {
        file << std::scientific << std::setw(col_width) << ouValuesHost(l,dir).real()
             << '+' << ouValuesHost(l,dir).imag() << 'j';
      }
    }
    file << std::endl;
    file.close();
  }
}

void OrnsteinUhlenbeckProcesses::ResetNormalValues() {
  if (idfx::prank==0) {
    file.open(normalFilename, std::ios::trunc);
    int col_width = 3*precision + 10;
    file << std::setw(col_width) << "t";
    for (int l=0; l<nSeries; l++) {
      for (int dir=0; dir<COMPONENTS; dir++) {
        std::string current_name = names[l][dir];
        file << std::setw(col_width) << current_name;
      }
    }
    file << std::endl;
    file.close();
  }
}

void OrnsteinUhlenbeckProcesses::WriteNormalValues(real time) {
  if (idfx::prank==0) {
    int col_width = precision + 10;
    file.open(normalFilename, std::ios::app);
    file.precision(precision);
    this->file << std::scientific << std::setw(col_width) << time;
    Kokkos::deep_copy(normalValuesRealHost, normalValuesReal);
    Kokkos::deep_copy(normalValuesImagHost, normalValuesImag);
    for (int l=0; l<nSeries; l++) {
      for (int dir=0; dir<COMPONENTS; dir++) {
        file << std::scientific << std::setw(col_width) << normalValuesRealHost(l,dir)
             << '+' << normalValuesImagHost(l,dir) << 'j';
      }
    }
    file << std::endl;
    file.close();
  }
}
