// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2022 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2021 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include "OrnsteinUhlenbeckProcess.hpp"

OrnsteinUhlenbeckProcesses::OrnsteinUhlenbeckProcesses()
{ // Default (empty) constructor
}

void OrnsteinUhlenbeckProcesses::InitProcesses(std::string folder, int seed, int nSeries, std::vector<std::vector<std::string>> modeNames) {
  this->names = modeNames;
  this->nSeries = nSeries;
  this->ouValues = IdefixArray2D<Kokkos::complex<real>> ("ouValues", nSeries, COMPONENTS);
  this->normalValuesReal = IdefixArray2D<real> ("normalValuesReal", nSeries, COMPONENTS);
  this->normalValuesImag = IdefixArray2D<real> ("normalValuesImag", nSeries, COMPONENTS);
  this->random_pool = Kokkos::Random_XorShift64_Pool<> (/*seed=*/seed);

  this->ouFilename = folder + "/ou_prank" + std::to_string(idfx::prank) + "_seed" + std::to_string(seed) + ".dat";
  this->normalFilename = folder + "/normal_prank" + std::to_string(idfx::prank) + "_seed" + std::to_string(seed) + ".dat";
  this->timestepFilename = folder + "/timestep.dat";
  this->precision = 10;
  this->ouValuesHost = IdefixHostArray2D<Kokkos::complex<real>> ("ouValuesHost", nSeries, COMPONENTS);
  this->normalValuesRealHost = IdefixHostArray2D<real> ("normalValuesRealHost", nSeries, COMPONENTS);
  this->normalValuesImagHost = IdefixHostArray2D<real> ("normalValuesImagHost", nSeries, COMPONENTS);
}

void OrnsteinUhlenbeckProcesses::SetProcesses(IdefixArray2D<real> mean, IdefixArray2D<real> tcorr, IdefixArray2D<real> epsilon) {
  this->epsilons = IdefixArray2D<real> ("ouEpsilons", nSeries, COMPONENTS);
  this->tcorrs = IdefixArray2D<real> ("ouTcorrs", nSeries, COMPONENTS);
  this->means = IdefixArray2D<real> ("ouMeans", nSeries, COMPONENTS);
  IdefixArray2D<real> means = this->means;
  IdefixArray2D<real> tcorrs = this->tcorrs;
  IdefixArray2D<real> epsilons = this->epsilons;
  IdefixArray2D<Kokkos::complex<real>> ouValues = this->ouValues;
  int nSeries = this->nSeries;
  idefix_for("SetProcesses", 0, nSeries, 0, COMPONENTS,
              KOKKOS_LAMBDA (int l, int dir) {
        means(l, dir) = mean(l, dir);
        tcorrs(l, dir) = tcorr(l, dir);
        epsilons(l, dir) = epsilon(l, dir)/nSeries; //so that the total amplitude is independent of the number of modes used to force
        ouValues(l, dir) = mean(l, dir);
  });
}

//void OrnsteinUhlenbeckProcesses::UpdateProcessesValues(real dt, IdefixArray1D<real> epsilons) {
void OrnsteinUhlenbeckProcesses::UpdateProcessesValues(real dt) {
  IdefixArray2D<real> means = this->means;
  IdefixArray2D<real> tcorrs = this->tcorrs;
  IdefixArray2D<real> epsilons = this->epsilons;
  IdefixArray2D<Kokkos::complex<real>> ouValues = this->ouValues;
  IdefixArray2D<real> normalValuesReal = this->normalValuesReal;
  IdefixArray2D<real> normalValuesImag = this->normalValuesImag;
  Kokkos::Random_XorShift64_Pool<> random_pool = this->random_pool;
  idefix_for("UpdateProcesses", 0, nSeries, 0, COMPONENTS,
              KOKKOS_LAMBDA (int l, int dir) {
      auto generator = random_pool.get_state();
      real normalReal = generator.normal(0., 1.);
      real normalImag = generator.normal(0., 1.);
      random_pool.free_state(generator);
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

//void OrnsteinUhlenbeckProcesses::AdvanceProcessesValues(std::vector<real> tabDt) {
//std::cout << tabDt.empty() << std::endl;
//  while (not tabDt.empty()) {
//    real dt = tabDt[0];
//    tabDt.erase(tabDt.begin());
//    UpdateProcessesValues(dt);
////std::cout << dt << std::endl;
//  }
//}

void OrnsteinUhlenbeckProcesses::AdvanceProcessesValues() {

  std::ifstream file(timestepFilename);
  std::string line;

  if (file.is_open()) {
    while (getline(file, line)) {
      std::string::size_type sz;
      real dt = std::stof(line, &sz);
      real t = std::stof(line.substr(sz));
      UpdateProcessesValues(dt);
    }
    file.close();
  }
  else {
      IDEFIX_WARNING("UNABLE TO READ THE TIMESTEP FILE TO ADVANCE OU PROCESSES");
  }
}

void OrnsteinUhlenbeckProcesses::ResetProcessesValues() {
  if(idfx::prank==0) {
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
  if(idfx::prank==0) {
    int col_width = 3*precision + 10;
    file.open(ouFilename, std::ios::app);
    file.precision(precision);
    this->file << std::scientific << time;
    Kokkos::deep_copy(ouValuesHost, ouValues);
    for (int l=0; l<nSeries; l++) {
      for (int dir=0; dir<COMPONENTS; dir++) {
//        file << std::scientific << std::setw(col_width) << ouValuesHost(l,dir);
        file << std::scientific << std::setw(col_width) << ouValuesHost(l,dir).real() << '+' << ouValuesHost(l,dir).imag() << 'j';
      }
    }
    file << std::endl;
    file.close();
  }
}

void OrnsteinUhlenbeckProcesses::ResetNormalValues() {
  if(idfx::prank==0) {
    file.open(normalFilename, std::ios::trunc);
    int col_width = 3*precision + 10;
    file << std::setw(col_width) << "t";
    for (int l=0; l<nSeries; l++) {
      for (int dir=0; dir<COMPONENTS; dir++) {
//        std::string current_name = std::to_string(l) + std::to_string(dir);
        std::string current_name = names[l][dir];
        file << std::setw(col_width) << current_name;
      }
    }
    file << std::endl;
    file.close();
  }
}

void OrnsteinUhlenbeckProcesses::WriteNormalValues(real time) {
  if(idfx::prank==0) {
    int col_width = precision + 10;
    file.open(normalFilename, std::ios::app);
    file.precision(precision);
    this->file << std::scientific << std::setw(col_width) << time;
    Kokkos::deep_copy(normalValuesRealHost, normalValuesReal);
    Kokkos::deep_copy(normalValuesImagHost, normalValuesImag);
    for (int l=0; l<nSeries; l++) {
      for (int dir=0; dir<COMPONENTS; dir++) {
//        file << std::scientific << std::setw(col_width) << normalValuesHost(l, dir);
        file << std::scientific << std::setw(col_width) << normalValuesRealHost(l,dir) << '+' << normalValuesImagHost(l,dir) << 'j';
      }
    }
    file << std::endl;
    file.close();
  }
}

void OrnsteinUhlenbeckProcesses::ResetTimestep() {
  if(idfx::prank==0) {
    file.open(timestepFilename, std::ios::trunc);
    file.close();
  }
}

void OrnsteinUhlenbeckProcesses::WriteTimestep(real time, real dt) {
  if(idfx::prank==0) {
    int col_width = precision + 10;
    file.open(timestepFilename, std::ios::app);
    file.precision(precision);
    this->file << std::scientific << std::setw(col_width) << dt;
    this->file << std::scientific << std::setw(col_width) << time;
    file << std::endl;
    file.close();
  }
}

