// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef FORCING_FORCING_HPP_
#define FORCING_FORCING_HPP_

#include "idefix.hpp"
#include "input.hpp"
#include "OrnsteinUhlenbeckProcess.hpp"

enum ForcingType {iso3D, iso2D};

class DataBlock;

class Forcing {
 public:
  Forcing(Input&, DataBlock*);

  bool stillHaveForcing;
  real stopTime;

  void InitForcingParameters();          ///< init forcing parameters

  void InitForcingModes();          ///< init forcing modes given its type

  void ComputeForcing(real);        ///< compute the required forcing field at current time t

  void ComputePristineForcing(real);///< compute the pristine forcing field at current time t

  void ComputeSolenoidalForcing(real); ///< compute the solenoidal part of the forcing field at current time t

//  void ComputeCompressiveForcing(real); ///< compute the compressive part of the forcing field at current time t

  void ResetForcingTerms();          ///< fill the forcing field with zeros.

  void ShowConfig();                ///< Show the forcing configuration

  // Forcing modes
//  IdefixArray4D<Kokkos::complex<real>> forcingModes;
  EXPAND(
  IdefixArray4D<Kokkos::complex<real>> forcingModesIdir; ,
  IdefixArray4D<Kokkos::complex<real>> forcingModesJdir; ,
  IdefixArray4D<Kokkos::complex<real>> forcingModesKdir; )
//  IdefixHostArray2D<std::string> modeNamesHost;
  std::vector<std::vector<std::string>> modeNames;
  // Forcing terms
  IdefixArray4D<real> forcingTerm;
  IdefixArray4D<real> pristineForcingTerm;
  IdefixArray4D<real> solenoidalForcingTerm;
//  IdefixArray4D<real> compressiveForcingTerm;
  OrnsteinUhlenbeckProcesses oUprocesses;

  void WriteNormalBasis(std::string);                ///< Write the functions of the normal basis in txt files

  int write;
 private:

  DataBlock *data;
  int seed;
  std::string folder;

  int nForcingModes;
  real machNumber;
  real cs;
  real tcorr;
  real epsilon;
  IdefixArray2D<real> tcorrs;
  IdefixArray2D<real> means;
  IdefixArray2D<real> epsilons;
  IdefixHostArray2D<real> hostTcorrs;
  IdefixHostArray2D<real> hostMeans;
  IdefixHostArray2D<real> hostEpsilons;

  ForcingType forcingType;
  int normal2Diso;
  std::string normal2DisoStr;
  int haveSolenoidalForcing;

  IdefixArray2D<real> k3Diso;
  IdefixHostArray2D<real> k3DisoHost;
  IdefixArray2D<real> k2Diso;
  IdefixHostArray2D<real> k2DisoHost;
  IdefixArray2D<real> k3Dani;
  IdefixHostArray2D<real> k3DaniHost;
  IdefixArray2D<int> ellmVsh;
  IdefixHostArray2D<int> ellmVshHost;
  real kmin;
  real kmax;
  real kx0;
  real ky0;
  real kz0;
  real xbeg;
  real ybeg;
  real zbeg;
  real xend;
  real yend;
  real zend;
  int ellmin;
  int ellmax;
  int mmin;
  int mmax;
//  #endif //GEOMETRY == SPHERICAL
};

#endif // FORCING_FORCING_HPP_
