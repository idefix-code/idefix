/ ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2022 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef DATABLOCK_STATECONTAINER_HPP_
#define DATABLOCK_STATECONTAINER_HPP_

#include <vector>
#include <string>
#include "idefix.hpp"

class State{
 public:
  enum TypeLocalisation{none, center, face, edge};

  IdefixArray4D<real> array;
  TypeLocalisation localisation;
  std::vector<std::string> arraysName; // Name of each 3D subarrays (when applicable)
  std::string stateName;               // Name of the fulle state (always applicable)
};

class StateContainer {
 public:
  StateContainer();
  StateContainer Copy();    // Return a deepcopy of the current state container
  void PushStateArray(&IdefixArray4D<real>, State::TypeLocalisation, std::string);
  void AddAndStore(&StateContainer, real, real);


 private:
  std::vector<State> stateVector;
};

#endif // DATABLOCK_STATECONTAINER_HPP_
