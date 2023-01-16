// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
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
  enum TypeLocation{undefined, center, face, edge};
  enum TypeState{none, idefixArray4D};

  TypeState type{none};           ///< type of data contained by this state
  IdefixArray4D<real> array;      ///< only defined if type==IdefixArray4D
  TypeLocation location{undefined};    ///< location of array when type==IdefixArray4D
                                       ///< (otherwise undefined)
  std::string name;               ///< Name of the full state (always applicable)
};

class StateContainer {
 public:
  StateContainer();
  void CopyFrom(StateContainer &);    // Return a deepcopy of the current state container
  void AllocateAs(StateContainer &);    // Return a deepcopy of the current state container
  void PushArray(IdefixArray4D<real> &, State::TypeLocation, std::string);
  void AddAndStore(const real, const real, StateContainer&);


 private:
  std::vector<State> stateVector;
};

#endif // DATABLOCK_STATECONTAINER_HPP_
