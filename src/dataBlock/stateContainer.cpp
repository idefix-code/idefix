// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2022 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include "stateContainer.hpp"
#include "idefix.hpp"
#include <vector>
#include <string>

StateContainer::StateContainer() {
  // do nothing
}

StateContainer StateContainer::Copy() {
  // Return a deepcopy of the current state container
  StateContainer sc;
  for(State stateIn : this->stateVector) {
    // We first do a shallow copy
    State stateOut = stateIn;
    // But then reinit the array
    stateOut.array = IdefixArray4D<real>(stateIn.stateName, stateIn.array.extent(0),
                                                            stateIn.array.extent(1),
                                                            stateIn.array.extent(2),
                                                            stateIn.array.extent(3));
    // A do a deep copy
    Kokkos::deep_copy(stateOut.array, stateIn.array);
    sc.stateVector.push_back(stateOut);
  }
}    

void StateContainer::PushArray(&IdefixArray4D<real> in, State::TypeLocalisation loc, std::string name) {
  idfx::pushRegion("StateContainer::PushArray");
  State state;
  state.array = in;
  state.localisation = loc;
  state.name = name;
  this->stateVector.push_back(state);
  idfx::popRegion();
}

void StateContainer::AddAndStore(&StateContainer, real, real) {

}