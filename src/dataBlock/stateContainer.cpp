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
    stateOut.array = IdefixArray4D<real>(stateIn.name, stateIn.array.extent(0),
                                                            stateIn.array.extent(1),
                                                            stateIn.array.extent(2),
                                                            stateIn.array.extent(3));
    // A do a deep copy
    Kokkos::deep_copy(stateOut.array, stateIn.array);
    sc.stateVector.push_back(stateOut);
  }
  return(sc);
}    

void StateContainer::PushArray(IdefixArray4D<real>& in, State::TypeLocalisation loc, std::string name) {
  idfx::pushRegion("StateContainer::PushArray");
  State state;
  state.array = in;
  state.localisation = loc;
  state.name = name;
  this->stateVector.push_back(state);
  idfx::popRegion();
}


void StateContainer::AddAndStore(real wl, real wr, StateContainer & in) {
  idfx::pushRegion("StateContainer::AddAndStore");

  if(in.stateVector.size() != this->stateVector.size()) {
    IDEFIX_ERROR("You cannot add two state containers of different sizes.");
  }
  for(int s = 0 ; s < this->stateVector.size() ; s++) {
    State stateIn = in.stateVector[s];
    State stateOut = this->stateVector[s];
    auto Vin = stateIn.array;
    auto Vout = stateOut.array;
    idefix_for("StateContainer::AddAndStore", 
                0, Vin.extent(0), 
                0, Vin.extent(1),
                0, Vin.extent(2),
                0, Vin.extent(3),
                KOKKOS_LAMBDA(int n, int k, int j, int i) {
                   Vout(n,k,j,i) = wl * Vout(n,k,j,i) + wr * Vin(n,k,j,i);
                } ); 
  }
  idfx::popRegion();

}