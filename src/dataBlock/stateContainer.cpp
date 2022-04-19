// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2022 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include "stateContainer.hpp"
#include <vector>
#include <string>
#include "idefix.hpp"


StateContainer::StateContainer() {
  // do nothing
}

void StateContainer::CopyFrom(StateContainer &in) {
  idfx::pushRegion("StateContainer::CopyFrom");
  // Make a deep copy from in
  if(in.stateVector.size() != this->stateVector.size()) {
    IDEFIX_ERROR("You cannot copy two state containers of different sizes.");
  }

  for(int s = 0 ; s < this->stateVector.size() ; s++) {
    Kokkos::deep_copy(this->stateVector[s].array, in.stateVector[s].array);
  }
  idfx::popRegion();
}

void StateContainer::AllocateAs(StateContainer &in) {
  idfx::pushRegion("StateContainer::AllocateAs");
  // Allocate arrays of current StateContainer with the same shape as in
  for(State stateIn : in.stateVector) {
    // We first do a shallow copy
    State stateOut = stateIn;
    // But then reinit the array
    stateOut.array = IdefixArray4D<real>(stateIn.name, stateIn.array.extent(0),
                                                            stateIn.array.extent(1),
                                                            stateIn.array.extent(2),
                                                            stateIn.array.extent(3));
    // And add the new state to our stateVector
    this->stateVector.push_back(stateOut);
  }
  idfx::popRegion();
}

void StateContainer::PushArray(IdefixArray4D<real>& in,
                               State::TypeLocation loc,
                               std::string name) {
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
