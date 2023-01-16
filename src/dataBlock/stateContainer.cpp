// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include "stateContainer.hpp"
#include <vector>
#include <string>
#include <algorithm>
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
    if(this->stateVector[s].type != in.stateVector[s].type) {
      IDEFIX_ERROR("Cannot copy from a stateContainer with different types");
    }
    if(this->stateVector[s].type == State::idefixArray4D) {
      Kokkos::deep_copy(this->stateVector[s].array, in.stateVector[s].array);
    } else {
      IDEFIX_ERROR("Cannot copy states which are undefined");
    }
  }
  idfx::popRegion();
}

void StateContainer::AllocateAs(StateContainer &in) {
  idfx::pushRegion("StateContainer::AllocateAs");
  // Allocate arrays of current StateContainer with the same shape as in
  for(State stateIn : in.stateVector) {
    // We first do a shallow copy
    State stateOut = stateIn;

    if(stateIn.type == State::idefixArray4D) {
      // But then reinit the array
      stateOut.array = IdefixArray4D<real>(stateIn.name, stateIn.array.extent(0),
                                                              stateIn.array.extent(1),
                                                              stateIn.array.extent(2),
                                                              stateIn.array.extent(3));
    } else {
      IDEFIX_ERROR("Cannot allocate a state with type none");
    }
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
  state.type = State::idefixArray4D;
  state.name = name;
  state.location = loc;
  this->stateVector.push_back(state);
  idfx::popRegion();
}


void StateContainer::AddAndStore(const real wl, const real wr, StateContainer & in) {
  idfx::pushRegion("StateContainer::AddAndStore");

  if(in.stateVector.size() != this->stateVector.size()) {
    IDEFIX_ERROR("You cannot add two state containers of different sizes.");
  }
  for(int s = 0 ; s < this->stateVector.size() ; s++) {
    State stateIn = in.stateVector[s];
    State stateOut = this->stateVector[s];

    if(stateIn.type != stateOut.type) {
      IDEFIX_ERROR("Cannot add and store states of different type");
    }
    if(stateIn.type == State::idefixArray4D) {
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
    } else {
      IDEFIX_ERROR("Cannot Add and store from state of unknown type");
    }
  }
  idfx::popRegion();
}
