#include "idefix.hpp"
#include "globals.hpp"

// Constructor for global variables
Globals::Globals() {
    // Allocate global variables on host and device
    this->inputParam = IdefixArray1D<real>("inputParam",32);
    this->inputParamH = create_mirror(this->inputParam);

    // Initialise inputParma on Host
    this->inputParamH(0) = HALF_F;

    // Copy data structure to device
    Kokkos::deep_copy(this->inputParam,this->inputParamH);

}