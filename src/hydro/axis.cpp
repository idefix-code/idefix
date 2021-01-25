// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2021 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include "axis.hpp"
#include "hydro.hpp"


void Init(Grid &grid, Hydro *h) {
  this->hydro = h;
  // Do we have a full circle?
  // Check that we have a fraction of 2PI:
  double should_be_integer = 2.0*M_PI/fabs(grid.xend[KDIR] - grid.xbeg[KDIR]);

  if(fabs(should_be_integer - round(should_be_integer))>1e-10) {
    IDEFIX_ERROR("The grid extent in X3 should be an integer fraction of 2Pi");
  }

  idfx::cout << "Hydro::Axis: Axis regularisation enabled ";

  if(fabs((grid.xend[JDIR] - grid.xbeg[IDIR] -2.0*M_PI)) < 1e-10) {
    this->isTwoPi = true;
    idfx::cout << "with full azimuthal extension"
  } else {
    this->isTwoPi = false;
    idfx::cout << "with partial azimuthal extension"
  }

  // Check where the axis is lying.
  if(hydro->data->lbound[JDIR] == axis) axisLeft = true;
  if(hydro->data->rbound[JDIR] == axis) axisRight = true;

  #if MHD == YES
  this->Ex1Avg = IdefixArray1D<real>("Axis:Ex1Avg",hydro->data->np_tot[IDIR]);
  #endif
  idfx::cout << std::endl;
}

// Average the Emf component along the axis
void SymmetrizeEx1() {
}

// enforce the boundary conditions on the ghost zone accross the axis
void EnforceAxisBoundary(int side) {
}
