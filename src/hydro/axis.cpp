// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2021 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include "axis.hpp"
#include "hydro.hpp"
#include "dataBlock.hpp"


void Axis::Init(Grid &grid, Hydro *h) {
  this->hydro = h;
  this->data = this->hydro->data;
  this->emf = & this->hydro->emf;

  // Do we have a full circle?
  // Check that we have a fraction of 2PI:
  double should_be_integer = 2.0*M_PI/fabs(grid.xend[KDIR] - grid.xbeg[KDIR]);

  if(fabs(should_be_integer - round(should_be_integer))>1e-10) {
    IDEFIX_ERROR("The grid extent in X3 should be an integer fraction of 2Pi");
  }

  idfx::cout << "Axis:: Axis regularisation enabled ";

  if(fabs((grid.xend[KDIR] - grid.xbeg[KDIR] -2.0*M_PI)) < 1e-10) {
    this->isTwoPi = true;
    idfx::cout << "with full (2pi) azimuthal extension";
  } else {
    this->isTwoPi = false;
    idfx::cout << "with partial (<2pi) azimuthal extension";
  }

  // Check where the axis is lying.
  if(hydro->data->lbound[JDIR] == axis) axisLeft = true;
  if(hydro->data->rbound[JDIR] == axis) axisRight = true;

  #if MHD == YES
  this->Ex1Avg = IdefixArray1D<real>("Axis:Ex1Avg",hydro->data->np_tot[IDIR]);
  #endif
  idfx::cout << std::endl;
}

void Axis::SymmetrizeEx1Side(int jref) {
  IdefixArray3D<real> Ex1 = emf->ex;
  IdefixAtomicArray1D<real> Ex1Avg = this->Ex1Avg;

  idefix_for("Ex1_ini",0,data->np_tot[IDIR],
      KOKKOS_LAMBDA(int i) {
        Ex1Avg(i) = ZERO_F;
      });

  idefix_for("Ex1_Symmetrize",data->beg[KDIR],data->end[KDIR],0,data->np_tot[IDIR],
    KOKKOS_LAMBDA(int k,int i) {
      Ex1Avg(i) += Ex1(k,jref,i);
    });

  int ncells=data->mygrid->np_int[KDIR];

  idefix_for("Ex1_Store",data->beg[KDIR],data->end[KDIR],0,data->np_tot[IDIR],
  KOKKOS_LAMBDA(int k,int i) {
    Ex1(k,jref,i) = Ex1Avg(i)/((real) ncells);
  });
}
// Average the Emf component along the axis
void Axis::SymmetrizeEx1() {
  idfx::pushRegion("Axis::SymmetrizeEx1");

  if(this->axisLeft) {
    int jref = hydro->data->beg[JDIR];
    this->SymmetrizeEx1Side(jref);
  }
  if(this->axisRight) {
    int jref = hydro->data->end[JDIR];
    this->SymmetrizeEx1Side(jref);
  }

  idfx::popRegion();
}

// enforce the boundary conditions on the ghost zone accross the axis
void Axis::EnforceAxisBoundary(int side) {
  idfx::pushRegion("Axis::EnforceAxisBoundary");

  idfx::popRegion();
}
