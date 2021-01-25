// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2021 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef HYDRO_AXIS_HPP_
#define HYDRO_AXIS_HPP_

#include "idefix.hpp"
#include "grid.hpp"
#include "electroMotiveForce.hpp"

// Forward class hydro declaration
class Hydro;
class DataBlock;

class Axis {
 public:
  void Init(Grid &, Hydro *);  // Initialisation
  void SymmetrizeEx1();                 // Symmetrize Emf component Ex1
  void EnforceAxisBoundary(int side);   // Enforce the boundary conditions (along X2)
 private:
  bool isTwoPi = false;
  bool axisRight = false;
  bool axisLeft = false;

  IdefixArray1D<real> Ex1Avg;
  Hydro *hydro;
  DataBlock *data;
  ElectroMotiveForce *emf;

  void SymmetrizeEx1Side(int);
};

#endif // HYDRO_AXIS_HPP_
