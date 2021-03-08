// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2021 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef HYDRO_FARGO_FARGO_HPP_
#define HYDRO_FARGO_FARGO_HPP_

#include "idefix.hpp"

// Forward class hydro declaration
class Hydro;
class DataBlock;

using FargoVelocityFunc = void (*) (DataBlock &, const real t, IdefixArray2D<real> &);
enum FargoType {none, userdef, shearingbox};

class Fargo {
 public:
  void Init(Input &, Grid &, Hydro *);  // Initialisation
  void ShiftSolution(const real t, const real dt);  // Effectively shift the solution
  void EnrollVelocity(FargoVelocityFunc);
 private:
  DataBlock *data;
  Hydro *hydro;
  FargoType type{none};                 // By default, Fargo is disabled

  IdefixArray2D<real> meanVelocity;
  IdefixArray4D<real> scratch;

  FargoVelocityFunc fargoVelocityFunc{NULL};  // The user-defined fargo velocity function
};

#endif // HYDRO_FARGO_FARGO_HPP_
