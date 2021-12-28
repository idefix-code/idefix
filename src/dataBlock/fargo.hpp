// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2022 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef DATABLOCK_FARGO_HPP_
#define DATABLOCK_FARGO_HPP_

#include "idefix.hpp"

// Forward class hydro declaration
class Hydro;
class DataBlock;

using FargoVelocityFunc = void (*) (DataBlock &, IdefixArray2D<real> &);

class Fargo {
 public:
  enum FargoType {none, userdef, shearingbox};
  void Init(Input &, DataBlock*);  // Initialisation
  void ShiftSolution(const real t, const real dt);  // Effectively shift the solution
  void SubstractVelocity(const real);
  void AddVelocity(const real);
  void EnrollVelocity(FargoVelocityFunc);
 private:
  friend class Hydro;
  DataBlock *data;
  Hydro *hydro;
  FargoType type{none};                 // By default, Fargo is disabled

  IdefixArray2D<real> meanVelocity;
  IdefixArray4D<real> scratch;

  bool velocityHasBeenComputed{false};
  void GetFargoVelocity(real);
  FargoVelocityFunc fargoVelocityFunc{NULL};  // The user-defined fargo velocity function
};

#endif // DATABLOCK_FARGO_HPP_
