// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef HYDRO_ELECTROMOTIVEFORCE_HPP_
#define  HYDRO_ELECTROMOTIVEFORCE_HPP_

#include "idefix.hpp"

// Forward declarations
class DataBlock;

class ElectroMotiveForce {
 public:
  // Face centered emf components
  IdefixArray3D<real>     exj;
  IdefixArray3D<real>     exk;
  IdefixArray3D<real>     eyi;
  IdefixArray3D<real>     eyk;
  IdefixArray3D<real>     ezi;
  IdefixArray3D<real>     ezj;

  // Edge centered emf components
  IdefixArray3D<real>     ex;
  IdefixArray3D<real>     ey;
  IdefixArray3D<real>     ez;

  IdefixArray3D<int>      svx;
  IdefixArray3D<int>      svy;
  IdefixArray3D<int>      svz;

  IdefixArray3D<real>     Ex1;
  IdefixArray3D<real>     Ex2;
  IdefixArray3D<real>     Ex3;

  // Range of existence
  int  ibeg, jbeg, kbeg;
  int  iend, jend, kend;

  // Constructor from datablock structure
  explicit ElectroMotiveForce(DataBlock *);

  // Default constructor
  ElectroMotiveForce();
};

#endif // HYDRO_ELECTROMOTIVEFORCE_HPP_
