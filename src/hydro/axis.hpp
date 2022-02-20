// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2022 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
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
  void SymmetrizeEx1Side(int);         // Symmetrize on a specific side (internal method)
  void EnforceAxisBoundary(int side);   // Enforce the boundary conditions (along X2)
  void ReconstructBx2s();


 private:
  bool isTwoPi = false;
  bool axisRight = false;
  bool axisLeft = false;
  bool needMPIExchange = false;

  enum {faceTop, faceBot};
#ifdef WITH_MPI
  MPI_Comm axisComm;
  std::vector<MPI_Datatype> typeVcSend;
  std::vector<MPI_Datatype> typeVcRecv;
#endif

  IdefixArray1D<real> Ex1Avg;
  IdefixArray1D<int> symmetryVc;
  IdefixArray1D<int> symmetryVs;


  Hydro *hydro;
  DataBlock *data;
  ElectroMotiveForce *emf;
};

#endif // HYDRO_AXIS_HPP_
