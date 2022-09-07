// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2022 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef HYDRO_AXIS_HPP_
#define HYDRO_AXIS_HPP_

#include <vector>
#include "idefix.hpp"
#include "grid.hpp"
#include "electroMotiveForce.hpp"

// Forward class hydro declaration
class Hydro;
class DataBlock;

class Axis {
 public:
  void Init(Grid &, Hydro *);  // Initialisation
  void RegularizeEMFs();                 // Regularize the EMF sitting on the axis
  void RegularizeCurrent();             // Regularize the currents along the axis
  void EnforceAxisBoundary(int side);   // Enforce the boundary conditions (along X2)
  void ReconstructBx2s();               // Reconstruct BX2s in the ghost zone using divB=0
  void ShowConfig();


  void SymmetrizeEx1Side(int);         // Symmetrize on a specific side (internal method)
  void RegularizeEx3side(int);         // Regularize Ex3 along the axis (internal method)
  void RegularizeCurrentSide(int);      // Regularize J along the axis (internal method)
  void FixBx2sAxis(int side);           // Fix BX2s on the axis using the field around it (internal)
  void ExchangeMPI(int side);           // Function has to be public for GPU, but its technically
                                        // a private function


 private:
  bool isTwoPi = false;
  bool axisRight = false;
  bool axisLeft = false;
  bool needMPIExchange = false;

  enum {faceTop, faceBot};
#ifdef WITH_MPI
  MPI_Request sendRequest;
  MPI_Request recvRequest;

  IdefixArray1D<real> bufferSend;
  IdefixArray1D<real> bufferRecv;

  int bufferSize;

  IdefixArray1D<int>  mapVars;
  int mapNVars{0};

#endif
  void InitMPI();

  IdefixArray1D<real> Ex1Avg;
  IdefixArray2D<real> BAvg;
  IdefixArray2D<real> JAvg;
  IdefixArray1D<int> symmetryVc;
  IdefixArray1D<int> symmetryVs;


  Hydro *hydro;
  DataBlock *data;
  ElectroMotiveForce *emf;
};

#endif // HYDRO_AXIS_HPP_
