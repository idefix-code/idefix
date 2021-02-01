// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2021 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef HYDRO_ELECTROMOTIVEFORCE_ELECTROMOTIVEFORCE_HPP_
#define  HYDRO_ELECTROMOTIVEFORCE_ELECTROMOTIVEFORCE_HPP_

#include "idefix.hpp"

// Forward declarations
class Hydro;
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

  // Init from Hydro class
  void Init(Hydro *);

#ifdef WITH_MPI
  // Exchange surface EMFs to remove interprocess round off errors
  void ExchangeAll();
  void ExchangeX1();
  void ExchangeX2();
  void ExchangeX3();
#endif

  // Default constructor
  ElectroMotiveForce();

  // Destructor
  ~ElectroMotiveForce();

 private:
  DataBlock *data;

#ifdef WITH_MPI
  enum {faceRight, faceLeft};

  // Buffers for MPI calls
  IdefixArray1D<real> BufferSendX1[2];
  IdefixArray1D<real> BufferSendX2[2];
  IdefixArray1D<real> BufferSendX3[2];
  IdefixArray1D<real> BufferRecvX1[2];
  IdefixArray1D<real> BufferRecvX2[2];
  IdefixArray1D<real> BufferRecvX3[2];

  IdefixArray1D<int>  mapVars;

  int bufferSizeX1;
  int bufferSizeX2;
  int bufferSizeX3;

  // Requests for MPI persistent communications
  MPI_Request sendRequestX1[2];
  MPI_Request sendRequestX2[2];
  MPI_Request sendRequestX3[2];
  MPI_Request recvRequestX1[2];
  MPI_Request recvRequestX2[2];
  MPI_Request recvRequestX3[2];
#endif
};

#endif // HYDRO_ELECTROMOTIVEFORCE_ELECTROMOTIVEFORCE_HPP_
