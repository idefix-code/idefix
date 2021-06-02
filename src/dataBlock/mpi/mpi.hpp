// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2021 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef DATABLOCK_MPI_MPI_HPP_
#define DATABLOCK_MPI_MPI_HPP_

#include "idefix.hpp"

class DataBlock;

class Mpi {
 public:
  // MPI Exchange functions
  void ExchangeAll();
  void ExchangeX1();
  void ExchangeX2();
  void ExchangeX3();

  // constructor from datablock
  Mpi(DataBlock *, IdefixArray1D<int>&, int, bool);

  // Destructor
  ~Mpi();

 private:
  static int nInstances;     // total number of mpi instances in the code
  int thisInstance;          // unique number of the current instance
  
  DataBlock *data;

  enum {faceRight, faceLeft};

  // Buffers for MPI calls
  IdefixArray1D<real> BufferSendX1[2];
  IdefixArray1D<real> BufferSendX2[2];
  IdefixArray1D<real> BufferSendX3[2];
  IdefixArray1D<real> BufferRecvX1[2];
  IdefixArray1D<real> BufferRecvX2[2];
  IdefixArray1D<real> BufferRecvX3[2];

  IdefixArray1D<int>  mapVars;
  int mapNVars{0};

  int bufferSizeX1;
  int bufferSizeX2;
  int bufferSizeX3;

  bool haveVs{false};

  // Requests for MPI persistent communications
  MPI_Request sendRequestX1[2];
  MPI_Request sendRequestX2[2];
  MPI_Request sendRequestX3[2];
  MPI_Request recvRequestX1[2];
  MPI_Request recvRequestX2[2];
  MPI_Request recvRequestX3[2];

  Grid *mygrid;

  // MPI throughput timer specific to this object
  double myTimer{0};
  long int bytesSentOrReceived{0};
};

#endif // DATABLOCK_MPI_MPI_HPP_
