// ***********************************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************************

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

  // Init from dataBlock
  void InitFromDataBlock(DataBlock *);

  // Destructor
  ~Mpi();

 private:
  DataBlock *data;

  enum {faceRight, faceLeft};

  // Buffers for MPI calls
  IdefixArray1D<real> BufferSendX1[2];
  IdefixArray1D<real> BufferSendX2[2];
  IdefixArray1D<real> BufferSendX3[2];
  IdefixArray1D<real> BufferRecvX1[2];
  IdefixArray1D<real> BufferRecvX2[2];
  IdefixArray1D<real> BufferRecvX3[2];

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

  Grid *mygrid;

  Kokkos::Timer timer;    // Internal MPI timer
};

#endif // DATABLOCK_MPI_MPI_HPP_
