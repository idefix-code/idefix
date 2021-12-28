// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2022 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef MPI_HPP_
#define MPI_HPP_

#include <vector>
#include "idefix.hpp"

#ifdef WITH_MPI
#include "mpi.hpp"
#endif

class DataBlock;

class Mpi {
 public:
  // MPI Exchange functions
  void ExchangeAll();   ///< Exchange boundary elements in all directions (todo)
  void ExchangeX1();    ///< Exchange boundary elements in the X1 direction
  void ExchangeX2();    ///< Exchange boundary elements in the X2 direction
  void ExchangeX3();    ///< Exchange boundary elements in the X3 direction

  // Init from datablock
  void Init(DataBlock *datain, IdefixArray4D<real> inputVc, std::vector<int> inputMap,
            bool inputHaveVs = false, IdefixArray4D<real> inputVs = IdefixArray4D<real>() );

  // Destructor
  ~Mpi();

 private:
  static int nInstances;     // total number of mpi instances in the code
  int thisInstance;          // unique number of the current instance
  bool isInitialized{false};

  DataBlock *data;          // pointer to datablock object
  IdefixArray4D<real> Vc;   // reference to cell-centered array on which this object works
  IdefixArray4D<real> Vs;   // reference to face-centered array on which this object works

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
  int64_t bytesSentOrReceived{0};
};

#endif // MPI_HPP_
