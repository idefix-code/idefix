// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2022 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef MPI_HPP_
#define MPI_HPP_

#include <signal.h>
#include <vector>
#include "idefix.hpp"
#include "grid.hpp"


class DataBlock;

class Mpi {
 public:
  // MPI Exchange functions
  void ExchangeAll();   ///< Exchange boundary elements in all directions (todo)
  void ExchangeX1(IdefixArray4D<real> inputVc,
                  IdefixArray4D<real> inputVs = IdefixArray4D<real>());
                                      ///< Exchange boundary elements in the X1 direction
  void ExchangeX2(IdefixArray4D<real> inputVc,
                IdefixArray4D<real> inputVs = IdefixArray4D<real>());
                                    ///< Exchange boundary elements in the X2 direction
  void ExchangeX3(IdefixArray4D<real> inputVc,
                IdefixArray4D<real> inputVs = IdefixArray4D<real>());
                                      ///< Exchange boundary elements in the X3 direction

  // Init from datablock
  void Init(Grid *grid, std::vector<int> inputMap,
            int nghost[3], int nint[3], bool inputHaveVs = false );

  // Check that MPI will work with the designated target (in particular GPU Direct)
  static void CheckConfig();


  // Destructor
  ~Mpi();

 private:
  static int nInstances;     // total number of mpi instances in the code
  int thisInstance;          // unique number of the current instance
  bool isInitialized{false};

  DataBlock *data;          // pointer to datablock object

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

  int nint[3];            //< number of internal elements of the arrays we treat
  int nghost[3];          //< number of ghost zone of the arrays we treat
  int ntot[3];            //< total number of cells of the arrays we treat
  int beg[3];             //< begining index of the active zone
  int end[3];             //< end index of the active zone

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

  // Error handler used by CheckConfig
  static void SigErrorHandler(int, siginfo_t* , void* );
};

#endif // MPI_HPP_
