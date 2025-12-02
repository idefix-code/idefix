// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef MPI_EXCHANGER_HPP_
#define MPI_EXCHANGER_HPP_

#include <mpi.h>

#include <signal.h>
#include <vector>
#include <utility>
#include "idefix.hpp"
#include "grid.hpp"
#include "buffer.hpp"
#include "arrays.hpp"

class Grid;

class Exchanger {
 public:
  Exchanger() = default;
  void Init(  Grid* grid,
              int direction,
              std::vector<int> inputMap,
              std::array<int, 3> nghost,
              std::array<int, 3> nint,
              bool inputHaveVs = false,
              bool overwriteBXn = true);

  void Exchange(IdefixArray4D<real> Vc, IdefixArray4D<real> Vs);
  ~Exchanger();

  bool isInitialized{false};

  // MPI throughput timer specific to this object
  double myTimer{0};
  int64_t bytesSentOrReceived{0};

  // Buffer sizes for throughput calculations
  int bufferSize[2];

 private:
  enum {faceRight, faceLeft};

  std::array<BoundingBox,2> boxSend, boxRecv; // bounding boxes for each face
  std::array<std::array<BoundingBox,2>,3> boxSendVs, boxRecvVs; // 3= 3 field components

  // Buffers for MPI calls
  Buffer BufferSend[2];
  Buffer BufferRecv[2];

  int procSend[2];  // MPI process to send to in X1 direction
  int procRecv[2];  // MPI process to receive from in X1 direction

  int direction;
  IdefixArray1D<int>  mapVars;
  int mapNVars{0};

  int nint[3];            //< number of internal elements of the arrays we treat
  int nghost[3];          //< number of ghost zone of the arrays we treat
  int ntot[3];            //< total number of cells of the arrays we treat
  int beg[3];             //< begining index of the active zone
  int end[3];             //< end index of the active zone

  bool haveVs{false};
  bool overwriteBXn{true};

  // Requests for MPI persistent communications
  MPI_Request sendRequest[2];
  MPI_Request recvRequest[2];

  Grid *grid;
};


#endif // MPI_EXCHANGER_HPP_
