// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef MPI_HPP_
#define MPI_HPP_

#include <signal.h>
#include <vector>
#include <utility>
#include "idefix.hpp"
#include "grid.hpp"


class DataBlock;
class Buffer {
 public:
  Buffer() = default;
  explicit Buffer(size_t size): pointer{0}, array{IdefixArray1D<real>("BufferArray",size)} { };

  void* data() {
    return(array.data());
  }

  int Size() {
    return(array.size());
  }

  void ResetPointer() {
    this->pointer = 0;
  }

  void Pack(IdefixArray3D<real>& in,
       std::pair<int,int> ib,
       std::pair<int,int> jb,
       std::pair<int,int> kb) {
    const int ni = ib.second-ib.first;
    const int ninj = (jb.second-jb.first)*ni;
    const int ninjnk = (kb.second-kb.first)*ninj;
    const int ibeg = ib.first;
    const int jbeg = jb.first;
    const int kbeg = kb.first;
    const int offset = this->pointer;

    auto arr = this->array;
    idefix_for("LoadBuffer3D",kb.first,kb.second,jb.first,jb.second,ib.first,ib.second,
      KOKKOS_LAMBDA (int k, int j, int i) {
      arr(i-ibeg + (j-jbeg)*ni + (k-kbeg)*ninj + offset ) = in(k,j,i);
    });

    // Update pointer
    this->pointer += ninjnk;
  }

  void Pack(IdefixArray4D<real>& in,
       const int var,
       std::pair<int,int> ib,
       std::pair<int,int> jb,
       std::pair<int,int> kb) {
    const int ni = ib.second-ib.first;
    const int ninj = (jb.second-jb.first)*ni;
    const int ninjnk = (kb.second-kb.first)*ninj;
    const int ibeg = ib.first;
    const int jbeg = jb.first;
    const int kbeg = kb.first;
    const int offset = this->pointer;

    auto arr = this->array;
    idefix_for("LoadBuffer4D",kb.first,kb.second,jb.first,jb.second,ib.first,ib.second,
      KOKKOS_LAMBDA (int k, int j, int i) {
      arr(i-ibeg + (j-jbeg)*ni + (k-kbeg)*ninj + offset ) = in(var, k,j,i);
    });

    // Update pointer
    this->pointer += ninjnk;
  }

  void Pack(IdefixArray4D<real>& in,
       IdefixArray1D<int>& map,
       std::pair<int,int> ib,
       std::pair<int,int> jb,
       std::pair<int,int> kb) {
    const int ni = ib.second-ib.first;
    const int ninj = (jb.second-jb.first)*ni;
    const int ninjnk = (kb.second-kb.first)*ninj;
    const int ibeg = ib.first;
    const int jbeg = jb.first;
    const int kbeg = kb.first;
    const int offset = this->pointer;
    auto arr = this->array;

    idefix_for("LoadBuffer4D",0,map.size(),
                             kb.first,kb.second,
                             jb.first,jb.second,
                             ib.first,ib.second,
      KOKKOS_LAMBDA (int n, int k, int j, int i) {
      arr(i-ibeg + (j-jbeg)*ni + (k-kbeg)*ninj + n*ninjnk + offset ) = in(map(n), k,j,i);
    });

    // Update pointer
    this->pointer += ninjnk*map.size();
  }

  void Unpack(IdefixArray3D<real>& out,
       std::pair<int,int> ib,
       std::pair<int,int> jb,
       std::pair<int,int> kb) {
    const int ni = ib.second-ib.first;
    const int ninj = (jb.second-jb.first)*ni;
    const int ninjnk = (kb.second-kb.first)*ninj;
    const int ibeg = ib.first;
    const int jbeg = jb.first;
    const int kbeg = kb.first;
    const int offset = this->pointer;
    auto arr = this->array;

    idefix_for("LoadBuffer3D",kb.first,kb.second,jb.first,jb.second,ib.first,ib.second,
      KOKKOS_LAMBDA (int k, int j, int i) {
        out(k,j,i) = arr(i-ibeg + (j-jbeg)*ni + (k-kbeg)*ninj + offset );
    });

    // Update pointer
    this->pointer += ninjnk;
  }

  void Unpack(IdefixArray4D<real>& out,
       const int var,
       std::pair<int,int> ib,
       std::pair<int,int> jb,
       std::pair<int,int> kb) {
    const int ni = ib.second-ib.first;
    const int ninj = (jb.second-jb.first)*ni;
    const int ninjnk = (kb.second-kb.first)*ninj;
    const int ibeg = ib.first;
    const int jbeg = jb.first;
    const int kbeg = kb.first;
    const int offset = this->pointer;

    auto arr = this->array;
    idefix_for("LoadBuffer3D",kb.first,kb.second,jb.first,jb.second,ib.first,ib.second,
      KOKKOS_LAMBDA (int k, int j, int i) {
        out(var,k,j,i) = arr(i-ibeg + (j-jbeg)*ni + (k-kbeg)*ninj + offset );
    });

    // Update pointer
    this->pointer += ninjnk;
  }

  void Unpack(IdefixArray4D<real>& out,
       IdefixArray1D<int>& map,
       std::pair<int,int> ib,
       std::pair<int,int> jb,
       std::pair<int,int> kb) {
    const int ni = ib.second-ib.first;
    const int ninj = (jb.second-jb.first)*ni;
    const int ninjnk = (kb.second-kb.first)*ninj;
    const int ibeg = ib.first;
    const int jbeg = jb.first;
    const int kbeg = kb.first;
    const int offset = this->pointer;

    auto arr = this->array;
    idefix_for("LoadBuffer4D",0,map.size(),
                              kb.first,kb.second,
                              jb.first,jb.second,
                              ib.first,ib.second,
      KOKKOS_LAMBDA (int n, int k, int j, int i) {
        out(map(n),k,j,i) = arr(i-ibeg + (j-jbeg)*ni + (k-kbeg)*ninj + n*ninjnk + offset );
    });

    // Update pointer
    this->pointer += ninjnk*map.size();
  }


 private:
  size_t pointer;
  IdefixArray1D<real> array;
};

class Mpi {
 public:
  Mpi() = default;
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

  // Check that MPI processes are synced
  static bool CheckSync(real);


  // Destructor
  ~Mpi();

 private:
  // Because the MPI class initialise internal pointers, we do not allow copies of this class
  // These lines should not be removed as they constitute a safeguard
  Mpi(const Mpi&);
  Mpi operator=(const Mpi&);

  static int nInstances;     // total number of mpi instances in the code
  int thisInstance;          // unique number of the current instance
  int nReferences;           // # of references to this instance
  bool isInitialized{false};

  DataBlock *data;          // pointer to datablock object

  enum {faceRight, faceLeft};

  // Buffers for MPI calls
  Buffer BufferSendX1[2];
  Buffer BufferSendX2[2];
  Buffer BufferSendX3[2];
  Buffer BufferRecvX1[2];
  Buffer BufferRecvX2[2];
  Buffer BufferRecvX3[2];

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
