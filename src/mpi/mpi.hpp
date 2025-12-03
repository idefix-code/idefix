// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef MPI_MPI_HPP_
#define MPI_MPI_HPP_

#include <mpi.h>

#include <signal.h>
#include <vector>
#include <utility>
#include "idefix.hpp"
#include "grid.hpp"
#include "buffer.hpp"
#include "exchanger.hpp"


class DataBlock;


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
  enum {faceRight, faceLeft};
  // Because the MPI class initialise internal pointers, we do not allow copies of this class
  // These lines should not be removed as they constitute a safeguard
  Mpi(const Mpi&);
  Mpi operator=(const Mpi&);

  static int nInstances;     // total number of mpi instances in the code
  int thisInstance;          // unique number of the current instance
  int nReferences;           // # of references to this instance
  bool isInitialized{false};

  std::array<Exchanger,3> exchanger;  ///< exchangers in each direction
  // Error handler used by CheckConfig
  static void SigErrorHandler(int, siginfo_t* , void* );
};

#endif // MPI_MPI_HPP_
