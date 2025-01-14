// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef UTILS_COLUMN_HPP_
#define UTILS_COLUMN_HPP_

#include <vector>
#include <string>

#include "idefix.hpp"
#include "dataBlock.hpp"



// A class to implement a parralel cumulative sum
class Column {
 public:
  ////////////////////////////////////////////////////////////////////////////////////
  /// @brief Constructor of a cumulative sum(=column density) from setup's arguments
  /// @param dir direction along which the integration is performed
  /// @param sign: +1 for an integration from left to right, -1 for an integration from right to
  ///              left i.e (backwards)
  ///////////////////////////////////////////////////////////////////////////////////
  Column(int dir, int sign, DataBlock *);

  ///////////////////////////////////////////////////////////////////////////////////
  /// @brief Effectively compute integral from the input array in argument
  /// @param in: 4D input array
  /// @param variable: index of the variable along which we do the integral (since the
  ///                   intput array is 4D)
  ///////////////////////////////////////////////////////////////////////////////////
  void ComputeColumn(IdefixArray4D<real> in, int variable);

    ///////////////////////////////////////////////////////////////////////////////////
  /// @brief Effectively compute integral from the input array in argument
  /// @param in: 3D input array
  ///////////////////////////////////////////////////////////////////////////////////
  void ComputeColumn(IdefixArray3D<real> in);

  ///////////////////////////////////////////////////////////////////////////////////
  /// @brief Get a reference to the computed column density array
  ///////////////////////////////////////////////////////////////////////////////////
  IdefixArray3D<real> GetColumn() {
    return (this->ColumnArray);
  }

 private:
  IdefixArray3D<real> ColumnArray;
  int direction; // The direction along which the column is computed
  int sign;      // whether we integrate from the left or from the right
  std::array<int,3> np_tot;
  std::array<int,3> np_int;
  std::array<int,3> beg;

  IdefixArray3D<real> Area;
  IdefixArray3D<real> Volume;

    IdefixArray2D<real> localSum;
  #ifdef WITH_MPI
  Mpi mpi;  // Mpi object when WITH_MPI is set
  MPI_Comm ColumnComm;
  int MPIrank;
  int MPIsize;

  std::array<int,3> nproc; // 3D size of the MPI cartesian geometry

  #endif
};

#endif // UTILS_COLUMN_HPP_
