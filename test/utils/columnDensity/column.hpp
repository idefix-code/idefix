#ifndef COLUMN_HPP_
#define COLUMN_HPP_

#include "idefix.hpp"
#include "input.hpp"
#include "output.hpp"
#include "grid.hpp"
#include "dataBlock.hpp"
#include "dataBlockHost.hpp"
#include <vector>
#include <string>


// A class to implement a parralel cumulative sum
class Column {
 public:
  // Constructor from Setup arguments
  Column(int dir, int sign, int variable, DataBlock *);
  // dir : direction along which the integration is performed
  // sign: +1 for an integration from left to right, -1 for an integration from right to left (backwards)
  // variable: index of the variable along which we do the integral (since the intput array is 4D)
  void ComputeColumn(IdefixArray4D<real> in);
  // Effectively compute integral from the input array in parameter
  IdefixArray3D<real> GetColumn() {
    return (this->ColumnArray);
  }// Return the column density
 private:
  IdefixArray3D<real> ColumnArray;
  int direction; // The direction along which the column is computed
  int sign;      // whether we integrate from the left or from the right
  int variable;     // The variable we use to compute the column
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

  IdefixArray4D<real> arr4D; // Dummy array 4D used for MPI xchanges
  std::array<int,3> nproc; // 3D size of the MPI cartesian geometry

  #endif
};

#endif // COLUMN_HPP_
