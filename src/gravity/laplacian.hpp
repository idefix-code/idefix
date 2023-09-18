// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef GRAVITY_LAPLACIAN_HPP_
#define GRAVITY_LAPLACIAN_HPP_

#include "idefix.hpp"


class Laplacian {
 public:
  Lapacian(DataBlock &, LaplacianBoundaryType &, LaplacianBoundaryType &, bool );

  // Types of boundary which can be treated
  enum LaplacianBoundaryType {internalgrav, periodic, nullgrad, nullpot, userdef, axis, origin};
  void InitPreconditionner();   // For preconditionning versions
  void PreComputeLaplacian();   // For faster Laplacian computation

  void InitInternalGrid(); // initialise the extra internal grid (for origin BCs)

  void SetBoundaries(IdefixArray3D<real> &);  // Set the proper boundaries for the given array

  void EnrollUserDefBoundary(UserDefBoundaryFunc);  // Enroll user-defined boundary conditions

  // The main laplacian operator
  void operator() (IdefixArray3D<real> in,  IdefixArray3D<real> laplacian);

  // Handling userdef boundary.
  using UserDefBoundaryFunc = void (*) (DataBlock &, int dir, BoundarySide side,
                                       const real t, IdefixArray3D<real> &arr);
  void EnforceBoundary(int dir, BoundarySide side, GravityBoundaryType type,
                       IdefixArray3D<real> &);
  bool haveUserDefBoundary{false};
  // User defined Boundary conditions
  UserDefBoundaryFunc userDefBoundaryFunc{NULL};

  // Local potential array size
  std::vector<int> np_tot;
  std::vector<int> np_int;

  std::vector<int> nghost;
  std::vector<int> beg;
  std::vector<int> end;

  // offset in the left and right direction between selfgravity grid and datablock grid
  std::vector<int> loffset;
  std::vector<int> roffset;

  // Grid for self-gravity solver (note that this grid may extend the grid of the current datablock)
  std::vector<IdefixArray1D<real>> x;    ///< geometrical central points
  std::vector<IdefixArray1D<real>> dx;   ///< cell width

  IdefixArray1D<real> sinx2;            ///< sinx2 (only in spherical)

  IdefixArray3D<real> dV;                ///< cell volume
  std::vector<IdefixArray3D<real>> A;    ///< cell left interface area

  bool isPeriodic; // Periodicity status of the density distribution
  std::array<LaplacianBoundaryType,3> lbound;  // Boundary condition to the left
  std::array<LaplacianBoundaryType,3> rbound;  // Boundary condition to the right
                           // Warning : might differ from (M)HD solver !

  IdefixArray3D<real> precond; // Diagonal preconditionner
  IdefixArray4D<real> Lx1; //< Laplacian operator in x1
  IdefixArray4D<real> Lx2; //< Laplacian operator in x2
  IdefixArray4D<real> Lx3; //< Laplacian operator in x3

  bool isTwoPi{false};
  bool havePreconditioner{false}; // Use of preconditionner (or not)


  #ifdef WITH_MPI
  Mpi mpi;  // Mpi object when WITH_MPI is set
  IdefixArray4D<real> arr4D; // Intermediate array for boundary handling
  #endif
};


#endif