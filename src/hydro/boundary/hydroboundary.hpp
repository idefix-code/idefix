// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2021 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef HYDRO_BOUNDARY_HYDROBOUNDARY_HPP_
#define HYDRO_BOUNDARY_HYDROBOUNDARY_HPP_
#include <string>
#include <vector>
#include "idefix.hpp"
#include "hydro_defs.hpp"
#include "grid.hpp"

#ifdef WITH_MPI
#include "mpi.hpp"
#endif

// forward class declaration
class DataBlock;
class Hydro;

class HydroBoundary {
 public:
  void Init(Input &, Grid &, Hydro* );
  void SetBoundaries(real);                         ///< Set the ghost zones in all directions
  void EnforceBoundaryDir(real, int);             ///< write in the ghost zone in specific direction
  void ReconstructVcField(IdefixArray4D<real> &);  ///< reconstruct cell-centered magnetic field
  void ReconstructNormalField(int dir);           ///< reconstruct normal field using divB=0

  void EnrollUserDefBoundary(UserDefBoundaryFunc);
  void EnrollInternalBoundary(InternalBoundaryFunc);

  void EnforcePeriodic(int, BoundarySide ); // Enforce periodic BC in direction and side

  #ifdef WITH_MPI
  Mpi mpi;                     ///< Mpi object when WITH_MPI is set
  #endif

    // User defined Boundary conditions
  UserDefBoundaryFunc userDefBoundaryFunc{NULL};
  bool haveUserDefBoundary{false};

  // Internal boundary function
  bool haveInternalBoundary{false};
  InternalBoundaryFunc internalBoundaryFunc{NULL};

  // specific for loops on ghost cells
  template <typename Function>
  void boundary_for(const std::string &,
                            const int &,
                            const BoundarySide &,
                            Function );

  template <typename Function>
  void boundary_for_all(const std::string &,
                            const int &,
                            const BoundarySide &,
                            Function );

  template <typename Function>
  void boundary_for_X1s(const std::string &,
                            const int &,
                            const BoundarySide &,
                            Function );

  template <typename Function>
  void boundary_for_X2s(const std::string &,
                            const int &,
                            const BoundarySide &,
                            Function );

  template <typename Function>
  void boundary_for_X3s(const std::string &,
                            const int &,
                            const BoundarySide &,
                            Function );

 private:
  Hydro *hydro;    // pointer to parent hydro object
  DataBlock *data;  // pointer to parent datablock
};

#endif // HYDRO_BOUNDARY_HYDROBOUNDARY_HPP_
