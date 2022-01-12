// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2022 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
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

  void EnforceFluxBoundaries(int dir);        ///< Apply boundary condition conditions to the fluxes

  void EnrollUserDefBoundary(UserDefBoundaryFunc);
  void EnrollInternalBoundary(InternalBoundaryFunc);
  void EnrollFluxBoundary(UserDefBoundaryFunc);

  void EnforcePeriodic(int, BoundarySide ); ///< Enforce periodic BC in direction and side
  void EnforceReflective(int, BoundarySide ); ///< Enforce reflective BC in direction and side
  void EnforceOutflow(int, BoundarySide ); ///< Enforce outflow BC in direction and side
  void EnforceShearingBox(real, int, BoundarySide ); ///< Enforce Shearing box BCs

  #ifdef WITH_MPI
  Mpi mpi;                     ///< Mpi object when WITH_MPI is set
  #endif

    // User defined Boundary conditions
  UserDefBoundaryFunc userDefBoundaryFunc{NULL};
  bool haveUserDefBoundary{false};

  // Internal boundary function
  bool haveInternalBoundary{false};
  InternalBoundaryFunc internalBoundaryFunc{NULL};

  // Flux boundary function
  bool haveFluxBoundary{false};
  UserDefBoundaryFunc fluxBoundaryFunc{NULL};

  // specific for loops on ghost cells
  template <typename Function>
  void BoundaryFor(const std::string &,
                            const int &,
                            const BoundarySide &,
                            Function );

  template <typename Function>
  void BoundaryForAll(const std::string &,
                            const int &,
                            const BoundarySide &,
                            Function );

  template <typename Function>
  void BoundaryForX1s(const std::string &,
                            const int &,
                            const BoundarySide &,
                            Function );

  template <typename Function>
  void BoundaryForX2s(const std::string &,
                            const int &,
                            const BoundarySide &,
                            Function );

  template <typename Function>
  void BoundaryForX3s(const std::string &,
                            const int &,
                            const BoundarySide &,
                            Function );
  IdefixArray4D<real> sBArray;    ///< Array use by shearingbox boundary conditions

 private:
  Hydro *hydro;    // pointer to parent hydro object
  DataBlock *data;  // pointer to parent datablock
};

#endif // HYDRO_BOUNDARY_HYDROBOUNDARY_HPP_
