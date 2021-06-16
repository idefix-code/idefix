// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2021 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef DATABLOCK_DATABLOCK_HPP_
#define DATABLOCK_DATABLOCK_HPP_

#include <vector>
#include <string>
#include "idefix.hpp"
#include "grid.hpp"
#include "gridHost.hpp"
#include "hydro.hpp"

//TODO(lesurg) What is this standing for?
#define BOUNDARY_

class DataBlock {
 public:
  // Local grid information
  IdefixArray1D<real> x[3];    ///< geometrical central points
  IdefixArray1D<real> xr[3];   ///< cell right interface
  IdefixArray1D<real> xl[3];   ///< cell left interface
  IdefixArray1D<real> dx[3];   ///< cell width
  IdefixArray1D<real> xgc[3];  ///< cell geometrical cell center
  IdefixArray1D<real> rt;      ///< In spherical coordinates, gives \tilde{r}
  IdefixArray1D<real> sinx2m;  ///< In spherical coordinates,
                               ///< gives sin(th) at a j-1/2 interface
  IdefixArray1D<real> tanx2m;  ///< In spherical coordinates,
                               ///< gives tan(th) at a j-1/2 interface
  IdefixArray1D<real> sinx2;   ///< In spherical coordinates, gives sin(th) at the cell center
  IdefixArray1D<real> tanx2;   ///< In spherical coordinates, gives tan(th) at the cell center
  IdefixArray1D<real> dmu;     ///< In spherical coordinates,
                               ///< gives the \theta volume = fabs(cos(th_m) - cos(th_p))

  real xbeg[3];                ///< Beginning of active domain in datablock
  real xend[3];                ///< End of active domain in datablock

  IdefixArray3D<real> dV;      ///< cell volume
  IdefixArray3D<real> A[3];    ///< cell right interface area

  int np_tot[3];               ///< total number of grid points in datablock
  int np_int[3];               ///< active number of grid points in datablock (excl. ghost cells)

  int nghost[3];               ///< number of ghost cells at each boundary
  BoundaryType lbound[3];      ///< Boundary condition to the left
  BoundaryType rbound[3];      ///< Boundary condition to the right

  bool haveAxis = false;       ///< DataBlock contains points on the axis and a special treatment
                               ///< has been required for these.

  int beg[3];                  ///< First local index of the active domain
  int end[3];                  ///< Last local index of the active domain+1

  int gbeg[3];                 ///< First global index of the active domain of this datablock
  int gend[3];                 ///< Last global index of the active domain of this datablock

  real dt;                     ///< Current timestep
  real t;                      ///< Current time

  Grid *mygrid;                ///< Parent grid object

  Hydro hydro;                  ///< The Hydro object attached to this datablock

  void InitFromGrid(Grid &, Input &); ///< init from a Grid object
  void MakeGeometry();                ///< Compute geometrical terms
  void DumpToFile(std::string);   ///< Dump current datablock to a file for inspection
  int CheckNan();                 ///< Return the number of cells which have Nans

  bool rklCycle{false};           ///<  // Set to true when we're inside a RKL call

  void EvolveStage();             ///< Evolve this DataBlock by dt
  void SetBoundaries();       ///< Enforce boundary conditions to this datablock


  void ResetStage();              ///< Reset the variables needed at each major integration Stage

  DataBlock();

 private:
  void WriteVariable(FILE* , int , int *, char *, void*);

  template<int dir> void LoopDir();     ///< // recursive loop on dimensions
};

#endif // DATABLOCK_DATABLOCK_HPP_
