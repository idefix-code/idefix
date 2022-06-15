// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2022 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef OUTPUT_VTK_HPP_
#define OUTPUT_VTK_HPP_
#include <string>
#include <map>
#include "idefix.hpp"
#include "input.hpp"
#include "dataBlock.hpp"

// Forward class declaration
class Output;

class Vtk {
  friend class Dump;

 public:
  void Init(Input &, DataBlock &);   // init VTK object
  int Write(DataBlock &, Output &);     // Create a VTK from the current DataBlock

 private:
  // define a mapping from global geometry flags defined in idefix.hpp
  // to the ones we write in vtk files
  std::map<int, int> VTKGeometryFlags = {
    {CARTESIAN, 0},
    {POLAR, 1},
    {SPHERICAL, 2},
    {CYLINDRICAL, 3},
  };
  int vtkFileNumber = 0;
  int geometry{VTKGeometryFlags[GEOMETRY]};
  int periodicity[3];

  // dimensions
  int64_t nx1,nx2,nx3;
  int64_t nx1loc,nx2loc,nx3loc;

  // number of ghost zones
  int64_t ngx1,ngx2,ngx3;

  // Coordinates needed by VTK outputs
  float *xnode, *ynode, *znode;

  IdefixHostArray4D<float> node_coord;

  // Array designed to store the temporary vector array
  float *vect3D;

  // Endianness swaping function and variable
  int doneEndianTest, shouldSwapEndian;

  // Timer
  Kokkos::Timer timer;

  // File offset
#ifdef WITH_MPI
  MPI_Offset offset;
  MPI_Datatype view;
  MPI_Datatype nodeView;
#endif

  void WriteHeader(IdfxFileHandler, real);
  void WriteScalar(IdfxFileHandler, float*,  const std::string &);
  template <typename T> T BigEndian(T);
  void WriteHeaderString(const char* , IdfxFileHandler );
  template <typename T> void WriteHeaderBinary(T* , int64_t, IdfxFileHandler);
  void WriteHeaderNodes(IdfxFileHandler);
};

#endif // OUTPUT_VTK_HPP_
