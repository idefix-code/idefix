// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2021 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef OUTPUT_VTK_HPP_
#define OUTPUT_VTK_HPP_
#include <string>
#include "../idefix.hpp"


class Vtk {
  friend class Dump;

 public:
  Vtk(Input &, DataBlock &);           // Create Vtk Object
  int Write(DataBlock &);     // Create a VTK from the current DataBlock

 private:
  int vtkFileNumber = 0;

  // dimensions
  int64_t nx1,nx2,nx3;
  int64_t nx1loc,nx2loc,nx3loc;

  // number of ghost zones
  int64_t ngx1,ngx2,ngx3;

  // Coordinates needed by VTK outputs
  float *node_coord, *xnode, *ynode, *znode, *Vwrite;

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
#endif

  void WriteHeader(IdfxFileHandler);
  void WriteScalar(IdfxFileHandler, float*,  std::string &);
  float BigEndian(float);
  void WriteHeaderString(const char* , IdfxFileHandler );
  void WriteHeaderFloat(float* , int64_t, IdfxFileHandler);
};

#endif // OUTPUT_VTK_HPP_
