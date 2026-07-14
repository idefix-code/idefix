// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef OUTPUT_XDMF_HPP_
#define OUTPUT_XDMF_HPP_
#include <string>
#if __has_include(<filesystem>)
  #include <filesystem> // NOLINT [build/c++17]
  namespace fs = std::filesystem;
#elif __has_include(<experimental/filesystem>)
  #include <experimental/filesystem>
  namespace fs = std::experimental::filesystem;
#else
  error "Missing the <filesystem> header."
#endif
#include <map>
#include "idefix.hpp"
#include "input.hpp"
#include "scalarField.hpp"

#define H5_USE_16_API
#include "hdf5.h"

#ifndef XDMF_DOUBLE
#define DUMP_DATATYPE float
#define H5_DUMP_DATATYPE H5T_NATIVE_FLOAT
#else
#define DUMP_DATATYPE double
#define H5_DUMP_DATATYPE H5T_NATIVE_DOUBLE
#endif

// Forward class declaration
class Output;
class DataBlock;

class Xdmf {
  friend class Dump;

 public:
  Xdmf(Input &, DataBlock *);   // init XDMF object
  int Write();  // Create a XDMF from the current DataBlock

  template<typename T>
  void RegisterVariable(T&, std::string, int var = -1);

 private:
  // Parent datablock
  DataBlock *data;

  // List of variables to be written to vtk files
  std::map<std::string, ScalarField> xdmfScalarMap;
  int xdmfFileNumber = 0;
  int periodicity[3];

  // dimensions
  int64_t nx1,nx2,nx3;
  int64_t nx1loc,nx2loc,nx3loc;
  int64_t nx1tot,nx2tot,nx3tot;
  int64_t nx1loctot,nx2loctot,nx3loctot;

  // number of ghost zones
  int64_t ngx1,ngx2,ngx3;

  // Coordinates needed by XDMF outputs
  DUMP_DATATYPE *xnode, *ynode, *znode;
  DUMP_DATATYPE *xcell, *ycell, *zcell;

  IdefixHostArray4D<DUMP_DATATYPE> node_coord;
  IdefixHostArray4D<DUMP_DATATYPE> cell_coord;
  // IdefixHostArray3D<DUMP_DATATYPE> field_data;

  // Array designed to store the temporary vector array
  DUMP_DATATYPE *vect3D;

  // Timer
  Kokkos::Timer timer;

  // Local offset in nodes and cells
  int nodestart[4];
  int nodesize[4];
  int nodesubsize[4];

  int cellstart[4];
  int cellsize[4];
  int cellsubsize[4];

  // output directory
  std::filesystem::path outputDirectory;

#ifdef WITH_MPI
  int mpi_data_start[3];
  int mpi_data_size[3];
  int mpi_data_subsize[3];
#endif

  void WriteHeader(
                       hid_t ,
                       const std::string ,
                       const std::string ,
                       real,
                       hid_t & ,
                       hid_t & );
  void WriteFooter(
                       const std::string ,
                       const std::string );
  void WriteScalar(
                       DUMP_DATATYPE* ,
                       const std::string &,
                       const hsize_t * ,
                       const std::string ,
                       const std::string ,
                       hid_t & ,
                       hid_t & ,
                       hid_t & ,
                       hid_t & );
};

template<typename T>
void Xdmf::RegisterVariable(T& in, std::string name, int var) {
  // if var>0, the caller provided explicitely an index
  if constexpr(std::is_same<T,IdefixArray3D<real>>::value ||
               std::is_same<T,IdefixHostArray3D<real>>::value) {
    xdmfScalarMap.emplace(name, ScalarField(in) );
  } else {
    xdmfScalarMap.emplace(name, ScalarField(in, var));
  }
}
#endif // OUTPUT_XDMF_HPP_
