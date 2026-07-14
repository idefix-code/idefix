// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef OUTPUT_VTK_HPP_
#define OUTPUT_VTK_HPP_
#include <string>
#include <map>
#if __has_include(<filesystem>)
  #include <filesystem> // NOLINT [build/c++17]
  namespace fs = std::filesystem;
#elif __has_include(<experimental/filesystem>)
  #include <experimental/filesystem>
  namespace fs = std::experimental::filesystem;
#else
  error "Missing the <filesystem> header."
#endif
#include <cstdio>
#include "idefix.hpp"
#include "input.hpp"
#include "dataBlock.hpp"
#include "bigEndian.hpp"
#include "scalarField.hpp"


// Forward class declaration
class Output;
class DataBlock;


class BaseVtk {
 protected:
  BaseVtk() {
    // Initialise the root tag (used for MPI non-collective I/Os)
    this->isRoot = idfx::prank == 0;
  }

  // File offset
#ifdef WITH_MPI
  MPI_Offset offset;
#endif
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

  bool isRoot;  // Whether this process is our root process for the current
  // Timer
  Kokkos::Timer timer;

  // DataBlock parent
  DataBlock *data;

  // BigEndian conversion
  BigEndian bigEndian;


  void WriteHeaderString(const char* header, IdfxFileHandler fvtk) {
  #ifdef WITH_MPI
    MPI_Status status;
    MPI_SAFE_CALL(MPI_File_set_view(fvtk, this->offset, MPI_BYTE,
                                    MPI_CHAR, "native", MPI_INFO_NULL ));
    if(this->isRoot) {
      MPI_SAFE_CALL(MPI_File_write(fvtk, header, strlen(header), MPI_CHAR, &status));
    }
    offset=offset+strlen(header);
  #else
    int rc = fprintf (fvtk, "%s", header);
    if(rc<0) {
      IDEFIX_ERROR("Unable to write to file. Check your filesystem permissions and disk quota.");
    }
  #endif
  }

  template <typename T>
  void WriteHeaderBinary(T* buffer, int64_t nelem, IdfxFileHandler fvtk) {
  #ifdef WITH_MPI
    MPI_Status status;
    MPI_SAFE_CALL(MPI_File_set_view(fvtk, this->offset, MPI_BYTE, MPI_CHAR,
                                    "native", MPI_INFO_NULL ));
    if(this->isRoot) {
      MPI_SAFE_CALL(MPI_File_write(fvtk, buffer, nelem*sizeof(T), MPI_CHAR, &status));
    }
    offset=offset+nelem*sizeof(T);
  #else
    if(fwrite(buffer, sizeof(T), nelem, fvtk) != nelem) {
      IDEFIX_ERROR("Unable to write to file. Check your filesystem permissions and disk quota.");
    }
  #endif
  }
};


class Vtk : public BaseVtk {
  friend class Dump;
  friend class Slice;

 public:
  explicit Vtk(Input &, DataBlock *, std::string filebase = "data");   // init VTK object
  int Write();     // Create a VTK from the current DataBlock

  template<typename T>
  void RegisterVariable(T&, std::string, int var = -1);

 private:
  // List of variables to be written to vtk files
  std::map<std::string, ScalarField> vtkScalarMap;

  // dimensions
  int64_t nx1,nx2,nx3;
  int64_t nx1loc,nx2loc,nx3loc;

  // number of ghost zones
  int64_t ngx1,ngx2,ngx3;

  // offset in each direction (used by the vtk grid)
  int ioffset, joffset, koffset;

  // Coordinates needed by VTK outputs
  float *xnode, *ynode, *znode;
  float *xcenter, *ycenter, *zcenter;

  IdefixHostArray4D<float> node_coord;

  // Array designed to store the temporary vector array
  float *vect3D;

  // File name
  std::string filebase;

  // File offset
#ifdef WITH_MPI
  MPI_Datatype view;
  MPI_Datatype nodeView;
  MPI_Comm comm;
#endif

  void WriteHeader(IdfxFileHandler, real);
  void WriteScalar(IdfxFileHandler, float*,  const std::string &);
  void WriteHeaderNodes(IdfxFileHandler);

  // output directory
  fs::path outputDirectory;
};

template<typename T>
void Vtk::RegisterVariable(T& in, std::string name, int var) {
  // if var>0, the caller provided explicitely an index
  if constexpr(std::is_same<T,IdefixArray3D<real>>::value ||
               std::is_same<T,IdefixHostArray3D<real>>::value) {
    vtkScalarMap.emplace(name, ScalarField(in) );
  } else {
    vtkScalarMap.emplace(name, ScalarField(in, var));
  }
}

#endif // OUTPUT_VTK_HPP_
