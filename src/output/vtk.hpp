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
#include "idefix.hpp"
#include "input.hpp"
#include "dataBlock.hpp"

// Forward class declaration
class Output;
class DataBlock;

class VtkScalarField {
 public:
  enum Type {Device3D, Device4D, Host3D, Host4D};

  VtkScalarField(IdefixArray4D<real>& in, const int varnum): d4Darray{in}, var{varnum}, type{Device4D} {};
  VtkScalarField(IdefixHostArray4D<real>& in, const int varnum): h4Darray{in}, var{varnum}, type{Host4D} {};
  VtkScalarField(IdefixArray3D<real>& in): d3Darray{in}, type{Device3D} {};
  VtkScalarField(IdefixHostArray3D<real>& in): h3Darray{in}, type{Host3D} {};

  IdefixHostArray3D<real> GetHostField() const {
    if(type==Host3D) {
      return(h3Darray);
    } else if(type==Host4D) {
      IdefixHostArray3D<real> arr3D = Kokkos::subview(h4Darray, var, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
      return(arr3D);
    } else if(type==Device3D) {
      IdefixHostArray3D<real> arr3D = Kokkos::create_mirror(d3Darray);
      Kokkos::deep_copy(arr3D,d3Darray);
      return(arr3D);
    } else if(type==Device4D) {
      IdefixArray3D<real> arrDev3D = Kokkos::subview(d4Darray, var, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
      IdefixHostArray3D<real> arr3D = Kokkos::create_mirror(arrDev3D);
      Kokkos::deep_copy(arr3D,arrDev3D);
      return(arr3D);
    } else {
      IDEFIX_ERROR("unknown field");
      return(h3Darray);
    }
  }

 private:
  IdefixArray4D<real> d4Darray;
  IdefixArray3D<real> d3Darray;
  IdefixHostArray4D<real> h4Darray;
  IdefixHostArray3D<real> h3Darray;
  int var;
  Type type;

};


class Vtk {
  friend class Dump;

 public:
  Vtk(Input &, DataBlock *);   // init VTK object
  int Write();     // Create a VTK from the current DataBlock

  template<typename T>
  void RegisterVariable(T&, std::string, int var = -1);

 private:
  // define a mapping from global geometry flags defined in idefix.hpp
  // to the ones we write in vtk files
  std::map<int, int> VTKGeometryFlags = {
    {CARTESIAN, 0},
    {POLAR, 1},
    {SPHERICAL, 2},
    {CYLINDRICAL, 3},
  };

  // List of variables to be written to vtk files
  std::map<std::string, VtkScalarField> vtkScalarMap;

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

  // DataBlock parent
  DataBlock *data;

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

template<typename T>
void Vtk::RegisterVariable(T& in, std::string name, int var) {
  // if var>0, the caller provided explicitely an index
  if constexpr(std::is_same<T,IdefixArray3D<real>>::value ||
               std::is_same<T,IdefixHostArray3D<real>>::value) {
    vtkScalarMap.emplace(name, VtkScalarField(in) );
  } else {
    vtkScalarMap.emplace(name, VtkScalarField(in, var));
  }
}

#endif // OUTPUT_VTK_HPP_
