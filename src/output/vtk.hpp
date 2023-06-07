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
#include <filesystem>
#include "idefix.hpp"
#include "input.hpp"
#include "dataBlock.hpp"

// Forward class declaration
class Output;
class DataBlock;

class VtkScalarField {
 public:
  enum Type {Device3D, Device4D, Host3D, Host4D};

  explicit VtkScalarField(IdefixArray4D<real>& in, const int varnum):
    d4Darray{in}, var{varnum}, type{Device4D} {};
  explicit VtkScalarField(IdefixHostArray4D<real>& in, const int varnum):
    h4Darray{in}, var{varnum}, type{Host4D} {};
  explicit VtkScalarField(IdefixArray3D<real>& in):
    d3Darray{in}, type{Device3D} {};
  explicit VtkScalarField(IdefixHostArray3D<real>& in):
    h3Darray{in}, type{Host3D} {};

  IdefixHostArray3D<real> GetHostField() const {
    if(type==Host3D) {
      return(h3Darray);
    } else if(type==Host4D) {
      IdefixHostArray3D<real> arr3D = Kokkos::subview(
                                      h4Darray, var, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
      return(arr3D);
    } else if(type==Device3D) {
      IdefixHostArray3D<real> arr3D = Kokkos::create_mirror(d3Darray);
      Kokkos::deep_copy(arr3D,d3Darray);
      return(arr3D);
    } else if(type==Device4D) {
      IdefixArray3D<real> arrDev3D = Kokkos::subview(
                                      d4Darray, var, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
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

class BaseVtk {
 private:
  // Endianness swaping function and variable
  bool shouldSwapEndian {true};

 protected:
  BaseVtk() {
    // Test endianness
    int tmp1 = 1;
    unsigned char *tmp2 = (unsigned char *) &tmp1;
    if (*tmp2 == 0)
      this->shouldSwapEndian = false;
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

  // Timer
  Kokkos::Timer timer;

  // DataBlock parent
  DataBlock *data;

  /* ****************************************************************************/
  /** Determines if the machine is little-endian.  If so,
    it will force the data to be big-endian.
  @param in_number floating point number to be converted in big endian */
  /* *************************************************************************** */

  template <typename T>
  T BigEndian(T in_number) {
    if (shouldSwapEndian) {
      unsigned char *bytes = (unsigned char*) &in_number;
      unsigned char tmp = bytes[0];
      bytes[0] = bytes[3];
      bytes[3] = tmp;
      tmp = bytes[1];
      bytes[1] = bytes[2];
      bytes[2] = tmp;
    }
    return(in_number);
  }

  void WriteHeaderString(const char* header, IdfxFileHandler fvtk) {
  #ifdef WITH_MPI
    MPI_Status status;
    MPI_SAFE_CALL(MPI_File_set_view(fvtk, this->offset, MPI_BYTE,
                                    MPI_CHAR, "native", MPI_INFO_NULL ));
    if(idfx::prank==0) {
      MPI_SAFE_CALL(MPI_File_write(fvtk, header, strlen(header), MPI_CHAR, &status));
    }
    offset=offset+strlen(header);
  #else
    fprintf (fvtk, "%s", header);
  #endif
  }

  template <typename T>
  void WriteHeaderBinary(T* buffer, int64_t nelem, IdfxFileHandler fvtk) {
  #ifdef WITH_MPI
    MPI_Status status;
    MPI_SAFE_CALL(MPI_File_set_view(fvtk, this->offset, MPI_BYTE, MPI_CHAR,
                                    "native", MPI_INFO_NULL ));
    if(idfx::prank==0) {
      MPI_SAFE_CALL(MPI_File_write(fvtk, buffer, nelem*sizeof(T), MPI_CHAR, &status));
    }
    offset=offset+nelem*sizeof(T);
  #else
    fwrite(buffer, sizeof(T), nelem, fvtk);
  #endif
  }
};


class Vtk : public BaseVtk {
  friend class Dump;

 public:
  Vtk(Input &, DataBlock *);   // init VTK object
  int Write();     // Create a VTK from the current DataBlock

  template<typename T>
  void RegisterVariable(T&, std::string, int var = -1);

 private:
  // List of variables to be written to vtk files
  std::map<std::string, VtkScalarField> vtkScalarMap;

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

  // File offset
#ifdef WITH_MPI
  MPI_Datatype view;
  MPI_Datatype nodeView;
#endif

  void WriteHeader(IdfxFileHandler, real);
  void WriteScalar(IdfxFileHandler, float*,  const std::string &);
  void WriteHeaderNodes(IdfxFileHandler);

  // output directory
  std::filesystem::path outputDirectory;
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
