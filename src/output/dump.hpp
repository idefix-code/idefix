// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef OUTPUT_DUMP_HPP_
#define OUTPUT_DUMP_HPP_
#include <string>
#include <map>
#include <array>
#if __has_include(<filesystem>)
  #include <filesystem>
  namespace fs = std::filesystem;
#elif __has_include(<experimental/filesystem>)
  #include <experimental/filesystem>
  namespace fs = std::experimental::filesystem;
#else
  error "Missing the <filesystem> header."
#endif
#include "idefix.hpp"
#include "input.hpp"
#include "dataBlock.hpp"


enum DataType {DoubleType, SingleType, IntegerType, BoolType};

// Define data descriptor used for distributed I/O when MPI is enabled
#ifdef WITH_MPI
  using IdfxDataDescriptor = MPI_Datatype;
#else
  using IdfxDataDescriptor = int;   // This is actually not used
#endif

// Forward class declaration
//class Vtk;
class Output;
class DataBlock;


class DumpField {
 public:
  enum Type {Int, Single, Double, Bool, IdefixArray};
  enum ArrayType {Device3D, Device4D, Host3D, Host4D};
  enum ArrayLocation {Center, Face, Edge};

  DumpField(IdefixArray4D<real>& in, const int varnum, const ArrayLocation loc, const int dir):
    d4Darray{in}, var{varnum}, arrayType{Device4D},
    type{IdefixArray}, arrayLocation{loc}, direction{dir} {};

  DumpField(IdefixHostArray4D<real>& in, const int varnum, const ArrayLocation loc, const int dir):
    h4Darray{in}, var{varnum}, arrayType{Host4D},
    type{IdefixArray}, arrayLocation{loc}, direction{dir} {};

  DumpField(IdefixArray3D<real>& in, const ArrayLocation loc, const int dir):
    d3Darray{in}, arrayType{Device3D},
    type{IdefixArray}, arrayLocation{loc}, direction{dir} {};

  DumpField(IdefixHostArray3D<real>& in, const ArrayLocation loc, const int dir):
    h3Darray{in}, arrayType{Host3D},
    type{IdefixArray}, arrayLocation{loc}, direction{dir} {};

  explicit DumpField(int * in, const int size = 1 ):
    rawData{static_cast<void*>(in)}, rawSize{size}, type{Int} {};

  explicit DumpField(float * in, const int size = 1 ):
    rawData{static_cast<void*>(in)}, rawSize{size}, type{Single} {};

  explicit DumpField(double * in, const int size = 1 ):
    rawData{static_cast<void*>(in)}, rawSize{size}, type{Double} {};

  explicit DumpField(bool * in, const int size = 1 ):
    rawData{static_cast<void*>(in)}, rawSize{size}, type{Bool} {};




  template <typename T>
  T GetHostField() const {
    // Check that return type is consistent with what is stored internally
    if constexpr(std::is_same<T,IdefixHostArray3D<real>>::value) {
      // User is expecting an IdefixArray, so let's convert what we have into an idefix array
      if(arrayType==Host3D) {
        return(h3Darray);
      } else if(arrayType==Host4D) {
        IdefixHostArray3D<real> arr3D = Kokkos::subview(
                                        h4Darray, var, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
        return(arr3D);
      } else if(arrayType==Device3D) {
        IdefixHostArray3D<real> arr3D = Kokkos::create_mirror(d3Darray);
        Kokkos::deep_copy(arr3D,d3Darray);
        return(arr3D);
      } else if(arrayType==Device4D) {
        IdefixArray3D<real> arrDev3D = Kokkos::subview(
                                       d4Darray, var, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
        IdefixHostArray3D<real> arr3D = Kokkos::create_mirror(arrDev3D);
        Kokkos::deep_copy(arr3D,arrDev3D);
        return(arr3D);
      } else {
        IDEFIX_ERROR("unknown field");
        return(h3Darray);
      }
    } else {
      // If it's not a raw pointer, then it must be a simple scalar
      return(static_cast<T>(rawData));
    }
  }

  // Synchronise field to Host
  template <typename T>
  void SyncFrom(T in) const {
    if constexpr(std::is_same<T,IdefixHostArray3D<real>>::value) {
      if(arrayType==Host3D) {
        Kokkos::deep_copy(h3Darray, in);
      } else if(arrayType==Host4D) {
        IdefixHostArray3D<real> arr3D = Kokkos::subview(
                                        h4Darray, var, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
        Kokkos::deep_copy(arr3D, in);
      } else if(arrayType==Device3D) {
        Kokkos::deep_copy(d3Darray,in);
      } else if(arrayType==Device4D) {
        IdefixArray3D<real> arrDev3D = Kokkos::subview(
                                       d4Darray, var, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
        Kokkos::deep_copy(arrDev3D,in);
      }
    }
    // Nothing to sync otherwise
  }

  Type GetType() const {
    return type;
  }

  ArrayType GetArrayType() const {
    return(arrayType);
  }

  ArrayLocation GetLocation() const {
    return arrayLocation;
  }

  int GetSize() const {
    return(rawSize);
  }

  int GetDirection() const {
    return(direction);
  }



 private:
  IdefixArray4D<real> d4Darray;
  IdefixArray3D<real> d3Darray;
  IdefixHostArray4D<real> h4Darray;
  IdefixHostArray3D<real> h3Darray;

  void *rawData;
  int rawSize;

  int var;
  int direction;
  ArrayLocation arrayLocation;
  ArrayType arrayType;
  Type type;
};

struct GridBox {
  std::array<int,3> start;
  std::array<int,3> size;
  std::array<int,3> sizeGlob;
};

class Dump {
  friend class DumpImage; // Allow dumpimag to have access to dump API
 public:
  explicit Dump(Input &, DataBlock *);               // Create Dump Object
  explicit Dump(DataBlock *);               // Create a dump object independent of input
  ~Dump();

  // Create a Dump file from the current state of the code
  int Write(Output&);
  // Read and load a dump file as current state of the code
  bool Read(Output&, int);

  // Register IdefixArrays
  void RegisterVariable(IdefixArray3D<real>&,
                        std::string,
                        int dir = -1,
                        DumpField::ArrayLocation loc = DumpField::ArrayLocation::Center );

  void RegisterVariable(IdefixHostArray3D<real>&,
                        std::string,
                        int dir = -1,
                        DumpField::ArrayLocation loc = DumpField::ArrayLocation::Center );

  void RegisterVariable(IdefixArray4D<real>&,
                        std::string,
                        int varnum,
                        int dir = -1,
                        DumpField::ArrayLocation loc = DumpField::ArrayLocation::Center );

  void RegisterVariable(IdefixHostArray4D<real>&,
                        std::string,
                        int varnum,
                        int dir = -1,
                        DumpField::ArrayLocation loc = DumpField::ArrayLocation::Center );

  // Register any other fundamental type
  template<typename T>
  void RegisterVariable(T*,
                        std::string,
                        int size = 1);

 private:
  void Init(DataBlock*);
  DataBlock *data;
  int dumpFileNumber;
  int geometry{GEOMETRY};
  int periodicity[3];

  real *scrch;                            // Scratch array in host space

  std::map<std::string, DumpField> dumpFieldMap;


  // Timer
  Kokkos::Timer timer;

  // File offset
#ifdef WITH_MPI
  MPI_Offset offset;
#endif
  // These descriptors are only useful with MPI
  IdfxDataDescriptor descCR;   // Descriptor for cell-centered fields (Read)
  IdfxDataDescriptor descCW;   // Descriptor for cell-centered fields (Write)
  IdfxDataDescriptor descSR[3]; // Descriptor for face-centered fields (Read)
  IdfxDataDescriptor descSW[3]; // Descriptor for face-centered fields (Write)
  IdfxDataDescriptor descER[3]; // Descriptor for edge-centered fields (Read)
  IdfxDataDescriptor descEW[3]; // Descriptor for edge-centered fields (Write)

  void WriteString(IdfxFileHandler, char *, int);
  void WriteSerial(IdfxFileHandler, int, int *, DataType, char*, void*);
  void WriteDistributed(IdfxFileHandler, int, int*, int*, char*, IdfxDataDescriptor&, real*);
  void ReadNextFieldProperties(IdfxFileHandler, int&, int*, DataType&, std::string&);
  void ReadSerial(IdfxFileHandler, int, int*, DataType, void*);
  void ReadDistributed(IdfxFileHandler, int, int*, int*, IdfxDataDescriptor&, void*);
  void Skip(IdfxFileHandler, int, int *, DataType);
  int GetLastDumpInDirectory(fs::path &);
  void CreateMPIDataType(GridBox, bool);

  fs::path outputDirectory;
};



template<typename T>
void Dump::RegisterVariable(T* in, std::string name, int size) {
    if(dumpFieldMap.count(name)>0) {
      std::stringstream msg;
      msg << "Error while registering variable in Dump I/O: " << std::endl
          << name << " has already been registered." << std::endl;
      IDEFIX_ERROR(msg);
    }
    dumpFieldMap.emplace(name, DumpField(in, size));
}



#endif // OUTPUT_DUMP_HPP_
