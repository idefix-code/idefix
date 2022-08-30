// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2022 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef OUTPUT_DUMP_HPP_
#define OUTPUT_DUMP_HPP_
#include <string>
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

class Dump {
  friend class DumpImage; // Allow dumpimag to have access to dump API
 public:
  void Init(Input &, DataBlock &);               // Create Dump Object
  // Create a Dump file from the current state of the code
  int Write(DataBlock &, Output&);
  // Read and load a dump file as current state of the code
  int Read(DataBlock &, Output&, int);

 private:
  int dumpFileNumber;
  int geometry{GEOMETRY};
  int periodicity[3];

  real *scrch;                            // Scratch array in host space

  // Timer
  Kokkos::Timer timer;

  // File offset
#ifdef WITH_MPI
  MPI_Offset offset;
#endif
  // These descriptors are only useful with MPI
  IdfxDataDescriptor descC;   // Descriptor for cell-centered fields (Read & write)
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
};

#endif // OUTPUT_DUMP_HPP_
