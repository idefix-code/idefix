// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef OUTPUT_OUTPUTDUMP_HPP_
#define OUTPUT_OUTPUTDUMP_HPP_
#include <string>
#include "../idefix.hpp"


enum DataType {DoubleType, SingleType, IntegerType};

// Define data descriptor used for distributed I/O when MPI is enabled
#ifdef WITH_MPI
  using IdfxDataDescriptor = MPI_Datatype;
#else
  using IdfxDataDescriptor = int;   // This is actually not used
#endif

// Forward class declaration
class OutputVTK;

class OutputDump {
 public:
  OutputDump(Input &, DataBlock &, real);               // Create Output Object
  // Create a Dump file from the current state of the code
  int Write(Grid&, DataBlock &, TimeIntegrator&, OutputVTK&);
  // Create a Dump file from the current state of the code
  int CheckForWrite(Grid&, DataBlock &, TimeIntegrator&, OutputVTK&);
  // Read and load a dump file as current state of the code
  int Read(Grid&, DataBlock &, TimeIntegrator&, OutputVTK&, int);

 private:
  int dumpFileNumber;
  real tperiod, tnext;

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



  void WriteString(IdfxFileHandler, char *);
  void WriteSerial(IdfxFileHandler, int, int *, DataType, char*, void*);
  void WriteDistributed(IdfxFileHandler, int, int*, int*, char*, IdfxDataDescriptor&, real*);
  void ReadNextFieldProperties(IdfxFileHandler, int&, int*, DataType&, std::string&);
  void ReadSerial(IdfxFileHandler, int, int*, DataType, void*);
  void ReadDistributed(IdfxFileHandler, int, int*, int*, IdfxDataDescriptor&, void*);
};

#endif // OUTPUT_OUTPUTDUMP_HPP_
