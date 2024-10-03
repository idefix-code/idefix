// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef OUTPUT_OUTPUT_HPP_
#define OUTPUT_OUTPUT_HPP_
#include <string>
#include <map>
#include <vector>
#include <memory>
#include "idefix.hpp"
#include "input.hpp"
#include "dataBlock.hpp"
#include "vtk.hpp"
#ifdef WITH_HDF5
#include "xdmf.hpp"
#endif
#ifdef WITH_PYTHON
#include "pydefix.hpp"
#endif
#include "dump.hpp"
#include "slice.hpp"

using AnalysisFunc = void (*) (DataBlock &);

using UserDefVariablesContainer = std::map<std::string,IdefixHostArray3D<real>>;
using UserDefVariablesFunc = void (*) (DataBlock &, UserDefVariablesContainer &);


class Output {
  friend class Dump;    // Allow dump to have R/W access to private variables
  friend class Vtk;     // Allow VTK to have access to user-defined variables
  #ifdef WITH_HDF5
  friend class Xdmf;    // Allow XDMF to have access to user-defined variables
  #endif
  friend class DumpImage; // Allow dumpimag to have access to dump API
 public:
  Output(Input &, DataBlock &);           // Create Output Object
  int CheckForWrites(DataBlock &);        // Check if outputs are needed at this stage
  bool RestartFromDump(DataBlock &, int);  // Restart from a dump file.
  void ForceWriteDump(DataBlock &);            // Force write dumps (needed during an abort)
  void ForceWriteVtk(DataBlock &);            // Force write vtks
  #ifdef WITH_HDF5
  void ForceWriteXdmf(DataBlock &);          // Force write xdmfs
  #endif
  void ResetTimer();                      // Reset internal timer
  double GetTimer();
  void EnrollAnalysis(AnalysisFunc);
  void EnrollUserDefVariables(UserDefVariablesFunc);

 private:
  bool forceNoWrite = false;    //< explicitely disable all writes
  bool vtkEnabled = false;
  real vtkPeriod = 0.0;   // periodicity of vtk outputs
  real vtkLast = 0.0;

  bool dumpEnabled = false;
  real dumpPeriod = 0.0;
  real dumpLast = 0.0;
  real dumpTimePeriod = 0.0;
  real dumpTimeLast = 0.0;

  bool xdmfEnabled = false;
  real xdmfPeriod = 0.0;   // periodicity of xdmf outputs
  real xdmfLast = 0.0;

  bool analysisEnabled = false;
  real analysisPeriod = 0.0;
  real analysisLast = 0.0;

  bool haveAnalysisFunc = false;
  AnalysisFunc analysisFunc;

  bool userDefVariablesEnabled = false;
  bool haveUserDefVariablesFunc = false;
  UserDefVariablesFunc userDefVariablesFunc;
  UserDefVariablesContainer userDefVariables;

  bool haveSlices = false;
  std::vector<std::unique_ptr<Slice>> slices;

  #ifdef WITH_PYTHON
    Pydefix pydefix;
    bool pythonEnabled = false;
    real pythonPeriod;
    real pythonLast;
    real pythonNumber;
  #endif

  Kokkos::Timer timer;
  double elapsedTime{0.0};
};
#endif // OUTPUT_OUTPUT_HPP_
