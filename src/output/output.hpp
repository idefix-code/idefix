// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2022 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef OUTPUT_OUTPUT_HPP_
#define OUTPUT_OUTPUT_HPP_
#include <string>
#include <map>
#include "idefix.hpp"
#include "input.hpp"
#include "dataBlock.hpp"
#include "vtk.hpp"
#include "dump.hpp"


using AnalysisFunc = void (*) (DataBlock &);

using UserDefVariablesContainer = std::map<std::string,IdefixHostArray3D<real>>;
using UserDefVariablesFunc = void (*) (DataBlock &, UserDefVariablesContainer &);


class Output {
  friend class Dump;    // Allow dump to have R/W access to private variables
  friend class Vtk;     // Allow VTK to have access to user-defined variables
  friend class DumpImage; // Allow dumpimag to have access to dump API
 public:
  Output(Input &, DataBlock &);           // Create Output Object
  int CheckForWrites(DataBlock &);        // Check if outputs are needed at this stage
  void RestartFromDump(DataBlock &, int);  // Restart from a dump file.
  void ForceWriteDump(DataBlock &);            // Force write dumps (needed during an abort)
  void ForceWriteVtk(DataBlock &);            // Force write vtks
  void ResetTimer();                      // Reset internal timer
  double GetTimer();
  void EnrollAnalysis(AnalysisFunc);
  void EnrollUserDefVariables(UserDefVariablesFunc);

 private:
  Vtk vtk;          // local instance of Vtk class
  Dump dump;        // local instance of Dump class

  bool forceNoWrite = false;    //< explicitely disable all writes
  bool vtkEnabled = false;
  real vtkPeriod = 0.0;   // periodicity of vtk outputs
  real vtkLast = 0.0;

  bool dumpEnabled = false;
  real dumpPeriod = 0.0;
  real dumpLast = 0.0;

  bool analysisEnabled = false;
  real analysisPeriod = 0.0;
  real analysisLast = 0.0;

  bool haveAnalysisFunc = false;
  AnalysisFunc analysisFunc;

  bool userDefVariablesEnabled = false;
  bool haveUserDefVariablesFunc = false;
  UserDefVariablesFunc userDefVariablesFunc;
  UserDefVariablesContainer userDefVariables;

  Kokkos::Timer timer;
  double elapsedTime{0.0};
};


#endif // OUTPUT_OUTPUT_HPP_
