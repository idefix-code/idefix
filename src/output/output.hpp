// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2021 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef OUTPUT_OUTPUT_HPP_
#define OUTPUT_OUTPUT_HPP_
#include <string>
#include "idefix.hpp"
#include "input.hpp"
#include "dataBlock.hpp"
#include "vtk.hpp"
#include "dump.hpp"


using AnalysisFunc = void (*) (DataBlock &);

using UserDefVariableFunc = void (*) (std::string &, DataBlock &, IdefixHostArray3D<real> &);

class Output {
  friend class Dump;    // Allow dump to have R/W access to private variables

 public:
  Output(Input &, DataBlock &);           // Create Output Object
  int CheckForWrites(DataBlock &);        // Check if outputs are needed at this stage
  void RestartFromDump(DataBlock &, int);  // Restart from a dump file.
  void ForceWrite(DataBlock &);            // Force write outputs (needed during an abort)

  void EnrollAnalysis(AnalysisFunc);
  void EnrollUserDefVariable(UserDefVariableFunc);

 private:
  Vtk vtk;          // local instance of Vtk class
  Dump dump;        // local instance of Dump class

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

  bool haveUserDefVariablesFunc = false;
  UserDefVariableFunc userDefVariableFunc;
};


#endif // OUTPUT_OUTPUT_HPP_
