// ***********************************************************************************
// Idefix MHD astrophysical code
//
// Source file src/output/output.hpp
//
// Last modified : 10/2024
//
// Copyright(C) by :
// - Geoffroy Lesur <geoffroy.lesur@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2023 - 2024)
// - Gaylor Wafflard <gaylor.wafflard@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2024)
// and other code contributors
//
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef OUTPUT_OUTPUT_HPP_
#define OUTPUT_OUTPUT_HPP_

#include <map>
#include <vector>
#include <string>
#include <memory>
#include <functional>
#include "idefix.hpp"
#include "dataBlock.hpp"
#include "input.hpp"
#include "vtk.hpp"
#include "dump.hpp"
#ifdef WITH_HDF5
#include "xdmf.hpp"
#endif
#include "slice.hpp"
#ifdef WITH_PYTHON
  #include "pydefix.hpp"
#endif

using AnalysisFunc = std::function<void(DataBlock &)>;

using UserDefVariablesContainer = std::map<std::string,IdefixHostArray3D<real>>;
using UserDefVariablesFunc = std::function<void(DataBlock &, UserDefVariablesContainer &)>;


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

  // Template enrollment methods (definitions inline)
  template<typename Func>
  void EnrollAnalysis(Func myFunc) {
    idfx::pushRegion("Output::EnrollAnalysis");
    if(!analysisEnabled) {
      IDEFIX_WARNING("You are enrolling an analysis function "
                    "but analysis are not enabled in the input file,"
                    "so no analysis will be performed");
    }
    analysisFunc = myFunc;
    haveAnalysisFunc = true;
    idfx::popRegion();
  }

  template<typename Func>
  void EnrollUserDefVariables(Func myFunc) {
    idfx::pushRegion("Output::EnrollUserDefVariable");
    if(!userDefVariablesEnabled) {
      IDEFIX_ERROR("You are enrolling a user-defined variables function "
                   "but the userdef variables are not set in the input file");
    }
    userDefVariablesFunc = myFunc;
    if(haveSlices) {
      for(int i = 0 ; i < slices.size() ; i++) {
        slices[i]->EnrollUserDefFunc(myFunc);
      }
    }
    haveUserDefVariablesFunc = true;
    idfx::popRegion();
  }

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
