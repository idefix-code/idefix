// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2021 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef OUTPUT_OUTPUT_HPP_
#define OUTPUT_OUTPUT_HPP_
#include <string>
#include "../idefix.hpp"


class Output {
 public:
  Output(Input &, DataBlock &);           // Create Vtk Object
  int CheckForWrites(DataBlock &);        // Check if outputs are needed at this stage
  int RestartFromDump(DataBlock &, int);  // Restart from a dump file.
  int ForceWrite(DataBlock &);            // Force write outputs (needed during an abort)

  Vtk vtk;          // local instance of Vtk class
  Dump dump;        // local instance of Dump class
 private:
  real vtkPeriod;   // periodicity of vtk outputs
  real vtkLast;

  real dumpPeriod;
  real dumpLast;
}


#endif // OUTPUT_OUTPUT_HPP_
