// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2022 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include "output.hpp"


Output::Output(Input &input, DataBlock &data) {
  idfx::pushRegion("Output::Output");
  // initialise the output objects for each format

  if(input.forceNoWrite) {
    this->forceNoWrite = true;
  }
  // Initialise vtk outputs
  if(input.CheckEntry("Output","vtk")>0) {
    vtkPeriod = input.Get<real>("Output","vtk",0);
    if(vtkPeriod>=0.0) {  // backward compatibility (negative value means no file)
      vtkLast = data.t - vtkPeriod; // write something in the next CheckForWrite()
      vtkEnabled = true;
      vtk.Init(input,data);
    }
  }

  // intialise dump outputs
  if(input.CheckEntry("Output","dmp")>0) {
    dumpPeriod = input.Get<real>("Output","dmp",0);
    dumpLast = data.t - dumpPeriod; // dump something in the next CheckForWrite()
    dumpEnabled = true;
    // Backwards compatibility: negative period means no dump
    if(dumpPeriod<0.0) dumpEnabled = false;
  }
  dump.Init(input,data);  // Always initialised since it is required on restarts

  // initialise analysis outputs
  if(input.CheckEntry("Output","analysis")>0) {
    analysisPeriod = input.Get<real>("Output","analysis",0);
    analysisLast = data.t - analysisPeriod; //dump something in the next CheckForWrite()
    analysisEnabled = true;
  }

  // Initialise userdefined outputs
  if(input.CheckEntry("Output","uservar")>0) {
    int nvars = input.CheckEntry("Output","uservar");
    for(int var = 0 ; var < nvars ; var++) {
      std::string arrayName = input.Get<std::string>("Output","uservar",var);
      // Create an array to store this user variable
      // and store the whole thing in the container
      userDefVariables[arrayName] = IdefixHostArray3D<real>("userVar-"+arrayName,
                                                                  data.np_tot[KDIR],
                                                                  data.np_tot[JDIR],
                                                                  data.np_tot[IDIR]);
    }
    userDefVariablesEnabled = true;
  }
  idfx::popRegion();
}

int Output::CheckForWrites(DataBlock &data) {
  idfx::pushRegion("Output::CheckForWrites");
  int nfiles=0;

  // skip everything if forced to disable all writes
  if(forceNoWrite) {
    idfx::popRegion();
    return(0);
  }
  // Do we need a restart dump?
  if(dumpEnabled) {
    if(data.t >= dumpLast + dumpPeriod) {
      elapsedTime -= timer.seconds();
      dumpLast += dumpPeriod;
      dump.Write(data,*this);
      nfiles++;
      elapsedTime += timer.seconds();

      // Check if our next predicted output should already have happened
      if((dumpLast+dumpPeriod <= data.t) && dumpPeriod>0.0) {
        // Move forward dumpLast
        while(dumpLast <= data.t - dumpPeriod) {
          dumpLast += dumpPeriod;
        }
      }
    }
  }

  // Do we need a VTK output?
  if(vtkEnabled) {
    if(data.t >= vtkLast + vtkPeriod) {
      elapsedTime -= timer.seconds();
      if(userDefVariablesEnabled) {
        if(haveUserDefVariablesFunc) {
          // Call user-def function to fill the userdefined variable arrays
          idfx::pushRegion("UserDef::User-defined variables function");
          userDefVariablesFunc(data, userDefVariables);
          idfx::popRegion();
        } else {
          IDEFIX_ERROR("Cannot output user-defined variables without "
                        "enrollment of your user-defined variables function");
        }
      }
      vtkLast += vtkPeriod;
      vtk.Write(data, *this);
      nfiles++;
      elapsedTime += timer.seconds();

      // Check if our next predicted output should already have happened
      if((vtkLast+vtkPeriod <= data.t) && vtkPeriod>0.0) {
        // Move forward vtkLast
        while(vtkLast <= data.t - vtkPeriod) {
          vtkLast += vtkPeriod;
        }
      }
    }
  }

  // Do we need an analysis ?
  if(analysisEnabled) {
    if(data.t >= analysisLast + analysisPeriod) {
      elapsedTime -= timer.seconds();
      if(!haveAnalysisFunc) {
        IDEFIX_ERROR("Cannot perform a user-defined analysis without "
                     "enrollment of your analysis function");
      }
      analysisLast += analysisPeriod;
      idfx::pushRegion("UserDef::User-defined analysis function");
      analysisFunc(data);
      idfx::popRegion();
      nfiles++;
      elapsedTime += timer.seconds();

      // Check if our next predicted output should already have happened
      if((analysisLast+analysisPeriod <= data.t) && analysisPeriod>0.0) {
        // Move forward analysisLast
        while(analysisLast <= data.t - analysisPeriod) {
          analysisLast += analysisPeriod;
        }
      }
    }
  }

  idfx::popRegion();

  return(nfiles);
}

void Output::RestartFromDump(DataBlock &data, int readNumber) {
  idfx::pushRegion("Output::RestartFromDump");

  dump.Read(data, *this, readNumber);

  idfx::popRegion();
}

void Output::ForceWriteDump(DataBlock &data) {
  idfx::pushRegion("Output::ForceWriteDump");

  if(!forceNoWrite) dump.Write(data,*this);

  idfx::popRegion();
}

void Output::ForceWriteVtk(DataBlock &data) {
  idfx::pushRegion("Output::ForceWriteVtk");

  if(!forceNoWrite) {
    if(userDefVariablesEnabled) {
        if(haveUserDefVariablesFunc) {
          // Call user-def function to fill the userdefined variable arrays
          idfx::pushRegion("UserDef::User-defined variables function");
          userDefVariablesFunc(data, userDefVariables);
          idfx::popRegion();
        } else {
          IDEFIX_ERROR("Cannot output user-defined variables without "
                        "enrollment of your user-defined variables function");
        }
      }
      vtkLast += vtkPeriod;
      vtk.Write(data, *this);
  }
  idfx::popRegion();
}

void Output::EnrollAnalysis(AnalysisFunc myFunc) {
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

void Output::EnrollUserDefVariables(UserDefVariablesFunc myFunc) {
  idfx::pushRegion("Output::EnrollUserDefVariable");
  if(!userDefVariablesEnabled) {
    IDEFIX_ERROR("You are enrolling a user-defined variables function "
                 "but the userdef variables are not set in the input file");
  }
  userDefVariablesFunc = myFunc;
  haveUserDefVariablesFunc = true;
  idfx::popRegion();
}

void Output::ResetTimer() {
  elapsedTime = 0.0;
}

double Output::GetTimer() {
  return(elapsedTime);
}
