// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include "output.hpp"
#include "dataBlock.hpp"
#include "fluid.hpp"
#include "slice.hpp"

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
    }
  }

  // Initialise xdmf outputs
  if(input.CheckEntry("Output","xdmf")>0) {
    xdmfPeriod = input.Get<real>("Output","xdmf",0);
    if(xdmfPeriod>=0.0) {  // backward compatibility (negative value means no file)
      xdmfLast = data.t - xdmfPeriod; // write something in the next CheckForWrite()
      #ifdef WITH_HDF5
      xdmfEnabled = true;
      #else
      xdmfEnabled = false;
      IDEFIX_ERROR("Attention: HDF5 library not linked when building Idefix!");
      #endif
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
      data.vtk->RegisterVariable(userDefVariables[arrayName],arrayName);
    }
    userDefVariablesEnabled = true;
  }

  // Look for slice outputs (in the form of VTK files)
  if(input.CheckEntry("Output","vtk_slice1")>0) {
    haveSlices = true;
    int n = 1;
    while(input.CheckEntry("Output","vtk_slice"+std::to_string(n))>0) {
      std::string sliceStr = "vtk_slice"+std::to_string(n);
      real period = input.Get<real>("Output", sliceStr,0);
      int direction = input.Get<int>("Output", sliceStr,1);
      real x0 = input.Get<real>("Output", sliceStr, 2);
      std::string typeStr = input.Get<std::string>("Output",sliceStr,3);
      SliceType type;
      if(typeStr.compare("cut")==0) {
        type = SliceType::Cut;
      } else if(typeStr.compare("average")==0) {
        type = SliceType::Average;
      } else {
        IDEFIX_ERROR("Unknown slice type "+typeStr);
      }
      slices.emplace_back(std::make_unique<Slice>(input, data, n, type, direction, x0, period));
      if(userDefVariablesEnabled) slices[n-1]->EnrollUserDefVariables(userDefVariables);
      // Next iteration
      n++;
    }
  }

  // Register variables that are needed in restart dumps
  data.dump->RegisterVariable(&dumpLast, "dumpLast");
  data.dump->RegisterVariable(&analysisLast, "analysisLast");
  data.dump->RegisterVariable(&vtkLast, "vtkLast");
  #ifdef WITH_HDF5
  data.dump->RegisterVariable(&xdmfLast, "xdmfLast");
  #endif

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
      data.vtk->Write();
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
  #ifdef WITH_HDF5
  // Do we need a XDMF output?
  if(xdmfEnabled) {
    if(data.t >= xdmfLast + xdmfPeriod) {
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
      xdmfLast += xdmfPeriod;
      data.xdmf->Write();
      nfiles++;
      elapsedTime += timer.seconds();

      // Check if our next predicted output should already have happened
      if((xdmfLast+xdmfPeriod <= data.t) && xdmfPeriod>0.0) {
        // Move forward xdmfLast
        while(xdmfLast <= data.t - xdmfPeriod) {
          xdmfLast += xdmfPeriod;
        }
      }
    }
  }
  #endif

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

  if(haveSlices) {
    for(int i = 0 ; i < slices.size() ; i++) {
      slices[i]->CheckForWrite(data);
    }
  }
  // Do we need a restart dump?
  if(dumpEnabled) {
    // Dumps contain metadata about the most recent outputs of other types,
    // so it's important that this part happens last.
    if(data.t >= dumpLast + dumpPeriod) {
      elapsedTime -= timer.seconds();
      dumpLast += dumpPeriod;
      data.dump->Write(*this);
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
  idfx::popRegion();

  return(nfiles);
}

bool Output::RestartFromDump(DataBlock &data, int readNumber) {
  idfx::pushRegion("Output::RestartFromDump");

  bool result = data.dump->Read(*this, readNumber);
  if(result) data.DeriveVectorPotential();

  idfx::popRegion();
  return(result);
}

void Output::ForceWriteDump(DataBlock &data) {
  idfx::pushRegion("Output::ForceWriteDump");

  if(!forceNoWrite) data.dump->Write(*this);

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
    data.vtk->Write();
    if(haveSlices) {
      for(int i = 0 ; i < slices.size() ; i++) {
        slices[i]->CheckForWrite(data,true);
      }
    }
  }
  idfx::popRegion();
}

#ifdef WITH_HDF5
void Output::ForceWriteXdmf(DataBlock &data) {
  idfx::pushRegion("Output::ForceWriteXdmf");

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
      xdmfLast += xdmfPeriod;
      data.xdmf->Write();
  }
  idfx::popRegion();
}
#endif

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
  if(haveSlices) {
    for(int i = 0 ; i < slices.size() ; i++) {
      slices[i]->EnrollUserDefFunc(myFunc);
    }
  }
  haveUserDefVariablesFunc = true;
  idfx::popRegion();
}

void Output::ResetTimer() {
  elapsedTime = 0.0;
}

double Output::GetTimer() {
  return(elapsedTime);
}
