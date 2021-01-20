// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2021 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include "output.hpp"


Output::Output(Input &input, DataBlock &data) {
  idfx::pushRegion("Output::Output");
  // initialise the output objects for each format

  // Initialise vtk outputs
  if(input.CheckEntry("Output","vtk")>0) {
    this->vtkPeriod = input.GetReal("Output","vtk",0);
    if(this->vtkPeriod>=0.0) {  // backward compatibility (negative value means no file)
      this->vtkLast = data.t - this->vtkPeriod; // write something in the next CheckForWrite()
      this->vtkEnabled = true;
      this->vtk.Init(input,data);
    }
  }


  // intialise dump outputs
  if(input.CheckEntry("Output","dmp")>0) {
    this->dumpPeriod = input.GetReal("Output","dmp",0);
    this->dumpLast = data.t - this->dumpPeriod; // dump something in the next CheckForWrite()
    this->dumpEnabled = true;
    // Backwards compatibility: negative period means no dump
    if(this->dumpPeriod<0.0) this->dumpEnabled = false;
  }
  this->dump.Init(input,data);  // Always initialised since it is required on restarts

  // initialise analysis outputs
  if(input.CheckEntry("Output","analysis")>0) {
    this->analysisPeriod = input.GetReal("Output","analysis",0);
    this->analysisLast = data.t - this->analysisPeriod; //dump something in the next CheckForWrite()
    this->analysisEnabled = true;
  }

  idfx::popRegion();
}

int Output::CheckForWrites(DataBlock &data) {
  idfx::pushRegion("Output::CheckForWrites");
  int nfiles=0;

  // Do we need a restart dump?
  if(this->dumpEnabled) {
    if(data.t >= this->dumpLast + this->dumpPeriod) {
      this->dumpLast += this->dumpPeriod;
      this->dump.Write(data,*this);
      nfiles++;
    }
  }

  // Do we need a VTK output?
  if(this->vtkEnabled) {
    if(data.t >= this->vtkLast + this->vtkPeriod) {
      this->vtkLast += this->vtkPeriod;
      this->vtk.Write(data);
      nfiles++;
    }
  }

  // Do we need an analysis ?
  if(this->analysisEnabled) {
    if(data.t >= this->analysisLast + this->analysisPeriod) {
      if(!this->haveAnalysisFunc) {
        IDEFIX_ERROR("Cannot perform a user-defined analysis without "
                     "enrollment of your analysis function");
      }
      this->analysisLast += this->analysisPeriod;
      analysisFunc(data);
      nfiles++;
    }
  }

  idfx::popRegion();

  return(nfiles);
}

void Output::RestartFromDump(DataBlock &data, int readNumber) {
  idfx::pushRegion("Output::RestartFromDump");

  this->dump.Read(data, *this, readNumber);

  idfx::popRegion();
}

void Output::ForceWrite(DataBlock &data) {
  idfx::pushRegion("Output::ForceWrite");

  this->dump.Write(data,*this);

  idfx::popRegion();
}

void Output::EnrollAnalysis(AnalysisFunc myFunc) {
  idfx::pushRegion("Output::EnrollAnalysis");
  if(!analysisEnabled) {
    IDEFIX_ERROR("You are enrolling an analysis function "
                  "but analysis are not enabled in the input file");
  }
  this->analysisFunc = myFunc;
  this->haveAnalysisFunc = true;

  idfx::popRegion();
}
