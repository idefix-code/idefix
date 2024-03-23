#include "idefix.hpp"
#include "setup.hpp"
#include "dumpImage.hpp"

/*********************************************/
/**
Customized random number generator
Allow one to have consistant random numbers
generators on different architectures.
**/
/*********************************************/

Output* myOutput;
int outnum;
// Analysis function
// This analysis checks that the restart routines are performing as they should
void Analysis(DataBlock& data) {

  // Mirror data on Host
  idfx::cout << "Analysis: Checking dumpImage routines" << std::endl;


  // Trigger dump creation
  myOutput->ForceWriteDump(data);

  // Create a dumpImage from the created dump
  char filename[20];
  std::snprintf(filename, 20, "dump.%04d.dmp", outnum);
  #ifdef WITH_MPI
  DumpImage image(filename,&data,true); // enable domain decomposition in dumpImage
  #else
  DumpImage image(filename,&data);
  #endif

  // Check that whetever is in the image matches the current state
  DataBlockHost d(data);
  d.SyncFromDevice();

  // increment outnum
  outnum++;
  int errornum;

  errornum = 0;
  idfx::cout << "Analysis: checking dumpImage consistency with current state..." << std::endl << std::flush;
  // Check that the save/load routines have left everything unchanged.
  char fieldName[20];
  for(auto const& [name, arr] : image.arrays) {
    idfx::cout << "Array: " << name;
    idfx::cout << " ; size: " << arr.extent(0) << " " << arr.extent(1) << " " << arr.extent(2) << std::endl;

  }
  for(int n = 0; n < NVAR ; n++) {
    std::snprintf(fieldName,20,"Vc-%s",data.hydro->VcName[n].c_str());
    if(auto it = image.arrays.find(std::string(fieldName)) ; it !=  image.arrays.end()) {
      idfx::cout << "doing " << std::string(fieldName) << std::endl;
      IdefixHostArray3D<real> arr = it->second;
        for(int k = d.beg[KDIR]; k < d.end[KDIR] ; k++) {
          for(int j = d.beg[JDIR]; j < d.end[JDIR] ; j++) {
            for(int i = d.beg[IDIR]; i < d.end[IDIR] ; i++) {
            if(arr(k-d.beg[KDIR], j-d.beg[JDIR], i-d.beg[IDIR]) != d.Vc(n,k,j,i)) {
              errornum++;
              idfx::cout << "-----------------------------------------" << std::endl
                        << " Error in Vc at (i,j,k,n) = ( " << i << ", " << j << ", " << k << ", " << n << ")" << std::endl
                        << " Coordinates (x1,x2,x3)   = ( " << d.x[IDIR](i) << ", " << d.x[JDIR](j) << ", " << d.x[KDIR](k) << ")" << std::endl;
            }
          }
        }
      }
    }
  }
  std::snprintf(fieldName,20,"Vs-%s",data.hydro->VsName[0].c_str());
  for(int k = d.beg[KDIR]; k < d.end[KDIR] ; k++) {
    for(int j = d.beg[JDIR]; j < d.end[JDIR] ; j++) {
      for(int i = d.beg[IDIR]; i < d.end[IDIR]+IOFFSET ; i++) {
        if(image.arrays[fieldName](k-d.beg[KDIR], j-d.beg[JDIR], i-d.beg[IDIR]) != d.Vs(BX1s,k,j,i)) {
          errornum++;
          idfx::cout << "-----------------------------------------" << std::endl
                     << " Error in Vs(BX1s) at (i,j,k) = ( " << i << ", " << j << ", " << k << ")" << std::endl
                     << " Coordinates (x1,x2,x3)   = ( " << d.x[IDIR](i) << ", " << d.x[JDIR](j) << ", " << d.x[KDIR](k) << ")" << std::endl;
        }

      }
    }
  }
  std::snprintf(fieldName,20,"Vs-%s",data.hydro->VsName[1].c_str());
  for(int k = d.beg[KDIR]; k < d.end[KDIR] ; k++) {
    for(int j = d.beg[JDIR]; j < d.end[JDIR]+JOFFSET ; j++) {
      for(int i = d.beg[IDIR]; i < d.end[IDIR] ; i++) {
        if(image.arrays[fieldName](k-d.beg[KDIR], j-d.beg[JDIR], i-d.beg[IDIR]) != d.Vs(BX2s,k,j,i)) {
          errornum++;
          idfx::cout << "-----------------------------------------" << std::endl
                     << " Error in Vs(BX2s) at (i,j,k) = ( " << i << ", " << j << ", " << k << ")" << std::endl
                     << " Coordinates (x1,x2,x3)   = ( " << d.x[IDIR](i) << ", " << d.x[JDIR](j) << ", " << d.x[KDIR](k) << ")" << std::endl;
        }

      }
    }
  }
  std::snprintf(fieldName,20,"Vs-%s",data.hydro->VsName[2].c_str());
  for(int k = d.beg[KDIR]; k < d.end[KDIR]+KOFFSET ; k++) {
    for(int j = d.beg[JDIR]; j < d.end[JDIR] ; j++) {
      for(int i = d.beg[IDIR]; i < d.end[IDIR] ; i++) {
        if(image.arrays[fieldName](k-d.beg[KDIR], j-d.beg[JDIR], i-d.beg[IDIR]) != d.Vs(BX3s,k,j,i)) {
          errornum++;
          idfx::cout << "-----------------------------------------" << std::endl
                     << " Error in Vs(BX3s) at (i,j,k) = ( " << i << ", " << j << ", " << k << ")" << std::endl
                     << " Coordinates (x1,x2,x3)   = ( " << d.x[IDIR](i) << ", " << d.x[JDIR](j) << ", " << d.x[KDIR](k) << ")" << std::endl;
        }

      }
    }
  }


  idfx::cout << "done with " << errornum << " errors " << std::endl;
  if(errornum>0) {
    IDEFIX_ERROR("DumpImage failed validation");
  }
}

// Initialisation routine. Can be used to allocate
// Arrays or variables which are used later on
Setup::Setup(Input &input, Grid &grid, DataBlock &data, Output &output) {
   if(input.CheckEntry("Output","analysis")>0) {
     output.EnrollAnalysis(&Analysis);
     myOutput = &output;
     outnum=0;
   }
}

// This routine initialize the flow
// Note that data is on the device.
// One can therefore define locally
// a datahost and sync it, if needed
void Setup::InitFlow(DataBlock &data) {
    // Create a host copy
    DataBlockHost d(data);
    real x,y,z;

    real B0=1.0/sqrt(4.0*M_PI);

    for(int k = 0; k < d.np_tot[KDIR] ; k++) {
        for(int j = 0; j < d.np_tot[JDIR] ; j++) {
            for(int i = 0; i < d.np_tot[IDIR] ; i++) {
                x=d.x[IDIR](i);
                y=d.x[JDIR](j);
                z=d.x[KDIR](k);

                d.Vc(RHO,k,j,i) = 25.0/(36.0*M_PI);
                d.Vc(PRS,k,j,i) = 5.0/(12.0*M_PI);
                d.Vc(VX1,k,j,i) = -sin(2.0*M_PI*y);
                d.Vc(VX2,k,j,i) = sin(2.0*M_PI*x)+cos(2.0*M_PI*z);
                d.Vc(VX3,k,j,i) = cos(2.0*M_PI*x);
                #ifdef EVOLVE_VECTOR_POTENTIAL
                  real xl=d.xl[IDIR](i);
                  real yl=d.xl[JDIR](j);
                  real zl=d.xl[KDIR](k);
                  d.Ve(AX1e,k,j,i) = B0/(2.0*M_PI)*(cos(2.0*M_PI*yl));
                  d.Ve(AX2e,k,j,i) = B0/(2.0*M_PI)*sin(2.0*M_PI*xl);
                  d.Ve(AX3e,k,j,i) = B0/(2.0*M_PI)*(
                                      cos(2.0*M_PI*yl) + cos(4.0*M_PI*xl)/2.0);
                #else
                  d.Vs(BX1s,k,j,i) = -B0*sin(2.0*M_PI*y);
                  d.Vs(BX2s,k,j,i) = B0*sin(4.0*M_PI*x);
                  d.Vs(BX3s,k,j,i) = B0*(cos(2.0*M_PI*x)+sin(2.0*M_PI*y));
                #endif

            }
        }
    }

    // Send it all, if needed
    d.SyncToDevice();
}

// Analyse data to produce an output
void MakeAnalysis(DataBlock & data) {

}



// Do a specifically designed user step in the middle of the integration
void ComputeUserStep(DataBlock &data, real t, real dt) {

}
