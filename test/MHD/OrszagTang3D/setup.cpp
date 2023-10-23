#include "idefix.hpp"
#include "setup.hpp"

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


  idfx::cout << "Analysis: Checking restart routines" << std::endl;

  // Trigger dump creation
  myOutput->ForceWriteDump(data);

    // Mirror data on Host
  DataBlockHost d(data);

  // Sync it
  d.SyncFromDevice();

  // Create local arrays to store the current physical state
  IdefixHostArray4D<real> myVc = IdefixHostArray4D<real>("myVc", d.Vc.extent(0), data.np_tot[KDIR], data.np_tot[JDIR],data.np_tot[IDIR]);
  IdefixHostArray4D<real> myVs = IdefixHostArray4D<real>("myVs", DIMENSIONS, data.np_tot[KDIR]+KOFFSET, data.np_tot[JDIR]+JOFFSET,data.np_tot[IDIR]+IOFFSET);
  #ifdef EVOLVE_VECTOR_POTENTIAL
  IdefixHostArray4D<real> myVe = IdefixHostArray4D<real>("myVe", AX3e+1, data.np_tot[KDIR]+KOFFSET, data.np_tot[JDIR]+JOFFSET,data.np_tot[IDIR]+IOFFSET);
  #endif
  // Transfer the datablock to myVc and myVs
  for(int n = 0; n < d.Vc.extent(0) ; n++) {
    for(int k = 0; k < d.np_tot[KDIR] ; k++) {
      for(int j = 0; j < d.np_tot[JDIR] ; j++) {
        for(int i = 0; i < d.np_tot[IDIR] ; i++) {
          myVc(n,k,j,i) = d.Vc(n,k,j,i);
          d.Vc(n,k,j,i) = 0.0;

        }
      }
    }
  }

  for(int n = 0; n < DIMENSIONS ; n++) {
    for(int k = 0; k < d.np_tot[KDIR] + KOFFSET; k++) {
      for(int j = 0; j < d.np_tot[JDIR] + JOFFSET; j++) {
        for(int i = 0; i < d.np_tot[IDIR] + IOFFSET; i++) {
          myVs(n,k,j,i) = d.Vs(n,k,j,i);
          d.Vs(n,k,j,i) = 0.0;
        }
      }
    }
  }
  #ifdef EVOLVE_VECTOR_POTENTIAL
    for(int n = 0; n < AX3e+1 ; n++) {
      for(int k = 0; k < d.np_tot[KDIR] + KOFFSET; k++) {
        for(int j = 0; j < d.np_tot[JDIR] + JOFFSET; j++) {
          for(int i = 0; i < d.np_tot[IDIR] + IOFFSET; i++) {
            myVe(n,k,j,i) = d.Ve(n,k,j,i);
            d.Ve(n,k,j,i) = 0.0;
          }
        }
      }
    }
  #endif

  // Push our datablockHost to erase everything
  d.SyncToDevice();
  // From this point, the dataBlock is full of zeros

  // Load back the restart dump
  myOutput->RestartFromDump(data, outnum);
  data.SetBoundaries();
  #ifdef EVOLVE_VECTOR_POTENTIAL
    data.hydro->emf->ComputeMagFieldFromA(data.hydro->Ve, data.hydro->Vs);
  #endif
  d.SyncFromDevice();

  // increment outnum
  outnum++;
  int errornum;

  errornum = 0;
  idfx::cout << "Analysis: checking consistency" << std::endl;
  // Check that the save/load routines have left everything unchanged.
  for(int n = 0; n < d.Vc.extent(0) ; n++) {
    for(int k = d.beg[KDIR]; k < d.end[KDIR] ; k++) {
      for(int j = d.beg[JDIR]; j < d.end[JDIR] ; j++) {
        for(int i = d.beg[IDIR]; i < d.end[IDIR] ; i++) {
          if(myVc(n,k,j,i) != d.Vc(n,k,j,i)) {
            errornum++;
            idfx::cout << "-----------------------------------------" << std::endl
                       << " Error in Vc at (i,j,k,n) = ( " << i << ", " << j << ", " << k << ", " << n << ")" << std::endl
                       << " Coordinates (x1,x2,x3)   = ( " << d.x[IDIR](i) << ", " << d.x[JDIR](j) << ", " << d.x[KDIR](k) << ")" << std::endl;
          }

        }
      }
    }
  }

    for(int k = d.beg[KDIR]; k < d.end[KDIR] ; k++) {
      for(int j = d.beg[JDIR]; j < d.end[JDIR] ; j++) {
        for(int i = d.beg[IDIR]; i < d.end[IDIR]+IOFFSET ; i++) {
          if(myVs(BX1s,k,j,i) != d.Vs(BX1s,k,j,i)) {
            errornum++;
            idfx::cout << "-----------------------------------------" << std::endl
                       << " Error in Vs(BX1s) at (i,j,k) = ( " << i << ", " << j << ", " << k << ")" << std::endl
                       << " Coordinates (x1,x2,x3)   = ( " << d.x[IDIR](i) << ", " << d.x[JDIR](j) << ", " << d.x[KDIR](k) << ")" << std::endl
                       << "Original= " << myVs(BX1s,k,j,i) << " New=" << d.Vs(BX1s,k,j,i) << " diff=" << myVs(BX1s,k,j,i)-d.Vs(BX1s,k,j,i) << std::endl;

          }

        }
      }
    }
    for(int k = d.beg[KDIR]; k < d.end[KDIR] ; k++) {
      for(int j = d.beg[JDIR]; j < d.end[JDIR]+JOFFSET ; j++) {
        for(int i = d.beg[IDIR]; i < d.end[IDIR] ; i++) {
          if(myVs(BX2s,k,j,i) != d.Vs(BX2s,k,j,i)) {
            errornum++;
            idfx::cout << "-----------------------------------------" << std::endl
                       << " Error in Vs(BX2s) at (i,j,k) = ( " << i << ", " << j << ", " << k << ")" << std::endl
                       << " Coordinates (x1,x2,x3)   = ( " << d.x[IDIR](i) << ", " << d.x[JDIR](j) << ", " << d.x[KDIR](k) << ")" << std::endl;
          }

        }
      }
    }
    for(int k = d.beg[KDIR]; k < d.end[KDIR]+KOFFSET ; k++) {
      for(int j = d.beg[JDIR]; j < d.end[JDIR] ; j++) {
        for(int i = d.beg[IDIR]; i < d.end[IDIR] ; i++) {
          if(myVs(BX3s,k,j,i) != d.Vs(BX3s,k,j,i)) {
            errornum++;
            idfx::cout << "-----------------------------------------" << std::endl
                       << " Error in Vs(BX3s) at (i,j,k) = ( " << i << ", " << j << ", " << k << ")" << std::endl
                       << " Coordinates (x1,x2,x3)   = ( " << d.x[IDIR](i) << ", " << d.x[JDIR](j) << ", " << d.x[KDIR](k) << ")" << std::endl
                       << "Original= " << myVs(BX3s,k,j,i) << " New=" << d.Vs(BX3s,k,j,i) << " diff=" << myVs(BX3s,k,j,i)-d.Vs(BX3s,k,j,i) << std::endl;
          }

        }
      }
    }
#ifdef EVOLVE_VECTOR_POTENTIAL
    for(int k = d.beg[KDIR]; k < d.end[KDIR]+KOFFSET ; k++) {
      for(int j = d.beg[JDIR]; j < d.end[JDIR]+JOFFSET ; j++) {
        for(int i = d.beg[IDIR]; i < d.end[IDIR] ; i++) {
          if(myVe(AX1e,k,j,i) != d.Ve(AX1e,k,j,i)) {
            errornum++;
            idfx::cout << "-----------------------------------------" << std::endl
                       << " Error in Ve(AX1e) at (i,j,k) = ( " << i << ", " << j << ", " << k << ")" << std::endl
                       << " Coordinates (x1,x2,x3)   = ( " << d.x[IDIR](i) << ", " << d.x[JDIR](j) << ", " << d.x[KDIR](k) << ")" << std::endl;

          }

        }
      }
    }
    for(int k = d.beg[KDIR]; k < d.end[KDIR]+KOFFSET ; k++) {
      for(int j = d.beg[JDIR]; j < d.end[JDIR] ; j++) {
        for(int i = d.beg[IDIR]; i < d.end[IDIR]+IOFFSET ; i++) {
          if(myVe(AX2e,k,j,i) != d.Ve(AX2e,k,j,i)) {
            errornum++;
            idfx::cout << "-----------------------------------------" << std::endl
                       << " Error in Ve(AX2e) at (i,j,k) = ( " << i << ", " << j << ", " << k << ")" << std::endl
                       << " Coordinates (x1,x2,x3)   = ( " << d.x[IDIR](i) << ", " << d.x[JDIR](j) << ", " << d.x[KDIR](k) << ")" << std::endl;
          }

        }
      }
    }
    for(int k = d.beg[KDIR]; k < d.end[KDIR] ; k++) {
      for(int j = d.beg[JDIR]; j < d.end[JDIR]+JOFFSET ; j++) {
        for(int i = d.beg[IDIR]; i < d.end[IDIR]+IOFFSET ; i++) {
          if(myVe(AX3e,k,j,i) != d.Ve(AX3e,k,j,i)) {
            errornum++;
            idfx::cout << "-----------------------------------------" << std::endl
                       << " Error in Ve(AX3e) at (i,j,k) = ( " << i << ", " << j << ", " << k << ")" << std::endl
                       << " Coordinates (x1,x2,x3)   = ( " << d.x[IDIR](i) << ", " << d.x[JDIR](j) << ", " << d.x[KDIR](k) << ")" << std::endl
                       << "Original= " << myVs(BX3s,k,j,i) << " New=" << d.Vs(BX3s,k,j,i) << " diff=" << myVs(BX3s,k,j,i)-d.Vs(BX3s,k,j,i) << std::endl;
          }

        }
      }
    }
#endif

  idfx::cout << "Analysis: consistency check done with " << errornum << " errors " << std::endl;
  if(errornum>0) {
    IDEFIX_ERROR("Restart from dump failed validation");
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
    IdefixHostArray4D<real> Ve;

    #ifndef EVOLVE_VECTOR_POTENTIAL
    Ve = IdefixHostArray4D<real>("Potential vector",3, d.np_tot[KDIR]+1, d.np_tot[JDIR]+1, d.np_tot[IDIR]+1);
    #else
    Ve = d.Ve;
    #endif

    bool haveTracer = data.hydro->haveTracer;

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

                real xl=d.xl[IDIR](i);
                real yl=d.xl[JDIR](j);
                real zl=d.xl[KDIR](k);
                Ve(IDIR,k,j,i) = B0/(2.0*M_PI)*(cos(2.0*M_PI*yl));
                Ve(JDIR,k,j,i) = B0/(2.0*M_PI)*sin(2.0*M_PI*xl);
                Ve(KDIR,k,j,i) = B0/(2.0*M_PI)*(
                                    cos(2.0*M_PI*yl) + cos(4.0*M_PI*xl)/2.0);

                if(haveTracer) {
                  d.Vc(TRG  ,k,j,i) = x>0.5?  1.0:0.0;
                  d.Vc(TRG+1,k,j,i) = z>0.5?  1.0:0.0;
                }
            }
        }
    }

    #ifndef EVOLVE_VECTOR_POTENTIAL
    d.MakeVsFromAmag(Ve);
    #endif
    // Send it all, if needed
    d.SyncToDevice();
}

// Analyse data to produce an output
void MakeAnalysis(DataBlock & data) {

}



// Do a specifically designed user step in the middle of the integration
void ComputeUserStep(DataBlock &data, real t, real dt) {

}
