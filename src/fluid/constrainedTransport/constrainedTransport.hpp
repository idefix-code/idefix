// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef FLUID_CONSTRAINEDTRANSPORT_CONSTRAINEDTRANSPORT_HPP_
#define FLUID_CONSTRAINEDTRANSPORT_CONSTRAINEDTRANSPORT_HPP_

#include "idefix.hpp"
#include "input.hpp"
#include "riemannSolver.hpp"

// Forward declarations
#include "physics.hpp"
template <typename Phys> class Fluid;

class DataBlock;

template <typename Phys>
class ConstrainedTransport {
 public:
  enum AveragingType {none, arithmetic, uct0, uct_contact, uct_hll, uct_hlld};

  // Type of averaging
  AveragingType averaging{none};

  // Face centered emf components
  IdefixArray3D<real>     exj;
  IdefixArray3D<real>     exk;
  IdefixArray3D<real>     eyi;
  IdefixArray3D<real>     eyk;
  IdefixArray3D<real>     ezi;
  IdefixArray3D<real>     ezj;

  // Edge centered emf components
  IdefixArray3D<real>     ex;
  IdefixArray3D<real>     ey;
  IdefixArray3D<real>     ez;

// Required by uct_contact averaging
  IdefixArray3D<real>      svx;
  IdefixArray3D<real>      svy;
  IdefixArray3D<real>      svz;

// required by uct_hll averaging
  IdefixArray3D<real> axL;
  IdefixArray3D<real> axR;
  IdefixArray3D<real> ayL;
  IdefixArray3D<real> ayR;
  #if DIMENSIONS == 3
  IdefixArray3D<real> azL;
  IdefixArray3D<real> azR;

  IdefixArray3D<real> dzL;
  IdefixArray3D<real> dzR;
  #endif
  IdefixArray3D<real> dxL;
  IdefixArray3D<real> dxR;
  IdefixArray3D<real> dyL;
  IdefixArray3D<real> dyR;


  IdefixArray3D<real>     Ex1;
  IdefixArray3D<real>     Ex2;
  IdefixArray3D<real>     Ex3;

  // Helper arrays for shearing box
  IdefixArray2D<real>     sbEyL;
  IdefixArray2D<real>     sbEyR;
  IdefixArray2D<real>     sbEyRL;

  // Range of existence

  // Init from Hydro class
  ConstrainedTransport(Input &, Fluid<Phys> *);
  ~ConstrainedTransport();

  void EvolveMagField(real, real, IdefixArray4D<real>&);
  void CalcCornerEMF(real );
  void ShowConfig();

  // Different flavors of EMF average schemes
  void CalcRiemannAverage();
  void CalcArithmeticAverage();
  void CalcCellCenteredEMF();
  void CalcUCT0Average();
  void CalcContactAverage();

  // Enforce boundary conditions on the EMFs.
  void EnforceEMFBoundary();
  void CalcNonidealEMF(real );

  // Specific routines to symmetrize EMFs with shearing box BCs
  void SymmetrizeEMFShearingBox();
  void ExtrapolateEMFShearingBox(BoundarySide,
                            IdefixArray2D<real>,
                            IdefixArray2D<real>);

  // Routines for evolving the magnetic potential (only available when EVOLVE_VECTOR_POTENTIAL)
  void EvolveVectorPotential(real, IdefixArray4D<real> &);
  void ComputeMagFieldFromA(IdefixArray4D<real> &Vein, IdefixArray4D<real> &Vsout);

#ifdef WITH_MPI
  // Exchange surface EMFs to remove interprocess round off errors
  void ExchangeAll();
  void ExchangeX1();
  void ExchangeX2();
  void ExchangeX3();
#endif

 private:
  DataBlock *data;
  Fluid<Phys> *hydro;

#ifdef WITH_MPI
  enum {faceRight, faceLeft};

  // Buffers for MPI calls
  IdefixArray1D<real> BufferSendX1[2];
  IdefixArray1D<real> BufferSendX2[2];
  IdefixArray1D<real> BufferSendX3[2];
  IdefixArray1D<real> BufferRecvX1[2];
  IdefixArray1D<real> BufferRecvX2[2];
  IdefixArray1D<real> BufferRecvX3[2];

  IdefixArray1D<int>  mapVars;

  int bufferSizeX1;
  int bufferSizeX2;
  int bufferSizeX3;

  // Requests for MPI persistent communications
  MPI_Request sendRequestX1[2];
  MPI_Request sendRequestX2[2];
  MPI_Request sendRequestX3[2];
  MPI_Request recvRequestX1[2];
  MPI_Request recvRequestX2[2];
  MPI_Request recvRequestX3[2];

  Kokkos::Timer timer;    // Internal MPI timer
#endif
};


#include <string>
#include "fluid.hpp"
#include "dataBlock.hpp"

// Init the emf from a datablock pointer
template<typename Phys>
ConstrainedTransport<Phys>::ConstrainedTransport(Input &input, Fluid<Phys> *hydro) {
  idfx::pushRegion("ConstrainedTransport::Init");

  if(input.CheckEntry("Hydro","emf")>=0) {
    std::string opType = input.Get<std::string>("Hydro","emf",0);
    if(opType.compare("arithmetic")==0) {
      this->averaging = arithmetic;
    } else if(opType.compare("uct0")==0) {
      this->averaging = uct0;
    } else if(opType.compare("uct_contact")==0) {
      this->averaging = uct_contact;
    } else if(opType.compare("uct_hll")==0) {
      this->averaging = uct_hll;
    } else if(opType.compare("uct_hlld")==0) {
      this->averaging = uct_hlld;
    } else {
      idfx::cout << "ConstrainedTransport: unknown averaging scheme " << opType << std::endl;
      IDEFIX_ERROR("Unknown EMF averaging scheme");
    }
  } else {
    if(hydro->hallStatus.status == HydroModuleStatus::Disabled) {
      // by default, use uct_contact
      this->averaging = uct_contact;
    } else {
      this->averaging = arithmetic;
    }
  }

  this->data = hydro->data;
  this->hydro = hydro;

  // Allocate shearing box arrays
  if(hydro->haveShearingBox == true) {
    sbEyL = IdefixArray2D<real>("EMF_sbEyL", data->np_tot[KDIR], data->np_tot[JDIR]);
    sbEyR = IdefixArray2D<real>("EMF_sbEyR", data->np_tot[KDIR], data->np_tot[JDIR]);
    sbEyRL = IdefixArray2D<real>("EMF_sbEyRL", data->np_tot[KDIR], data->np_tot[JDIR]);
  }

  D_EXPAND( ez = IdefixArray3D<real>("EMF_ez",
                              data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);  ,
                                                                                            ,
            ex = IdefixArray3D<real>("EMF_ex",
                              data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
            ey = IdefixArray3D<real>("EMF_ey",
                              data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);  )

  D_EXPAND( ezi = IdefixArray3D<real>("EMF_ezi",
                              data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
            ezj = IdefixArray3D<real>("EMF_ezj",
                              data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);  ,
                                                                                            ,
            exj = IdefixArray3D<real>("EMF_exj",
                              data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
            exk = IdefixArray3D<real>("EMF_exj",
                              data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
            eyi = IdefixArray3D<real>("EMF_eyi",
                              data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
            eyk = IdefixArray3D<real>("EMF_eyi",
                              data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]); )

  if(averaging==uct_contact) {
    D_EXPAND( svx = IdefixArray3D<real>("EMF_svx",
                                data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);  ,
              svy = IdefixArray3D<real>("EMF_svy",
                                data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);  ,
              svz = IdefixArray3D<real>("EMF_svz",
                                data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);  )
  }


  if(averaging==uct_hll || averaging==uct_hlld) {
    D_EXPAND( axL = IdefixArray3D<real>("EMF_axL",
                                data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
              axR = IdefixArray3D<real>("EMF_axR",
                                data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);  ,

              ayL = IdefixArray3D<real>("EMF_ayL",
                                data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
              ayR = IdefixArray3D<real>("EMF_ayR",
                                data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);  ,

              azL = IdefixArray3D<real>("EMF_azL",
                                data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
              azR = IdefixArray3D<real>("EMF_azR",
                                data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);  )

    D_EXPAND( dxL = IdefixArray3D<real>("EMF_dxL",
                                data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
              dxR = IdefixArray3D<real>("EMF_dxR",
                                data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);  ,

              dyL = IdefixArray3D<real>("EMF_dyL",
                                data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
              dyR = IdefixArray3D<real>("EMF_dyR",
                                data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);  ,

              dzL = IdefixArray3D<real>("EMF_dzL",
                                data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
              dzR = IdefixArray3D<real>("EMF_dzR",
                                data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);  )
  }
  if(averaging==uct_hlld) {
    if(   hydro->rSolver->GetSolver() == RiemannSolver<Phys>::Solver::HLL_MHD
       || hydro->rSolver->GetSolver() == RiemannSolver<Phys>::Solver::TVDLF_MHD) {
      IDEFIX_ERROR("HLLD EMF reconstruction is only compatible with HLLD or ROE Riemann solvers");
    }
  }

  Ex1 = IdefixArray3D<real>("EMF_Ex1", data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
  Ex2 = IdefixArray3D<real>("EMF_Ex2", data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
  Ex3 = IdefixArray3D<real>("EMF_Ex3", data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);

  // MPI initialisation
  #ifdef WITH_MPI

  Grid *mygrid = data->mygrid;
  // init timer
  this->timer.reset();
  /////////////////////////////////////////////////////////////////////////////
  // Init exchange datasets
  bufferSizeX1 = 0;
  bufferSizeX2 = 0;
  bufferSizeX3 = 0;

  // in X1, we exchage ez
  bufferSizeX1 = (data->np_int[JDIR]+JOFFSET) * data->np_int[KDIR];

  #if DIMENSIONS == 3
  // We also have ey
  bufferSizeX1 += data->np_int[JDIR] * (data->np_int[KDIR]+KOFFSET);
  #endif

  BufferRecvX1[faceLeft ] = IdefixArray1D<real>("EmfRecvX1Left", bufferSizeX1);
  BufferRecvX1[faceRight] = IdefixArray1D<real>("EmfRecvX1Right",bufferSizeX1);
  BufferSendX1[faceLeft ] = IdefixArray1D<real>("EmfSendX1Left", bufferSizeX1);
  BufferSendX1[faceRight] = IdefixArray1D<real>("EmfSendX1Right",bufferSizeX1);

  // Number of cells in X2 boundary condition (only required when problem >2D):
#if DIMENSIONS >= 2

  bufferSizeX2 = (data->np_int[IDIR]+IOFFSET) * data->np_int[KDIR];

  #if DIMENSIONS == 3
  // We also have ex
  bufferSizeX2 += data->np_int[IDIR] * (data->np_int[KDIR]+KOFFSET);
  #endif

  BufferRecvX2[faceLeft ] = IdefixArray1D<real>("EmfRecvX2Left", bufferSizeX2);
  BufferRecvX2[faceRight] = IdefixArray1D<real>("EmfRecvX2Right",bufferSizeX2);
  BufferSendX2[faceLeft ] = IdefixArray1D<real>("EmfSendX2Left", bufferSizeX2);
  BufferSendX2[faceRight] = IdefixArray1D<real>("EmfSendX2Right",bufferSizeX2);

#endif
// Number of cells in X3 boundary condition (only required when problem is 3D):
#if DIMENSIONS ==3
  // ex
  bufferSizeX3 = data->np_int[IDIR] * (data->np_int[JDIR]+JOFFSET);
  // ey
  bufferSizeX3 += (data->np_int[IDIR]+IOFFSET) * data->np_int[JDIR];

  BufferRecvX3[faceLeft ] = IdefixArray1D<real>("EmfRecvX3Left", bufferSizeX3);
  BufferRecvX3[faceRight] = IdefixArray1D<real>("EmfRecvX3Right",bufferSizeX3);
  BufferSendX3[faceLeft ] = IdefixArray1D<real>("EmfSendX3Left", bufferSizeX3);
  BufferSendX3[faceRight] = IdefixArray1D<real>("EmfSendX3Right",bufferSizeX3);
#endif // DIMENSIONS

  // Init persistent MPI communications
  int procSend, procRecv;

  // X1-dir exchanges
  // We receive from procRecv, and we send to procSend
  MPI_SAFE_CALL(MPI_Cart_shift(mygrid->CartComm,0,1,&procRecv,&procSend ));
  // Avoid the exchange when shearing box is enabled (so that EMFs don't get averaged)
  if(data->rbound[IDIR] == shearingbox ) procSend = MPI_PROC_NULL;
  if(data->lbound[IDIR] == shearingbox ) procRecv = MPI_PROC_NULL;

  MPI_SAFE_CALL(MPI_Send_init(BufferSendX1[faceRight].data(), bufferSizeX1, realMPI, procSend, 100,
                mygrid->CartComm, &sendRequestX1[faceRight]));

  MPI_SAFE_CALL(MPI_Recv_init(BufferRecvX1[faceLeft].data(), bufferSizeX1, realMPI, procRecv, 100,
                mygrid->CartComm, &recvRequestX1[faceLeft]));

  // Send to the left
  // We receive from procRecv, and we send to procSend
  MPI_SAFE_CALL(MPI_Cart_shift(mygrid->CartComm,0,-1,&procRecv,&procSend ));
  // Avoid the exchange when shearing box is enabled (so that EMFs don't get averaged)
  if(data->lbound[IDIR] == shearingbox ) procSend = MPI_PROC_NULL;
  if(data->rbound[IDIR] == shearingbox ) procRecv = MPI_PROC_NULL;

  MPI_SAFE_CALL(MPI_Send_init(BufferSendX1[faceLeft].data(), bufferSizeX1, realMPI, procSend, 101,
                mygrid->CartComm, &sendRequestX1[faceLeft]));

  MPI_SAFE_CALL(MPI_Recv_init(BufferRecvX1[faceRight].data(), bufferSizeX1, realMPI, procRecv, 101,
                mygrid->CartComm, &recvRequestX1[faceRight]));

  #if DIMENSIONS >= 2
  // We receive from procRecv, and we send to procSend
  MPI_SAFE_CALL(MPI_Cart_shift(mygrid->CartComm,1,1,&procRecv,&procSend ));

  MPI_SAFE_CALL(MPI_Send_init(BufferSendX2[faceRight].data(), bufferSizeX2, realMPI, procSend, 200,
                mygrid->CartComm, &sendRequestX2[faceRight]));

  MPI_SAFE_CALL(MPI_Recv_init(BufferRecvX2[faceLeft].data(), bufferSizeX2, realMPI, procRecv, 200,
                mygrid->CartComm, &recvRequestX2[faceLeft]));

  // Send to the left
  // We receive from procRecv, and we send to procSend
  MPI_SAFE_CALL(MPI_Cart_shift(mygrid->CartComm,1,-1,&procRecv,&procSend ));

  MPI_SAFE_CALL(MPI_Send_init(BufferSendX2[faceLeft].data(), bufferSizeX2, realMPI, procSend, 201,
                mygrid->CartComm, &sendRequestX2[faceLeft]));

  MPI_SAFE_CALL(MPI_Recv_init(BufferRecvX2[faceRight].data(), bufferSizeX2, realMPI, procRecv, 201,
                mygrid->CartComm, &recvRequestX2[faceRight]));
  #endif

  #if DIMENSIONS == 3
  // We receive from procRecv, and we send to procSend
  MPI_SAFE_CALL(MPI_Cart_shift(mygrid->CartComm,2,1,&procRecv,&procSend ));

  MPI_SAFE_CALL(MPI_Send_init(BufferSendX3[faceRight].data(), bufferSizeX3, realMPI, procSend, 300,
                mygrid->CartComm, &sendRequestX3[faceRight]));

  MPI_SAFE_CALL(MPI_Recv_init(BufferRecvX3[faceLeft].data(), bufferSizeX3, realMPI, procRecv, 300,
                mygrid->CartComm, &recvRequestX3[faceLeft]));

  // Send to the left
  // We receive from procRecv, and we send to procSend
  MPI_SAFE_CALL(MPI_Cart_shift(mygrid->CartComm,2,-1,&procRecv,&procSend ));

  MPI_SAFE_CALL(MPI_Send_init(BufferSendX3[faceLeft].data(), bufferSizeX3, realMPI, procSend, 301,
                mygrid->CartComm, &sendRequestX3[faceLeft]));

  MPI_SAFE_CALL(MPI_Recv_init(BufferRecvX3[faceRight].data(), bufferSizeX3, realMPI, procRecv, 301,
                mygrid->CartComm, &recvRequestX3[faceRight]));
  #endif

#endif  // WITH_MPI

  idfx::popRegion();
}





// Destructor (clean up persistent communication channels)
template<typename Phys>
ConstrainedTransport<Phys>::~ConstrainedTransport() {
  #if MHD == YES
    #ifdef WITH_MPI
    // Properly clean up the mess
    idfx::cout << "Emf: Cleaning up MPI persistent communication channels" << std::endl;
    for(int i=0 ; i< 2; i++) {
      MPI_Request_free( &sendRequestX1[i]);
      MPI_Request_free( &recvRequestX1[i]);

    #if DIMENSIONS >= 2
      MPI_Request_free( &sendRequestX2[i]);
      MPI_Request_free( &recvRequestX2[i]);
    #endif

    #if DIMENSIONS == 3
      MPI_Request_free( &sendRequestX3[i]);
      MPI_Request_free( &recvRequestX3[i]);
    #endif
    }
    #endif
  #endif
}

template<typename Phys>
void ConstrainedTransport<Phys>::ShowConfig() {
  switch(averaging) {
    case arithmetic:
      idfx::cout << "ConstrainedTransport: Using ARITHMETIC averaging scheme." << std::endl;
      break;
    case uct0:
      idfx::cout << "ConstrainedTransport: Using UCT0 averaging scheme." << std::endl;
      break;
    case uct_contact:
      idfx::cout << "ConstrainedTransport: Using UCT_CONTACT averaging scheme." << std::endl;
      break;
    case uct_hll:
      idfx::cout << "ConstrainedTransport: Using 2D-HLL averaging scheme." << std::endl;
      break;
    case uct_hlld:
      idfx::cout << "ConstrainedTransport: Using 2D-HLLD averaging scheme." << std::endl;
      break;
    default:
      IDEFIX_ERROR("Unknown averaging scheme");
  }
}

#include "calcCornerEmf.hpp"
#include "calcNonidealEMF.hpp"
#include "calcRiemannEmf.hpp"
#include "EMFexchange.hpp"
#include "enforceEMFBoundary.hpp"
#include "evolveMagField.hpp"
#include "evolveVectorPotential.hpp"

#endif // FLUID_CONSTRAINEDTRANSPORT_CONSTRAINEDTRANSPORT_HPP_
