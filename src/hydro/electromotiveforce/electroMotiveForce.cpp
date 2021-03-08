// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2021 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include "electroMotiveForce.hpp"
#include "hydro.hpp"
#include "dataBlock.hpp"
// Initialisation of electromotive force object
ElectroMotiveForce::ElectroMotiveForce() {
  // Do nothing
}

// Init the emf from a datablock pointer
void ElectroMotiveForce::Init(Hydro *hydro) {
  idfx::pushRegion("ElectroMotiveForce::Init");

#if MHD == YES
  #if EMF_AVERAGE == UCT_CONTACT
    idfx::cout << "ElectroMotiveForce: Using UCT_CONTACT averaging scheme." << std::endl;
  #elif EMF_AVERAGE == UCT0
    idfx::cout << "ElectroMotiveForce: Using UCT0 averaging scheme." << std::endl;
  #elif EMF_AVERAGE == ARITHMETIC
    idfx::cout << "ElectroMotiveForce: Using ARITHMETIC averaging scheme." << std::endl;
  #elif EMF_AVERAGE == UCT_HLL
    idfx::cout << "ElectroMotiveForce: Using 2D-HLL averaging scheme." << std::endl;
  #else
    IDEFIX_ERROR("Unknown EMF averaging scheme in definitions.hpp");
  #endif
  this->data = hydro->data;
  this->hydro = hydro;

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

  #if EMF_AVERAGE == UCT_CONTACT
  D_EXPAND( svx = IdefixArray3D<int>("EMF_svx",
                              data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);  ,
            svy = IdefixArray3D<int>("EMF_svy",
                              data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);  ,
            svz = IdefixArray3D<int>("EMF_svz",
                              data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);  )

  #elif EMF_AVERAGE == UCT_HLL
  D_EXPAND( SxL = IdefixArray3D<real>("EMF_SxL",
                              data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
            SxR = IdefixArray3D<real>("EMF_SxR",
                              data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);  ,

            SyL = IdefixArray3D<real>("EMF_SyL",
                              data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
            SyR = IdefixArray3D<real>("EMF_SyR",
                              data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);  ,

            SzL = IdefixArray3D<real>("EMF_SzL",
                              data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
            SzR = IdefixArray3D<real>("EMF_SzR",
                              data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);  )

  dvx_dx = IdefixArray3D<real>("EMF_dvx_dx",
                               data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
  dvx_dy = IdefixArray3D<real>("EMF_dvx_dy",
                               data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);

  dvy_dx = IdefixArray3D<real>("EMF_dvy_dx",
                               data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
  dvy_dy = IdefixArray3D<real>("EMF_dvy_dy",
                               data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);

  #if DIMENSIONS == 3
  dvx_dz = IdefixArray3D<real>("EMF_dvx_dz",
                               data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
  dvy_dz = IdefixArray3D<real>("EMF_dvy_dz",
                               data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);

  dvz_dx = IdefixArray3D<real>("EMF_dvz_dx",
                               data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
  dvz_dy = IdefixArray3D<real>("EMF_dvz_dy",
                               data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
  dvz_dz = IdefixArray3D<real>("EMF_dvz_dz",
                               data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
  #endif

  dbx_dy = IdefixArray3D<real>("EMF_dbx_dy",
                               data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
  dby_dx = IdefixArray3D<real>("EMF_dby_dx",
                               data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
  #if DIMENSIONS == 3
  dbx_dz = IdefixArray3D<real>("EMF_dbx_dz",
                               data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
  dby_dz = IdefixArray3D<real>("EMF_dby_dz",
                               data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
  dbz_dx = IdefixArray3D<real>("EMF_dbz_dx",
                               data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
  dbz_dy = IdefixArray3D<real>("EMF_dbz_dy",
                               data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
  #endif // DIMENSIONS
#endif // EMF_AVERAGE

  Ex1 = IdefixArray3D<real>("EMF_Ex1", data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
  Ex2 = IdefixArray3D<real>("EMF_Ex2", data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
  Ex3 = IdefixArray3D<real>("EMF_Ex3", data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);

  Grid *mygrid = data->mygrid;

  // MPI initialisation
  #ifdef WITH_MPI

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

  MPI_SAFE_CALL(MPI_Send_init(BufferSendX1[faceRight].data(), bufferSizeX1, realMPI, procSend, 100,
                mygrid->CartComm, &sendRequestX1[faceRight]));

  MPI_SAFE_CALL(MPI_Recv_init(BufferRecvX1[faceLeft].data(), bufferSizeX1, realMPI, procRecv, 100,
                mygrid->CartComm, &recvRequestX1[faceLeft]));

  // Send to the left
  // We receive from procRecv, and we send to procSend
  MPI_SAFE_CALL(MPI_Cart_shift(mygrid->CartComm,0,-1,&procRecv,&procSend ));

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

  #endif // MHD==YES
  idfx::popRegion();
}





// Destructor (clean up persistent communication channels)
ElectroMotiveForce::~ElectroMotiveForce() {
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
