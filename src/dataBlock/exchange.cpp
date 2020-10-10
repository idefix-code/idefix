// ********************************************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ********************************************************************************************************

#include "../idefix.hpp"
#include "dataBlock.hpp"

//#define MPI_NON_BLOCKING
// MPI Routines exchange
void DataBlock::ExchangeAll(){
#ifdef WITH_MPI

#else
    IDEFIX_ERROR("You should not call DataBlock exchange routines without MPI");
#endif
}

void DataBlock::ExchangeX1() {
    idfx::pushRegion("DataBlock::ExchangeX1");
#ifdef WITH_MPI
    // init MPI Timer
    idfx::mpiTimer -= timer.seconds();
    // Load  the buffers with data
    int ibeg,iend,jbeg,jend,kbeg,kend,offset;
    int nx,ny,nz;
    IdefixArray1D<real> BufferLeft=BufferSendX1[faceLeft];
    IdefixArray1D<real> BufferRight=BufferSendX1[faceRight];
    IdefixArray4D<real> Vc=this->Vc;
#if MHD==YES
    IdefixArray4D<real> Vs=this->Vs;
    int VsIndex;
#endif
    // Coordinates of the ghost region which needs to be transfered
    ibeg = 0;
    iend = nghost[IDIR]; 
    nx = nghost[IDIR];                    // Number of points in x
    offset = end[IDIR];   // Distance between beginning of left and right ghosts
    jbeg = beg[JDIR];
    jend = end[JDIR];
    ny = jend - jbeg;
    kbeg = beg[KDIR];
    kend = end[KDIR];   
    nz = kend - kbeg;

    idefix_for("LoadBufferX1Vc",0,NVAR,kbeg,kend,jbeg,jend,ibeg,iend,
                    KOKKOS_LAMBDA (int n, int k, int j, int i) {
                        BufferLeft(i + (j-jbeg)*nx + (k-kbeg)*nx*ny + n*nx*ny*nz ) = Vc(n,k,j,i+nx);
                        BufferRight(i + (j-jbeg)*nx + (k-kbeg)*nx*ny + n*nx*ny*nz ) = Vc(n,k,j,i+offset-nx);
                    });

#if MHD==YES
    // Load face-centered field in the buffer
    
    VsIndex = NVAR*nx*ny*nz;

    idefix_for("LoadBufferX1BX1s",kbeg,kend,jbeg,jend,ibeg,iend,
                    KOKKOS_LAMBDA (int k, int j, int i) {
                        BufferLeft(i + (j-jbeg)*nx + (k-kbeg)*nx*ny + VsIndex ) = Vs(BX1s,k,j,i+nx+1);
                        BufferRight(i + (j-jbeg)*nx + (k-kbeg)*nx*ny + VsIndex ) = Vs(BX1s,k,j,i+offset-nx);
                    });
    
    #if DIMENSIONS >= 2
    VsIndex = (NVAR+1)*nx*ny*nz;

    idefix_for("LoadBufferX1BX2s",kbeg,kend,jbeg,jend+1,ibeg,iend,
                    KOKKOS_LAMBDA (int k, int j, int i) {
                        BufferLeft(i + (j-jbeg)*nx + (k-kbeg)*nx*(ny+1) + VsIndex ) = Vs(BX2s,k,j,i+nx);
                        BufferRight(i + (j-jbeg)*nx + (k-kbeg)*nx*(ny+1) + VsIndex ) = Vs(BX2s,k,j,i+offset-nx);
                    });

    #endif

    #if DIMENSIONS == 3
    VsIndex = (NVAR+1)*nx*ny*nz + nx*(ny+1)*nz;

    idefix_for("LoadBufferX1BX3s",kbeg,kend+1,jbeg,jend,ibeg,iend,
                    KOKKOS_LAMBDA (int k, int j, int i) {
                        BufferLeft(i + (j-jbeg)*nx + (k-kbeg)*nx*ny + VsIndex ) = Vs(BX3s,k,j,i+nx);
                        BufferRight(i + (j-jbeg)*nx + (k-kbeg)*nx*ny + VsIndex ) = Vs(BX3s,k,j,i+offset-nx);
                    });

    #endif
#endif

   
    int procSend, procRecv;
    
    
    Kokkos::fence();
    #ifdef MPI_NON_BLOCKING
        MPI_Status sendStatus[2];
        MPI_Status recvStatus[2];
        MPI_Request sendRequest[2];
        MPI_Request recvRequest[2];


        MPI_SAFE_CALL(MPI_Cart_shift(mygrid->CartComm,0,1,&procRecv,&procSend ));   // We receive from procRecv, and we send to procSend
        MPI_SAFE_CALL(MPI_Isend(BufferSendX1[faceRight].data(), bufferSizeX1, realMPI, procSend, 100,
                    mygrid->CartComm, &sendRequest[0]));
        
        MPI_SAFE_CALL(MPI_Irecv(BufferRecvX1[faceLeft].data(), bufferSizeX1, realMPI, procRecv, 100,
                    mygrid->CartComm, &recvRequest[0]));

        // Send to the left
        MPI_SAFE_CALL(MPI_Cart_shift(mygrid->CartComm,0,-1,&procRecv,&procSend ));   // We receive from procRecv, and we send to procSend
        MPI_SAFE_CALL(MPI_Isend(BufferSendX1[faceLeft].data(), bufferSizeX1, realMPI, procSend, 101,
                    mygrid->CartComm, &sendRequest[1]));
        
        MPI_SAFE_CALL(MPI_Irecv(BufferRecvX1[faceRight].data(), bufferSizeX1, realMPI, procRecv, 101,
                    mygrid->CartComm, &recvRequest[1]));

                
        // Wait for recv to complete (we don't care about the sends)
        MPI_Waitall(2, recvRequest, recvStatus);
    #else
        MPI_Status status;
         // Send to the right
        MPI_SAFE_CALL(MPI_Cart_shift(mygrid->CartComm,0,1,&procRecv,&procSend ));   // We receive from procRecv, and we send to procSend
        MPI_SAFE_CALL(MPI_Sendrecv(BufferSendX1[faceRight].data(), bufferSizeX1, realMPI, procSend, 100,
                    BufferRecvX1[faceLeft].data(), bufferSizeX1, realMPI, procRecv, 100,
                    mygrid->CartComm, &status));
        
        // Send to the left
        MPI_SAFE_CALL(MPI_Cart_shift(mygrid->CartComm,0,-1,&procRecv,&procSend ));   // We receive from procRecv, and we send to procSend
        MPI_SAFE_CALL(MPI_Sendrecv(BufferSendX1[faceLeft].data(), bufferSizeX1, realMPI, procSend, 101,
                    BufferRecvX1[faceRight].data(), bufferSizeX1, realMPI, procRecv, 101,
                    mygrid->CartComm, &status));
    #endif
    
    // Unpack
    BufferLeft=BufferRecvX1[faceLeft];
    BufferRight=BufferRecvX1[faceRight];

    // We fill the ghost zones

    idefix_for("StoreBufferX1Vc",0,NVAR,kbeg,kend,jbeg,jend,ibeg,iend,
                    KOKKOS_LAMBDA (int n, int k, int j, int i) {
                        Vc(n,k,j,i) = BufferLeft(i + (j-jbeg)*nx + (k-kbeg)*nx*ny + n*nx*ny*nz );
                        Vc(n,k,j,i+offset) = BufferRight(i + (j-jbeg)*nx + (k-kbeg)*nx*ny + n*nx*ny*nz );
                    });

#if MHD==YES
    // Load face-centered field in the buffer
    
    VsIndex = NVAR*nx*ny*nz;

    idefix_for("StoreBufferX1BX1s",kbeg,kend,jbeg,jend,ibeg,iend,
                    KOKKOS_LAMBDA (int k, int j, int i) {
                        Vs(BX1s,k,j,i) = BufferLeft(i + (j-jbeg)*nx + (k-kbeg)*nx*ny + VsIndex );
                        Vs(BX1s,k,j,i+offset+1) = BufferRight(i + (j-jbeg)*nx + (k-kbeg)*nx*ny + VsIndex );
                    });
    
    #if DIMENSIONS >= 2
    VsIndex = (NVAR+1)*nx*ny*nz;

    idefix_for("StoreBufferX1BX2s",kbeg,kend,jbeg,jend+1,ibeg,iend,
                    KOKKOS_LAMBDA (int k, int j, int i) {
                        Vs(BX2s,k,j,i) = BufferLeft(i + (j-jbeg)*nx + (k-kbeg)*nx*(ny+1) + VsIndex );
                        Vs(BX2s,k,j,i+offset) = BufferRight(i + (j-jbeg)*nx + (k-kbeg)*nx*(ny+1) + VsIndex );
                    });

    #endif

    #if DIMENSIONS == 3
    VsIndex = (NVAR+1)*nx*ny*nz + nx*(ny+1)*nz;

    idefix_for("StoreBufferX1BX3s",kbeg,kend+1,jbeg,jend,ibeg,iend,
                    KOKKOS_LAMBDA ( int k, int j, int i) {
                        Vs(BX3s,k,j,i) = BufferLeft(i + (j-jbeg)*nx + (k-kbeg)*nx*ny + VsIndex );
                        Vs(BX3s,k,j,i+offset) = BufferRight(i + (j-jbeg)*nx + (k-kbeg)*nx*ny + VsIndex );
                    });

    #endif
#endif
    #ifdef MPI_NON_BLOCKING
        // Wait for the sends if they have not yet completed
        MPI_Waitall(2, sendRequest, sendStatus);
    #endif

    // Stop MPI Timer
    idfx::mpiTimer += timer.seconds();

#else
    IDEFIX_ERROR("You should not call DataBlock exchange routines without MPI");
#endif
    idfx::popRegion();
}


void DataBlock::ExchangeX2() {
    idfx::pushRegion("DataBlock::ExchangeX2");
#ifdef WITH_MPI
    // init MPI Timer
    idfx::mpiTimer -= timer.seconds();
    // Load  the buffers with data
    int ibeg,iend,jbeg,jend,kbeg,kend,offset;
    int nx,ny,nz;
    IdefixArray1D<real> BufferLeft=BufferSendX2[faceLeft];
    IdefixArray1D<real> BufferRight=BufferSendX2[faceRight];
    IdefixArray4D<real> Vc=this->Vc;
#if MHD==YES
    IdefixArray4D<real> Vs=this->Vs;
    int VsIndex;
#endif
    // Coordinates of the ghost region which needs to be transfered
    ibeg = 0;
    iend = np_tot[IDIR]; 
    nx = np_tot[IDIR];                    // Number of points in x
    jbeg = 0;
    jend = nghost[JDIR];
    offset = end[JDIR];   // Distance between beginning of left and right ghosts
    ny = nghost[JDIR];
    kbeg = beg[KDIR];
    kend = end[KDIR];   
    nz = kend - kbeg;

    idefix_for("LoadBufferX2Vc",0,NVAR,kbeg,kend,jbeg,jend,ibeg,iend,
                    KOKKOS_LAMBDA (int n, int k, int j, int i) {
                        BufferLeft(i + (j-jbeg)*nx + (k-kbeg)*nx*ny + n*nx*ny*nz ) = Vc(n,k,j+ny,i);
                        BufferRight(i + (j-jbeg)*nx + (k-kbeg)*nx*ny + n*nx*ny*nz ) = Vc(n,k,j+offset-ny,i);
                    });

#if MHD==YES
    // Load face-centered field in the buffer
    
    VsIndex = NVAR*nx*ny*nz;

    idefix_for("LoadBufferX2BX1s",kbeg,kend,jbeg,jend,ibeg,iend+1,
                    KOKKOS_LAMBDA (int k, int j, int i) {
                        BufferLeft(i + (j-jbeg)*(nx+1) + (k-kbeg)*(nx+1)*ny + VsIndex ) = Vs(BX1s,k,j+ny,i);
                        BufferRight(i + (j-jbeg)*(nx+1) + (k-kbeg)*(nx+1)*ny + VsIndex ) = Vs(BX1s,k,j+offset-ny,i);
                    });
    
    #if DIMENSIONS >= 2
    VsIndex = NVAR*nx*ny*nz + (nx+1)*ny*nz;

    idefix_for("LoadBufferX2BX2s",kbeg,kend,jbeg,jend,ibeg,iend,
                    KOKKOS_LAMBDA (int k, int j, int i) {
                        BufferLeft(i + (j-jbeg)*nx + (k-kbeg)*nx*ny + VsIndex ) = Vs(BX2s,k,j+ny+1,i);
                        BufferRight(i + (j-jbeg)*nx + (k-kbeg)*nx*ny + VsIndex ) = Vs(BX2s,k,j+offset-ny,i);
                    });

    #endif

    #if DIMENSIONS == 3
    VsIndex = (NVAR+1)*nx*ny*nz + (nx+1)*ny*nz;

    idefix_for("LoadBufferX2BX3s",kbeg,kend+1,jbeg,jend,ibeg,iend,
                    KOKKOS_LAMBDA (int k, int j, int i) {
                        BufferLeft(i + (j-jbeg)*nx + (k-kbeg)*nx*ny + VsIndex ) = Vs(BX3s,k,j+ny,i);
                        BufferRight(i + (j-jbeg)*nx + (k-kbeg)*nx*ny + VsIndex ) = Vs(BX3s,k,j+offset-ny,i);
                    });

    #endif
#endif

    
    // Send to the right
    int procSend, procRecv;

    Kokkos::fence();
    #ifdef MPI_NON_BLOCKING
        MPI_Status sendStatus[2];
        MPI_Status recvStatus[2];
        MPI_Request sendRequest[2];
        MPI_Request recvRequest[2];


        MPI_SAFE_CALL(MPI_Cart_shift(mygrid->CartComm,1,1,&procRecv,&procSend ));   // We receive from procRecv, and we send to procSend
        MPI_SAFE_CALL(MPI_Isend(BufferSendX2[faceRight].data(), bufferSizeX2, realMPI, procSend, 100,
                    mygrid->CartComm, &sendRequest[0]));
        
        MPI_SAFE_CALL(MPI_Irecv(BufferRecvX2[faceLeft].data(), bufferSizeX2, realMPI, procRecv, 100,
                    mygrid->CartComm, &recvRequest[0]));

        // Send to the left
        MPI_SAFE_CALL(MPI_Cart_shift(mygrid->CartComm,1,-1,&procRecv,&procSend ));   // We receive from procRecv, and we send to procSend
        MPI_SAFE_CALL(MPI_Isend(BufferSendX2[faceLeft].data(), bufferSizeX2, realMPI, procSend, 101,
                    mygrid->CartComm, &sendRequest[1]));
        
        MPI_SAFE_CALL(MPI_Irecv(BufferRecvX2[faceRight].data(), bufferSizeX2, realMPI, procRecv, 101,
                    mygrid->CartComm, &recvRequest[1]));

                
        // Wait for recv to complete (we don't care about the sends)
        MPI_Waitall(2, recvRequest, recvStatus);
    #else
        MPI_Status status;
        MPI_SAFE_CALL(MPI_Cart_shift(mygrid->CartComm,1,1,&procRecv,&procSend ));   // We receive from procRecv, and we send to procSend
        
        MPI_SAFE_CALL(MPI_Sendrecv(BufferSendX2[faceRight].data(), bufferSizeX2, realMPI, procSend, 200,
                    BufferRecvX2[faceLeft].data(), bufferSizeX2, realMPI, procRecv, 200,
                    mygrid->CartComm, &status));
        
        // Send to the left
        MPI_SAFE_CALL(MPI_Cart_shift(mygrid->CartComm,1,-1,&procRecv,&procSend ));   // We receive from procRecv, and we send to procSend
        MPI_SAFE_CALL(MPI_Sendrecv(BufferSendX2[faceLeft].data(), bufferSizeX2, realMPI, procSend, 201,
                    BufferRecvX2[faceRight].data(), bufferSizeX2, realMPI, procRecv, 201,
                    mygrid->CartComm, &status));
    #endif

    
    // Unpack
    BufferLeft=BufferRecvX2[faceLeft];
    BufferRight=BufferRecvX2[faceRight];

    // We fill the ghost zones

    idefix_for("StoreBufferX2Vc",0,NVAR,kbeg,kend,jbeg,jend,ibeg,iend,
                    KOKKOS_LAMBDA (int n, int k, int j, int i) {
                        Vc(n,k,j,i) = BufferLeft(i + (j-jbeg)*nx + (k-kbeg)*nx*ny + n*nx*ny*nz );
                        Vc(n,k,j+offset,i) = BufferRight(i + (j-jbeg)*nx + (k-kbeg)*nx*ny + n*nx*ny*nz );
                    });

#if MHD==YES
    // Load face-centered field in the buffer
    
    VsIndex = NVAR*nx*ny*nz;

    idefix_for("StoreBufferX2BX1s",kbeg,kend,jbeg,jend,ibeg,iend+1,
                    KOKKOS_LAMBDA (int k, int j, int i) {
                        Vs(BX1s,k,j,i) = BufferLeft(i + (j-jbeg)*(nx+1) + (k-kbeg)*(nx+1)*ny + VsIndex );
                        Vs(BX1s,k,j+offset,i) = BufferRight(i + (j-jbeg)*(nx+1) + (k-kbeg)*(nx+1)*ny + VsIndex );
                    });
    
    #if DIMENSIONS >= 2
    VsIndex = NVAR*nx*ny*nz + (nx+1)*ny*nz;

    idefix_for("StoreBufferX2BX2s",kbeg,kend,jbeg,jend,ibeg,iend,
                    KOKKOS_LAMBDA (int k, int j, int i) {
                        Vs(BX2s,k,j,i) = BufferLeft(i + (j-jbeg)*nx + (k-kbeg)*nx*ny + VsIndex );
                        Vs(BX2s,k,j+offset+1,i) = BufferRight(i + (j-jbeg)*nx + (k-kbeg)*nx*ny + VsIndex );
                    });

    #endif

    #if DIMENSIONS == 3
    VsIndex = (NVAR+1)*nx*ny*nz + (nx+1)*ny*nz;

    idefix_for("StoreBufferX2BX3s",kbeg,kend+1,jbeg,jend,ibeg,iend,
                    KOKKOS_LAMBDA ( int k, int j, int i) {
                        Vs(BX3s,k,j,i) = BufferLeft(i + (j-jbeg)*nx + (k-kbeg)*nx*ny + VsIndex );
                        Vs(BX3s,k,j+offset,i) = BufferRight(i + (j-jbeg)*nx + (k-kbeg)*nx*ny + VsIndex );
                    });

    #endif

#endif
    #ifdef MPI_NON_BLOCKING
        // Wait for the sends if they have not yet completed
        MPI_Waitall(2, sendRequest, sendStatus);
    #endif
    // Stop MPI Timer
    idfx::mpiTimer += timer.seconds();

#else
    IDEFIX_ERROR("You should not call DataBlock exchange routines without MPI");
#endif
    idfx::popRegion();
}


void DataBlock::ExchangeX3(){
    idfx::pushRegion("DataBlock::ExchangeX3");
#ifdef WITH_MPI
    // init MPI Timer
    idfx::mpiTimer -= timer.seconds();
    // Load  the buffers with data
    int ibeg,iend,jbeg,jend,kbeg,kend,offset;
    int nx,ny,nz;
    IdefixArray1D<real> BufferLeft=BufferSendX3[faceLeft];
    IdefixArray1D<real> BufferRight=BufferSendX3[faceRight];
    IdefixArray4D<real> Vc=this->Vc;
#if MHD==YES
    IdefixArray4D<real> Vs=this->Vs;
    int VsIndex;
#endif
    // Coordinates of the ghost region which needs to be transfered
    ibeg = 0;
    iend = np_tot[IDIR]; 
    nx = np_tot[IDIR];                    // Number of points in x
    jbeg = 0;
    jend = np_tot[JDIR];  
    ny = np_tot[JDIR];
    kbeg = 0;
    kend = nghost[KDIR];
    offset = end[KDIR];   // Distance between beginning of left and right ghosts
    nz = nghost[KDIR];

    idefix_for("LoadBufferX3Vc",0,NVAR,kbeg,kend,jbeg,jend,ibeg,iend,
                    KOKKOS_LAMBDA (int n, int k, int j, int i) {
                        BufferLeft((i-ibeg) + (j-jbeg)*nx + (k-kbeg)*nx*ny + n*nx*ny*nz ) = Vc(n,k+nz,j,i);
                        BufferRight((i-ibeg) + (j-jbeg)*nx + (k-kbeg)*nx*ny + n*nx*ny*nz ) = Vc(n,k+offset-nz,j,i);
                    });

#if MHD==YES
    // Load face-centered field in the buffer
    
    VsIndex = NVAR*nx*ny*nz;

    idefix_for("LoadBufferX3BX1s",kbeg,kend,jbeg,jend,ibeg,iend+1,
                    KOKKOS_LAMBDA (int k, int j, int i) {
                        BufferLeft(i + (j-jbeg)*(nx+1) + (k-kbeg)*(nx+1)*ny + VsIndex ) = Vs(BX1s,k+nz,j,i);
                        BufferRight(i + (j-jbeg)*(nx+1) + (k-kbeg)*(nx+1)*ny + VsIndex ) = Vs(BX1s,k+offset-nz,j,i);
                    });
    
    #if DIMENSIONS >= 2
    VsIndex = NVAR*nx*ny*nz + (nx+1)*ny*nz;

    idefix_for("LoadBufferX3BX2s",kbeg,kend,jbeg,jend+1,ibeg,iend,
                    KOKKOS_LAMBDA (int k, int j, int i) {
                        BufferLeft(i + (j-jbeg)*nx + (k-kbeg)*nx*(ny+1) + VsIndex ) = Vs(BX2s,k+nz,j,i);
                        BufferRight(i + (j-jbeg)*nx + (k-kbeg)*nx*(ny+1) + VsIndex ) = Vs(BX2s,k+offset-nz,j,i);
                    });

    #endif

    #if DIMENSIONS == 3
    VsIndex = NVAR*nx*ny*nz + (nx+1)*ny*nz + nx*(ny+1)*nz;

    idefix_for("LoadBufferX3BX3s",kbeg,kend,jbeg,jend,ibeg,iend,
                    KOKKOS_LAMBDA (int k, int j, int i) {
                        BufferLeft(i + (j-jbeg)*nx + (k-kbeg)*nx*ny + VsIndex ) = Vs(BX3s,k+nz+1,j,i);
                        BufferRight(i + (j-jbeg)*nx + (k-kbeg)*nx*ny + VsIndex ) = Vs(BX3s,k+offset-nz,j,i);
                    });

    #endif
#endif

    // Send to the right
    int procSend, procRecv;
    
    Kokkos::fence();
    #ifdef MPI_NON_BLOCKING
        MPI_Status sendStatus[2];
        MPI_Status recvStatus[2];
        MPI_Request sendRequest[2];
        MPI_Request recvRequest[2];


        MPI_SAFE_CALL(MPI_Cart_shift(mygrid->CartComm,2,1,&procRecv,&procSend ));   // We receive from procRecv, and we send to procSend
        MPI_SAFE_CALL(MPI_Isend(BufferSendX3[faceRight].data(), bufferSizeX3, realMPI, procSend, 100,
                    mygrid->CartComm, &sendRequest[0]));
        
        MPI_SAFE_CALL(MPI_Irecv(BufferRecvX3[faceLeft].data(), bufferSizeX3, realMPI, procRecv, 100,
                    mygrid->CartComm, &recvRequest[0]));

        // Send to the left
        MPI_SAFE_CALL(MPI_Cart_shift(mygrid->CartComm,2,-1,&procRecv,&procSend ));   // We receive from procRecv, and we send to procSend
        MPI_SAFE_CALL(MPI_Isend(BufferSendX3[faceLeft].data(), bufferSizeX3, realMPI, procSend, 101,
                    mygrid->CartComm, &sendRequest[1]));
        
        MPI_SAFE_CALL(MPI_Irecv(BufferRecvX3[faceRight].data(), bufferSizeX3, realMPI, procRecv, 101,
                    mygrid->CartComm, &recvRequest[1]));

                
        // Wait for recv to complete (we don't care about the sends)
        MPI_Waitall(2, recvRequest, recvStatus);
    #else
        MPI_Status status;
        MPI_SAFE_CALL(MPI_Cart_shift(mygrid->CartComm,2,1,&procRecv,&procSend ));   // We receive from procRecv, and we send to procSend

        MPI_SAFE_CALL(MPI_Sendrecv(BufferSendX3[faceRight].data(), bufferSizeX3, realMPI, procSend, 300,
                    BufferRecvX3[faceLeft].data(), bufferSizeX3, realMPI, procRecv, 300,
                    mygrid->CartComm, &status));
        
        // Send to the left
        MPI_SAFE_CALL(MPI_Cart_shift(mygrid->CartComm,2,-1,&procRecv,&procSend ));   // We receive from procRecv, and we send to procSend
        MPI_SAFE_CALL(MPI_Sendrecv(BufferSendX3[faceLeft].data(), bufferSizeX3, realMPI, procSend, 301,
                    BufferRecvX3[faceRight].data(), bufferSizeX3, realMPI, procRecv, 301,
                    mygrid->CartComm, &status));
    #endif

    
    // Unpack
    BufferLeft=BufferRecvX3[faceLeft];
    BufferRight=BufferRecvX3[faceRight];

    // We fill the ghost zones

    idefix_for("StoreBufferX3Vc",0,NVAR,kbeg,kend,jbeg,jend,ibeg,iend,
                    KOKKOS_LAMBDA (int n, int k, int j, int i) {
                        Vc(n,k,j,i) = BufferLeft((i-ibeg) + (j-jbeg)*nx + (k-kbeg)*nx*ny + n*nx*ny*nz );
                        Vc(n,k+offset,j,i) = BufferRight((i-ibeg) + (j-jbeg)*nx + (k-kbeg)*nx*ny + n*nx*ny*nz );
                    });

#if MHD==YES
    // Load face-centered field in the buffer
    
    VsIndex = NVAR*nx*ny*nz;

    idefix_for("StoreBufferX3BX1s",kbeg,kend,jbeg,jend,ibeg,iend+1,
                    KOKKOS_LAMBDA (int k, int j, int i) {
                        Vs(BX1s,k,j,i) = BufferLeft(i + (j-jbeg)*(nx+1) + (k-kbeg)*(nx+1)*ny + VsIndex );
                        Vs(BX1s,k+offset,j,i) = BufferRight(i + (j-jbeg)*(nx+1) + (k-kbeg)*(nx+1)*ny + VsIndex );
                    });
    
    #if DIMENSIONS >= 2
    VsIndex = NVAR*nx*ny*nz + (nx+1)*ny*nz;

    idefix_for("StoreBufferX3BX2s",kbeg,kend,jbeg,jend+1,ibeg,iend,
                    KOKKOS_LAMBDA (int k, int j, int i) {
                        Vs(BX2s,k,j,i) = BufferLeft(i + (j-jbeg)*nx + (k-kbeg)*nx*(ny+1) + VsIndex );
                        Vs(BX2s,k+offset,j,i) = BufferRight(i + (j-jbeg)*nx + (k-kbeg)*nx*(ny+1) + VsIndex );
                    });

    #endif

    #if DIMENSIONS == 3
    VsIndex = NVAR*nx*ny*nz + (nx+1)*ny*nz + nx*(ny+1)*nz;

    idefix_for("StoreBufferX3BX3s",kbeg,kend,jbeg,jend,ibeg,iend,
                    KOKKOS_LAMBDA ( int k, int j, int i) {
                        Vs(BX3s,k,j,i) = BufferLeft(i + (j-jbeg)*nx + (k-kbeg)*nx*ny + VsIndex );
                        Vs(BX3s,k+offset+1,j,i) = BufferRight(i + (j-jbeg)*nx + (k-kbeg)*nx*ny + VsIndex );
                    });

    #endif

#endif
    #ifdef MPI_NON_BLOCKING
        // Wait for the sends if they have not yet completed
        MPI_Waitall(2, sendRequest, sendStatus);
    #endif
    // Stop MPI Timer
    idfx::mpiTimer += timer.seconds();

#else
    IDEFIX_ERROR("You should not call DataBlock exchange routines without MPI");
#endif
    idfx::popRegion();
}
