#include "idefix.hpp"
#include "dataBlock.hpp"

DataBlock::DataBlock() {
    // Do nothing
}

void DataBlock::InitFromGrid(Grid &grid) {
    // This initialisation is only valid for *serial*
    // MPI initialisation will involve domain decomposition of grids into DataBlocks

    Kokkos::Profiling::pushRegion("DataBlock::InitFromGrid");

    this->mygrid=&grid;

    // Make a local copy of the grid for future usage.
    GridHost gridHost(grid);

    // Copy the number of points from grid since DataBlock=Grid in serial
    for(int dir = 0 ; dir < 3 ; dir++) {

        // Check that the dimension is effectively divisible by number of procs
        if(grid.np_int[dir] % grid.nproc[dir]) {
            IDEFIX_ERROR("Grid size must be a multiple of the domain decomposition");
        }

        nghost[dir] = grid.nghost[dir];
        np_int[dir] = grid.np_int[dir]/grid.nproc[dir];
        np_tot[dir] = np_int[dir]+2*nghost[dir];
        

        // Boundary conditions
        if(grid.xproc[dir]==0) lbound[dir] = grid.lbound[dir];
        else lbound[dir] = internal;

        if(grid.xproc[dir] == grid.nproc[dir]-1) rbound[dir] = grid.rbound[dir];
        else rbound[dir] = internal;

        beg[dir] = grid.nghost[dir];
        end[dir] = grid.nghost[dir]+np_int[dir];

        // Where does this datablock starts and end in the grid?
        gbeg[dir] = grid.nghost[dir] + grid.xproc[dir]*np_int[dir];  // This assumes even distribution of points between procs
        gend[dir] = grid.nghost[dir] + (grid.xproc[dir]+1)*np_int[dir];

        // Local start and end of current datablock
        xbeg[dir] = gridHost.xl[dir](gbeg[dir]);
        xend[dir] = gridHost.xl[dir](gend[dir]);

    }

    if(idfx::psize>1) {
        idfx::cout << "DataBlock::initFromGrid local size is " << std::endl;
        
        for(int dir = 0 ; dir < DIMENSIONS ; dir++) {
            idfx::cout << "\t Direction X" << (dir+1) << ": " << xbeg[dir] << "...." << np_int[dir] << "...." << xend[dir] << std::endl;
        }
    }
    
    // Allocate the required fields
    for(int dir = 0 ; dir < 3 ; dir++) {
        x[dir] = IdefixArray1D<real>("DataBlock_x",np_tot[dir]);
        xr[dir] = IdefixArray1D<real>("DataBlock_xr",np_tot[dir]);
        xl[dir] = IdefixArray1D<real>("DataBlock_xl",np_tot[dir]);
        dx[dir] = IdefixArray1D<real>("DataBlock_dx",np_tot[dir]);
        
        A[dir] = IdefixArray3D<real>("DataBlock_A",np_tot[KDIR],np_tot[JDIR],np_tot[IDIR]);
    }

    dV = IdefixArray3D<real>("DataBlock_dV",np_tot[KDIR],np_tot[JDIR],np_tot[IDIR]);
    Vc = IdefixArray4D<real>("DataBlock_Vc", NVAR, np_tot[KDIR], np_tot[JDIR], np_tot[IDIR]);
    Uc = IdefixArray4D<real>("DataBlock_Uc", NVAR, np_tot[KDIR], np_tot[JDIR], np_tot[IDIR]);
    Vc0 = IdefixArray4D<real>("DataBlock_Vc0", NVAR, np_tot[KDIR], np_tot[JDIR], np_tot[IDIR]);
#if MHD == YES

    Vs = IdefixArray4D<real>("DataBlock_Vs", DIMENSIONS, np_tot[KDIR]+KOFFSET, np_tot[JDIR]+JOFFSET, np_tot[IDIR]+IOFFSET);
    Vs0 = IdefixArray4D<real>("DataBlock_Vs0", DIMENSIONS, np_tot[KDIR]+KOFFSET, np_tot[JDIR]+JOFFSET, np_tot[IDIR]+IOFFSET);

    this->emf = ElectroMotiveForce(this);
#endif

    InvDtHyp = IdefixArray3D<real>("DataBlock_InvDtHyp", np_tot[KDIR], np_tot[JDIR], np_tot[IDIR]);
    InvDtPar = IdefixArray3D<real>("DataBlock_InvDtPar", np_tot[KDIR], np_tot[JDIR], np_tot[IDIR]);
    PrimL =  IdefixArray4D<real>("DataBlock_PrimL", NVAR, np_tot[KDIR], np_tot[JDIR], np_tot[IDIR]);
    PrimR =  IdefixArray4D<real>("DataBlock_PrimR", NVAR, np_tot[KDIR], np_tot[JDIR], np_tot[IDIR]);
    FluxRiemann =  IdefixArray4D<real>("DataBlock_FluxRiemann", NVAR, np_tot[KDIR], np_tot[JDIR], np_tot[IDIR]);

// Init MPI Buffer Arrays
#ifdef WITH_MPI
    bufferSizeX1 = 0;
    bufferSizeX2 = 0;
    bufferSizeX3 = 0;

    // Number of cells in X1 boundary condition:
    bufferSizeX1 = nghost[IDIR] * np_int[JDIR] * np_int[KDIR] * NVAR;

    #if MHD == YES
    // BX1s
    bufferSizeX1 += nghost[IDIR] * np_int[JDIR] * np_int[KDIR];
    #if DIMENSIONS>=2
    bufferSizeX1 += nghost[IDIR] * (np_int[JDIR]+1) * np_int[KDIR];
    #endif
    #if DIMENSIONS==3
    bufferSizeX1 += nghost[IDIR] * np_int[JDIR] * (np_int[KDIR]+1);
    #endif  // DIMENSIONS
    #endif  // MHD

    BufferRecvX1[faceLeft ] = IdefixArray1D<real>("BufferRecvX1Left", bufferSizeX1);
    BufferRecvX1[faceRight] = IdefixArray1D<real>("BufferRecvX1Right",bufferSizeX1);
    BufferSendX1[faceLeft ] = IdefixArray1D<real>("BufferSendX1Left", bufferSizeX1);
    BufferSendX1[faceRight] = IdefixArray1D<real>("BufferSendX1Right",bufferSizeX1);

    // Number of cells in X2 boundary condition (only required when problem >2D):
#if DIMENSIONS >= 2 
    bufferSizeX2 = np_tot[IDIR] * nghost[JDIR] * np_int[KDIR] * NVAR;
    #if MHD == YES
    // BX1s
    bufferSizeX2 += (np_tot[IDIR]+1) * nghost[JDIR] * np_int[KDIR];
    // BX2s
    bufferSizeX2 += np_tot[IDIR] * nghost[JDIR] * np_int[KDIR];
    #if DIMENSIONS==3
    bufferSizeX2 += np_tot[IDIR] * nghost[JDIR] * (np_int[KDIR]+1);
    #endif  // DIMENSIONS
    #endif  // MHD

    BufferRecvX2[faceLeft ] = IdefixArray1D<real>("BufferRecvX2Left", bufferSizeX2);
    BufferRecvX2[faceRight] = IdefixArray1D<real>("BufferRecvX2Right",bufferSizeX2);
    BufferSendX2[faceLeft ] = IdefixArray1D<real>("BufferSendX2Left", bufferSizeX2);
    BufferSendX2[faceRight] = IdefixArray1D<real>("BufferSendX2Right",bufferSizeX2);

#endif
// Number of cells in X3 boundary condition (only required when problem is 3D):
#if DIMENSIONS ==3  
    bufferSizeX3 = np_tot[IDIR] * np_tot[JDIR] * nghost[KDIR] * NVAR;

    #if MHD == YES
    // BX1s
    bufferSizeX3 += (np_tot[IDIR]+1) * np_tot[JDIR] * nghost[KDIR];
    // BX2s
    bufferSizeX3 += np_tot[IDIR] * (np_tot[JDIR]+1) * nghost[KDIR];
    // BX3s (not needed because reconstructed)
    bufferSizeX3 += np_tot[IDIR] * np_tot[JDIR] * nghost[KDIR];
    #endif  // DIMENSIONS
    #endif  // MHD

    BufferRecvX3[faceLeft ] = IdefixArray1D<real>("BufferRecvX3Left", bufferSizeX3);
    BufferRecvX3[faceRight] = IdefixArray1D<real>("BufferRecvX3Right",bufferSizeX3);
    BufferSendX3[faceLeft ] = IdefixArray1D<real>("BufferSendX3Left", bufferSizeX3);
    BufferSendX3[faceRight] = IdefixArray1D<real>("BufferSendX3Right",bufferSizeX3);
#endif


    // We need to copy the relevant part of the coordinate system to the datablock
    for(int dir = 0 ; dir < 3 ; dir++) {
        int offset=gbeg[dir]-beg[dir];

        IdefixArray1D<real> x_input = grid.x[dir];
        IdefixArray1D<real> x_output= x[dir];
        IdefixArray1D<real> xr_input = grid.xr[dir];
        IdefixArray1D<real> xr_output= xr[dir];
        IdefixArray1D<real> xl_input = grid.xl[dir];
        IdefixArray1D<real> xl_output= xl[dir];
        IdefixArray1D<real> dx_input = grid.dx[dir];
        IdefixArray1D<real> dx_output= dx[dir];
        

        idefix_for("coordinates",0,np_tot[dir],
                        KOKKOS_LAMBDA (int i) {
                            x_output(i) = x_input(i+offset);
                            xr_output(i) = xr_input(i+offset);
                            xl_output(i) = xl_input(i+offset);
                            dx_output(i) = dx_input(i+offset); 
                        });
    }


    // Fill the names of the fields
    for(int i = 0 ; i < NVAR ;  i++) {
        switch(i) {
            case RHO:
                VcName.push_back("RHO");
                break;
            case VX1:
                VcName.push_back("VX1");
                break;
            case VX2:
                VcName.push_back("VX2");
                break;
            case VX3:
                VcName.push_back("VX3");
                break;
            case BX1:
                VcName.push_back("BX1");
                break;
            case BX2:
                VcName.push_back("BX2");
                break;
            case BX3:
                VcName.push_back("BX3");
                break;
            #if HAVE_ENERGY
            case PRS:
                VcName.push_back("PRS");
                break;
            #endif
            default:
                VcName.push_back("Vc_"+std::to_string(i));
        }
    }

    for(int i = 0 ; i < DIMENSIONS ; i++) {
        switch(i) {
            case 0:
                VsName.push_back("BX1s");
                break;
            case 1:
                VsName.push_back("BX2s");
                break;
            case 2:
                VsName.push_back("BX3s");
                break;
            default:
                VsName.push_back("Vs_"+std::to_string(i));
        }
    }

    Kokkos::Profiling::popRegion();
    // TODO: A proper initialisation of A and dV should be done at this stage
}

// MPI Routines exchange
void DataBlock::ExchangeAll(){
#ifdef WITH_MPI

#else
    IDEFIX_ERROR("You should not call DataBlock exchange routines without MPI");
#endif
}

void DataBlock::ExchangeX1() {
    Kokkos::Profiling::pushRegion("DataBlock::ExchangeX1");
#ifdef WITH_MPI
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

    // Send to the right
    int procSend, procRecv;
    MPI_Status status;
    MPI_SAFE_CALL(MPI_Cart_shift(mygrid->CartComm,0,1,&procRecv,&procSend ));   // We receive from procRecv, and we send to procSend
    Kokkos::fence();
    MPI_SAFE_CALL(MPI_Sendrecv(BufferSendX1[faceRight].data(), bufferSizeX1, realMPI, procSend, 100,
                 BufferRecvX1[faceLeft].data(), bufferSizeX1, realMPI, procRecv, 100,
                 mygrid->CartComm, &status));
    
    // Send to the left
    MPI_SAFE_CALL(MPI_Cart_shift(mygrid->CartComm,0,-1,&procRecv,&procSend ));   // We receive from procRecv, and we send to procSend
    MPI_SAFE_CALL(MPI_Sendrecv(BufferSendX1[faceLeft].data(), bufferSizeX1, realMPI, procSend, 101,
                 BufferRecvX1[faceRight].data(), bufferSizeX1, realMPI, procRecv, 101,
                 mygrid->CartComm, &status));

    
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


#else
    IDEFIX_ERROR("You should not call DataBlock exchange routines without MPI");
#endif
    Kokkos::Profiling::popRegion();
}


void DataBlock::ExchangeX2() {
    Kokkos::Profiling::pushRegion("DataBlock::ExchangeX2");
#ifdef WITH_MPI
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
    MPI_Status status;
    MPI_SAFE_CALL(MPI_Cart_shift(mygrid->CartComm,1,1,&procRecv,&procSend ));   // We receive from procRecv, and we send to procSend
    Kokkos::fence();
    MPI_SAFE_CALL(MPI_Sendrecv(BufferSendX2[faceRight].data(), bufferSizeX2, realMPI, procSend, 200,
                 BufferRecvX2[faceLeft].data(), bufferSizeX2, realMPI, procRecv, 200,
                 mygrid->CartComm, &status));
    
    // Send to the left
    MPI_SAFE_CALL(MPI_Cart_shift(mygrid->CartComm,1,-1,&procRecv,&procSend ));   // We receive from procRecv, and we send to procSend
    MPI_SAFE_CALL(MPI_Sendrecv(BufferSendX2[faceLeft].data(), bufferSizeX2, realMPI, procSend, 201,
                 BufferRecvX2[faceRight].data(), bufferSizeX2, realMPI, procRecv, 201,
                 mygrid->CartComm, &status));

    
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


#else
    IDEFIX_ERROR("You should not call DataBlock exchange routines without MPI");
#endif
    Kokkos::Profiling::popRegion();
}


void DataBlock::ExchangeX3(){
    Kokkos::Profiling::pushRegion("DataBlock::ExchangeX3");
#ifdef WITH_MPI
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
    MPI_Status status;
    MPI_SAFE_CALL(MPI_Cart_shift(mygrid->CartComm,2,1,&procRecv,&procSend ));   // We receive from procRecv, and we send to procSend
    Kokkos::fence();
    MPI_SAFE_CALL(MPI_Sendrecv(BufferSendX3[faceRight].data(), bufferSizeX3, realMPI, procSend, 300,
                 BufferRecvX3[faceLeft].data(), bufferSizeX3, realMPI, procRecv, 300,
                 mygrid->CartComm, &status));
    
    // Send to the left
    MPI_SAFE_CALL(MPI_Cart_shift(mygrid->CartComm,2,-1,&procRecv,&procSend ));   // We receive from procRecv, and we send to procSend
    MPI_SAFE_CALL(MPI_Sendrecv(BufferSendX3[faceLeft].data(), bufferSizeX3, realMPI, procSend, 301,
                 BufferRecvX3[faceRight].data(), bufferSizeX3, realMPI, procRecv, 301,
                 mygrid->CartComm, &status));

    
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


#else
    IDEFIX_ERROR("You should not call DataBlock exchange routines without MPI");
#endif
    Kokkos::Profiling::popRegion();
}


// Initialisation of electromotive force object
ElectroMotiveForce::ElectroMotiveForce() {
    // Do nothing
}

// Init the emf from a datablock pointer
ElectroMotiveForce::ElectroMotiveForce(DataBlock *data) {
    Kokkos::Profiling::pushRegion("ElectroMotiveForce::Constructor(Datablock*)");

    #if MHD == YES
    D_EXPAND( ez = IdefixArray3D<real>("EMF_ez", data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]) ;     ,
                                                                                                             ;     ,
              ex = IdefixArray3D<real>("EMF_ex", data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]) ;
              ey = IdefixArray3D<real>("EMF_ey", data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]) ; )
    
    D_EXPAND( ezi = IdefixArray3D<real>("EMF_ezi", data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]) ;
              ezj = IdefixArray3D<real>("EMF_ezj", data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]) ;     ,
                                                                                                               ;     ,
              exj = IdefixArray3D<real>("EMF_exj", data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]) ;
              exk = IdefixArray3D<real>("EMF_exj", data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]) ;
              eyi = IdefixArray3D<real>("EMF_eyi", data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]) ;
              eyk = IdefixArray3D<real>("EMF_eyi", data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]) ; )

    #endif
    Kokkos::Profiling::popRegion();

}
