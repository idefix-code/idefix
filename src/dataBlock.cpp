#include "idefix.hpp"
#include "dataBlock.hpp"

DataBlock::DataBlock() {
    // Do nothing
}

void DataBlock::InitFromGrid(Grid &grid) {
    // This initialisation is only valid for *serial*
    // MPI initialisation will involve domain decomposition of grids into DataBlocks

    Kokkos::Profiling::pushRegion("DataBlock::InitFromGrid");

    // Copy the number of points from grid since DataBlock=Grid in serial
    for(int dir = 0 ; dir < 3 ; dir++) {
        nghost[dir] = grid.nghost[dir];
        np_tot[dir] = grid.np_tot[dir];
        np_int[dir] = grid.np_int[dir];

        // Boundary conditions are yet to be defined
        lbound[dir] = grid.lbound[dir]; 
        rbound[dir] = grid.rbound[dir];

        beg[dir] = grid.nghost[dir];
        end[dir] = grid.nghost[dir]+np_int[dir];

        // Where does this datablock starts and end in the grid?
        gbeg[dir] = beg[dir];
        gend[dir] = end[dir];

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

    Kokkos::Profiling::popRegion();
    // TODO: A proper initialisation of A and dV should be done at this stage
}

DataBlock::DataBlock(const DataBlock &data) {

    Kokkos::Profiling::pushRegion("DataBlock::DataBlock(copy)");
    
    // Shallow copy of a datablock
    for(int dir = 0 ; dir < 3; dir++) {
        nghost[dir]=data.nghost[dir];
        np_tot[dir]=data.np_tot[dir];
        np_int[dir]=data.np_int[dir];
        lbound[dir]=data.lbound[dir];
        rbound[dir]=data.rbound[dir];
        beg[dir]=data.beg[dir];
        end[dir]=data.end[dir];
        gbeg[dir]=data.gbeg[dir];
        gend[dir]=data.gend[dir];
        x[dir]=data.x[dir];
        xr[dir]=data.xr[dir];
        xl[dir]=data.xl[dir];
        dx[dir]=data.dx[dir];
        A[dir]=data.A[dir];
    }

    dV=data.dV;
    Vc=data.Vc;
    Vc0=data.Vc0;
    Vs0=data.Vs0;
    Uc=data.Uc;
    InvDtHyp=data.InvDtHyp;
    InvDtPar=data.InvDtPar;

    #if MHD == YES
    emf=data.emf;
    Vs=data.Vs;
    #endif

    Kokkos::Profiling::popRegion();
}


DataBlock& DataBlock::operator=(const DataBlock& data) {
    Kokkos::Profiling::pushRegion("DataBlock::operator=");
    // Shallow reference operator
    if (this != & data) {
        for(int dir = 0 ; dir < 3; dir++) {
            nghost[dir]=data.nghost[dir];
            np_tot[dir]=data.np_tot[dir];
            np_int[dir]=data.np_int[dir];
            lbound[dir]=data.lbound[dir];
            rbound[dir]=data.rbound[dir];
            beg[dir]=data.beg[dir];
            end[dir]=data.end[dir];
            gbeg[dir]=data.gbeg[dir];
            gend[dir]=data.gend[dir];
            x[dir]=data.x[dir];
            xr[dir]=data.xr[dir];
            xl[dir]=data.xl[dir];
            dx[dir]=data.dx[dir];
            A[dir]=data.A[dir];
        }
        
        dV=data.dV;
        Vc=data.Vc;
        Vc0=data.Vc0;
        Vs0=data.Vs0;
        Uc=data.Uc;
        InvDtHyp=data.InvDtHyp;
        InvDtPar=data.InvDtPar;
        
        #if MHD == YES
        emf=data.emf;
        Vs=data.Vs;
        #endif
    }

    Kokkos::Profiling::popRegion();

    return *this;
}
