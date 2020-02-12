#include "idefix.hpp"
#include "dataBlock.hpp"

DataBlock::DataBlock() {
    // Do nothing
}

void DataBlock::InitFromGrid(Grid &grid) {
    // This initialisation is only valid for *serial*
    // MPI initialisation will involve domain decomposition of grids into DataBlocks

    
    // Copy the number of points from grid since DataBlock=Grid in serial
    for(int dir = 0 ; dir < 3 ; dir++) {
        nghost[dir] = grid.nghost[dir];
        np_tot[dir] = grid.np_tot[dir];
        np_int[dir] = grid.np_int[dir];

        // Boundary conditions are yet to be defined
        lbound[dir] = 0; 
        rbound[dir] = 0;

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
    Vc = IdefixArray4D<real>("DataBlock_Vc", NVAR, np_tot[KDIR], grid.np_tot[JDIR], grid.np_tot[IDIR]);
    Uc = IdefixArray4D<real>("DataBlock_Uc", NVAR, np_tot[KDIR], grid.np_tot[JDIR], grid.np_tot[IDIR]);

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
    
    // TODO: A proper initialisation of A and dV should be done at this stage
}

DataBlock::DataBlock(const DataBlock &data) {
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
    Uc=data.Uc;
}

DataBlock& DataBlock::operator=(const DataBlock& data) {
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
        Uc=data.Uc;
    }

    return *this;
}