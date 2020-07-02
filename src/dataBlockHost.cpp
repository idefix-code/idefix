#include "idefix.hpp"
#include "dataBlockHost.hpp"

DataBlockHost::DataBlockHost() {
    // Do nothing
}
DataBlockHost::DataBlockHost(DataBlock& datain) {
    
    idfx::pushRegion("DataBlockHost::DataBlockHost(DataBlock)");

    // copy the dataBlock object for later use
    this->data=&datain; 

    // Create mirrors (should be mirror_view)
    for(int dir = 0 ; dir < 3 ; dir++) {
        x[dir] = Kokkos::create_mirror_view(data->x[dir]);
        xr[dir] = Kokkos::create_mirror_view(data->xr[dir]);
        xl[dir] = Kokkos::create_mirror_view(data->xl[dir]);
        dx[dir] = Kokkos::create_mirror_view(data->dx[dir]);
        A[dir] = Kokkos::create_mirror_view(data->A[dir]);
        np_tot[dir] = data->np_tot[dir];
        np_int[dir] = data->np_int[dir];
        nghost[dir] = data->nghost[dir];

        xbeg[dir] = data->xbeg[dir];
        xend[dir] = data->xend[dir];
        beg[dir] = data->beg[dir];
        end[dir] = data->end[dir];

        // TO BE COMPLETED...
    }

    dV = Kokkos::create_mirror_view(data->dV);
    Vc = Kokkos::create_mirror_view(data->Vc);
#if MHD == YES
    Vs = Kokkos::create_mirror_view(data->Vs);
#endif
    Uc = Kokkos::create_mirror_view(data->Uc);

    // Store the grid informations from the dataBlock
    for(int dir = 0 ; dir < 3 ; dir++) {
        Kokkos::deep_copy(x[dir],data->x[dir]);
        Kokkos::deep_copy(xr[dir],data->xr[dir]);
        Kokkos::deep_copy(xl[dir],data->xl[dir]);
        Kokkos::deep_copy(dx[dir],data->dx[dir]);
        Kokkos::deep_copy(A[dir],data->A[dir]);
    }

    Kokkos::deep_copy(dV,data->dV);

    idfx::popRegion();

}

// Synchronisation routines of Data (*Only*)
void DataBlockHost::SyncToDevice() {

    idfx::pushRegion("DataBlockHost::SyncToDevice()");

    Kokkos::deep_copy(data->Vc,Vc);
#if MHD == YES
    Kokkos::deep_copy(data->Vs,Vs);
#endif
    Kokkos::deep_copy(data->Uc,Uc);

    idfx::popRegion();
}

void DataBlockHost::SyncFromDevice() {

    idfx::pushRegion("DataBlockHost::SyncFromDevice()");
    Kokkos::deep_copy(Vc,data->Vc);
#if MHD == YES
    Kokkos::deep_copy(Vs,data->Vs);
#endif
    Kokkos::deep_copy(Uc,data->Uc);

    idfx::popRegion();
}

void DataBlockHost::MakeVsFromAmag(IdefixHostArray4D<real> &Ain) {
    IdefixHostArray1D<real> dx1 = this->dx[IDIR];
    IdefixHostArray1D<real> dx2 = this->dx[JDIR];
    IdefixHostArray1D<real> dx3 = this->dx[KDIR];

    #if MHD == YES

    for(int k = data->beg[KDIR] ; k < data->end[KDIR] +1 ; k++) {
        for(int j = data->beg[JDIR] ; j < data->end[JDIR] +1 ; j++) {
            for(int i = data->beg[IDIR] ; i < data->end[IDIR] +1; i++) {
                Vs(BX1s,k,j,i) = D_EXPAND( ZERO_F                                     ,
                                                - 1/dx2(j) * (Ain(KDIR,k,j+1,i) - Ain(KDIR,k,j,i) )  ,
                                                + 1/dx3(k) * (Ain(JDIR,k+1,j,i) - Ain(JDIR,k,j,i) ) );
                #if DIMENSIONS >= 2
                Vs(BX2s,k,j,i) =  D_EXPAND(   1/dx1(i) * (Ain(KDIR,k,j,i+1) - Ain(KDIR,k,j,i) )  ,
                                                                                                ,
                                            - 1/dx3(k) * (Ain(IDIR,k+1,j,i) - Ain(IDIR,k,j,i) ) );
                #endif
                #if DIMENSIONS == 3
                Vs(BX3s,k,j,i) =  - 1/dx1(i) * (Ain(JDIR,k,j,i+1) - Ain(JDIR,k,j,i) )  
                                  + 1/dx2(j) * (Ain(IDIR,k,j+1,i) - Ain(IDIR,k,j,i) ) ;
                #endif
            }
        }
    }
    #else
    IDEFIX_ERROR("This function cannot be used without MHD enabled");
    #endif


}
