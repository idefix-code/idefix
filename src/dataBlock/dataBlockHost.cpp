#include "idefix.hpp"
#include "dataBlockHost.hpp"

DataBlockHost::DataBlockHost() {
    // Do nothing
}
DataBlockHost::DataBlockHost(DataBlock& datain) {
    
    idfx::pushRegion("DataBlockHost::DataBlockHost(DataBlock)");

    // copy the dataBlock object for later use
    this->data=&datain; 

    // By default, no current
    this->haveCurrent = false;
    
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
    Uc = Kokkos::create_mirror_view(data->Uc);
#if MHD == YES
    Vs = Kokkos::create_mirror_view(data->Vs);
    if(datain.haveCurrent) {
        this->haveCurrent = datain.haveCurrent;
        J = Kokkos::create_mirror_view(data->J);
    }
#endif
    

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
    if(this->haveCurrent && data->haveCurrent) Kokkos::deep_copy(data->J,J);
#endif
    Kokkos::deep_copy(data->Uc,Uc);

    idfx::popRegion();
}

void DataBlockHost::SyncFromDevice() {

    idfx::pushRegion("DataBlockHost::SyncFromDevice()");
    Kokkos::deep_copy(Vc,data->Vc);
#if MHD == YES
    Kokkos::deep_copy(Vs,data->Vs);
    if(this->haveCurrent && data->haveCurrent) Kokkos::deep_copy(J,data->J);
#endif
    Kokkos::deep_copy(Uc,data->Uc);

    idfx::popRegion();
}

void DataBlockHost::MakeVsFromAmag(IdefixHostArray4D<real> &Ain) {
    IdefixHostArray1D<real> dx1 = this->dx[IDIR];
    IdefixHostArray1D<real> dx2 = this->dx[JDIR];
    IdefixHostArray1D<real> dx3 = this->dx[KDIR];

    IdefixHostArray1D<real> x1m = this->xl[IDIR];
    IdefixHostArray1D<real> x2m = this->xl[JDIR];
    IdefixHostArray1D<real> x1 = this->x[IDIR];
    IdefixHostArray1D<real> x2 = this->x[JDIR];


    #if MHD == YES

    for(int k = data->beg[KDIR] ; k < data->end[KDIR] +KOFFSET ; k++) {
        for(int j = data->beg[JDIR] ; j < data->end[JDIR] +JOFFSET ; j++) {
            for(int i = data->beg[IDIR] ; i < data->end[IDIR] +IOFFSET; i++) {
                #if GEOMETRY == CARTESIAN
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
                #endif
                #if GEOMETRY == CYLINDRICAL
                    IDEFIX_ERROR("Not yet defined");
                #endif
                #if GEOMETRY == POLAR
                    Vs(BX1s,k,j,i) = D_EXPAND( ZERO_F                                     ,
                                                    - 1/(x1m(i)*dx2(j)) * (Ain(KDIR,k,j+1,i) - Ain(KDIR,k,j,i) )  ,
                                                    + 1/dx3(k) * (Ain(JDIR,k+1,j,i) - Ain(JDIR,k,j,i) ) );

                    Vs(BX2s,k,j,i) =  D_EXPAND(   1/dx1(i) * (Ain(KDIR,k,j,i+1) - Ain(KDIR,k,j,i) )  ,
                                                                                                    ,
                                                - 1/dx3(k) * (Ain(IDIR,k+1,j,i) - Ain(IDIR,k,j,i) ) );
                    
                    #if DIMENSIONS == 3
                    Vs(BX3s,k,j,i) =  - 1/(x1(i)*dx1(i)) * (x1m(i+1)*Ain(JDIR,k,j,i+1) - x1m(i)*Ain(JDIR,k,j,i) )  
                                    + 1/(x1(i)*dx2(j)) * (Ain(IDIR,k,j+1,i) - Ain(IDIR,k,j,i) ) ;
                    #endif

                #endif
                #if GEOMETRY == SPHERICAL
                    Vs(BX1s,k,j,i) = D_EXPAND( ZERO_F                                     ,
                                                    + 1/(x1m(i)*(cos(x2m(j))-cos(x2m(j+1)))) * (sin(x2m(j+1))*Ain(KDIR,k,j+1,i) - sin(x2m(j))*Ain(KDIR,k,j,i) )  ,
                                                    - 1/(x1m(i)*sin(x2(j))*dx3(k)) * (Ain(JDIR,k+1,j,i) - Ain(JDIR,k,j,i) ) );

                    Vs(BX2s,k,j,i) =  D_EXPAND(  - 1/(x1(i)*dx1(i)) * (x1m(i+1)*Ain(KDIR,k,j,i+1) - x1m(i)*Ain(KDIR,k,j,i) )  ,
                                                                                                    ,
                                                 + 1/(x1m(i)*sin(x2(j))*dx3(k)) * (Ain(IDIR,k+1,j,i) - Ain(IDIR,k,j,i) ) );
                    
                    #if DIMENSIONS == 3
                    Vs(BX3s,k,j,i) =  1/(x1(i)*dx1(i)) * (x1m(i+1)*Ain(JDIR,k,j,i+1) - x1m(i)*Ain(JDIR,k,j,i) )  
                                    - 1/(x1(i)*dx2(j)) * (Ain(IDIR,k,j+1,i) - Ain(IDIR,k,j,i) ) ;
                    #endif
                #endif
            }
        }
    }
    #else
    IDEFIX_ERROR("This function cannot be used without MHD enabled");
    #endif


}
