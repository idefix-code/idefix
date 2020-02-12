#include "idefix.hpp"
#include "dataBlockHost.hpp"

DataBlockHost::DataBlockHost() {
    // Do nothing
}
DataBlockHost::DataBlockHost(DataBlock& datain) {
    
    // copy the dataBlock object for later use
    this->data=datain; 

    // Create mirrors (should be mirror_view)
    for(int dir = 0 ; dir < 3 ; dir++) {
        x[dir] = Kokkos::create_mirror_view(data.x[dir]);
        xr[dir] = Kokkos::create_mirror_view(data.xr[dir]);
        xl[dir] = Kokkos::create_mirror_view(data.xl[dir]);
        dx[dir] = Kokkos::create_mirror_view(data.dx[dir]);
        A[dir] = Kokkos::create_mirror_view(data.A[dir]);

    }
    dV = Kokkos::create_mirror_view(data.dV);
    Vc = Kokkos::create_mirror_view(data.Vc);
    Uc = Kokkos::create_mirror_view(data.Uc);

    // Store the grid informations from the dataBlock
    for(int dir = 0 ; dir < 3 ; dir++) {
        Kokkos::deep_copy(x[dir],data.x[dir]);
        Kokkos::deep_copy(xr[dir],data.xr[dir]);
        Kokkos::deep_copy(xl[dir],data.xl[dir]);
        Kokkos::deep_copy(dx[dir],data.dx[dir]);
        Kokkos::deep_copy(A[dir],data.A[dir]);
    }

    Kokkos::deep_copy(dV,data.dV);

}

// Synchronisation routines of Data (*Only*)
void DataBlockHost::SyncToDevice() {

    Kokkos::deep_copy(data.Vc,Vc);
    Kokkos::deep_copy(data.Uc,Uc);
}

void DataBlockHost::SyncFromDevice() {

    Kokkos::deep_copy(Vc,data.Vc);
    Kokkos::deep_copy(Uc,data.Uc);
}
