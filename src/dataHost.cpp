#include "idefix.hpp"
#include "dataHost.hpp"

DataHost::DataHost() {
    // Do nothing
}
DataHost::DataHost(Data& datain) {
    std::cout << "Building dataHost Array\n";
    
    // copy the data object for later use
    this->data=datain; 

    // Create mirrors (should be mirror_view)
    Vc = Kokkos::create_mirror(data.Vc);
    Uc = Kokkos::create_mirror(data.Uc);

}

// Synchronisation routines
void DataHost::SyncToDevice() {
    Kokkos::deep_copy(data.Vc,Vc);
    Kokkos::deep_copy(data.Uc,Uc);
}

void DataHost::SyncFromDevice() {
    Kokkos::deep_copy(Vc,data.Vc);
    Kokkos::deep_copy(Uc,data.Uc);
}
