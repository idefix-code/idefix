#include "idefix.hpp"
#include "electromotiveforce.hpp"

ElectroMotiveForce::ElectroMotiveForce() {
    // Do nothing
}

// Init the emf from a datablock pointer
ElectroMotiveForce::ElectroMotiveForce(DataBlock *data) {
    Kokkos::Profiling::pushRegion("ElectroMotiveForce::Constructor(Datablock*)");

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

    Kokkos::Profiling::popRegion();

}