// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2021 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include "electroMotiveForce.hpp"
#include "dataBlock.hpp"
// Initialisation of electromotive force object
ElectroMotiveForce::ElectroMotiveForce() {
  // Do nothing
}

// Init the emf from a datablock pointer
ElectroMotiveForce::ElectroMotiveForce(DataBlock *data) {
  idfx::pushRegion("ElectroMotiveForce::Constructor(Datablock*)");

  #if MHD == YES
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

  EXPAND(   svx = IdefixArray3D<int>("EMF_svx",
                              data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);  ,
            svy = IdefixArray3D<int>("EMF_svy",
                              data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);  ,
            svz = IdefixArray3D<int>("EMF_svz",
                              data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);  )

  Ex1 = IdefixArray3D<real>("EMF_Ex1", data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
  Ex2 = IdefixArray3D<real>("EMF_Ex2", data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
  Ex3 = IdefixArray3D<real>("EMF_Ex3", data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);

  #endif
  idfx::popRegion();
}
