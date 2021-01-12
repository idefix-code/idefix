// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include "../idefix.hpp"
#include "dataBlock.hpp"

// dump the current dataBlock to a file (mainly used for debug purposes)
void DataBlock::DumpToFile(std::string filebase)  {
  FILE *fileHdl;

  real nfield;
  real nx1;
  real nx2;
  real nx3;
  real nv;
  static int n=0;

  IdefixArray4D<real>::HostMirror locVc = Kokkos::create_mirror_view(this->hydro.Vc);
  Kokkos::deep_copy(locVc,this->hydro.Vc);

  // TODO(lesurg) Make datablock a friend of hydro to get the Riemann flux?
  //IdefixArray4D<real>::HostMirror locFlux = Kokkos::create_mirror_view(this->hydro.FluxRiemann);
  //Kokkos::deep_copy(locFlux, this->FluxRiemann);
#if MHD == YES
  IdefixArray4D<real>::HostMirror locVs = Kokkos::create_mirror_view(this->hydro.Vs);
  Kokkos::deep_copy(locVs,this->hydro.Vs);

  IdefixArray4D<real>::HostMirror locJ;
  if(hydro.haveCurrent) {
    locJ = Kokkos::create_mirror_view(this->hydro.J);
    Kokkos::deep_copy(locJ, this->hydro.J);
  }
#endif

  // Data format:
  // First groupe is nfield (number of 4D arrays written)
  // Then for each field:
  // nx1, nx2, nx3, nvar then the data, then off you go
  std::string dot = std::string(".");
  std::string ext = std::string("idfx");
  std::string filename = filebase + dot + std::to_string(n)
                        + dot + std::to_string(idfx::prank) + dot + ext;
  n++;
  fileHdl = fopen(filename.c_str(),"wb");

#if MHD== YES
  nfield = 2;
  if(hydro.haveCurrent) nfield++;
#else
  nfield = 1;
#endif

  fwrite(&nfield, sizeof(real), 1, fileHdl);

  // Write Vc
  nx1=this->np_tot[IDIR];
  nx2=this->np_tot[JDIR];
  nx3=this->np_tot[KDIR];
  nv=NVAR;

  fwrite(&nx1, sizeof(real),1,fileHdl);
  fwrite(&nx2, sizeof(real),1,fileHdl);
  fwrite(&nx3, sizeof(real),1,fileHdl);
  fwrite(&nv, sizeof(real),1,fileHdl);

  fwrite(locVc.data(), sizeof(real), nx1*nx2*nx3*nv, fileHdl);

  // Write Flux
  /*
  nx1=this->np_tot[IDIR];
  nx2=this->np_tot[JDIR];
  nx3=this->np_tot[KDIR];
  nv=NVAR;

  fwrite(&nx1, sizeof(real),1,fileHdl);
  fwrite(&nx2, sizeof(real),1,fileHdl);
  fwrite(&nx3, sizeof(real),1,fileHdl);
  fwrite(&nv, sizeof(real),1,fileHdl);

  fwrite(locFlux.data(), sizeof(real), nx1*nx2*nx3*nv, fileHdl);
  */

  // Write Vs
#if MHD == YES
  nx1=nx1+IOFFSET;
  nx2=nx2+JOFFSET;
  nx3=nx3+KOFFSET;
  nv=DIMENSIONS;

  fwrite(&nx1, sizeof(real),1,fileHdl);
  fwrite(&nx2, sizeof(real),1,fileHdl);
  fwrite(&nx3, sizeof(real),1,fileHdl);
  fwrite(&nv, sizeof(real),1,fileHdl);

  fwrite(locVs.data(), sizeof(real), nx1*nx2*nx3*nv, fileHdl);

  if(hydro.haveCurrent) {
    nx1=this->np_tot[IDIR];
    nx2=this->np_tot[JDIR];
    nx3=this->np_tot[KDIR];
    nv=3;

    fwrite(&nx1, sizeof(real),1,fileHdl);
    fwrite(&nx2, sizeof(real),1,fileHdl);
    fwrite(&nx3, sizeof(real),1,fileHdl);
    fwrite(&nv, sizeof(real),1,fileHdl);

    fwrite(locJ.data(), sizeof(real), nx1*nx2*nx3*nv, fileHdl);
  }

#endif

  fclose(fileHdl);
}
