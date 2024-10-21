// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include <string>
#include <cstdio>
#include "../idefix.hpp"
#include "dataBlock.hpp"
#include "version.hpp"
#include "fluid.hpp"

#define  NAMESIZE     16
#define  HEADERSIZE 128

void DataBlock::WriteVariable(FILE* fileHdl, int ndim, int *dim,
                                char *name, void* data) {
  int ntot = 1;   // Number of elements to be written
  int size = sizeof(real);  // Elements always assumed to be real

  // Write field name
  fwrite (name, sizeof(char), NAMESIZE, fileHdl);

  // Write dimensions of array
  fwrite(&ndim, 1, sizeof(int), fileHdl);
  for(int n = 0 ; n < ndim ; n++) {
    fwrite(dim+n, 1, sizeof(int), fileHdl);
    ntot = ntot * dim[n];
  }
  // Write raw data
  fwrite(data, ntot, size, fileHdl);
}

// dump the current dataBlock to a file (mainly used for debug purposes)
void DataBlock::DumpToFile(std::string filebase)  {
  FILE *fileHdl;

  int dims[5];
  char fieldName[NAMESIZE+1]; // +1 is just in case

  static int n=0;



  // TODO(lesurg) Make datablock a friend of hydro to get the Riemann flux?
  //IdefixArray4D<real>::HostMirror locFlux = Kokkos::create_mirror_view(this->hydro->FluxRiemann);
  //Kokkos::deep_copy(locFlux, this->FluxRiemann);
#if MHD == YES


  IdefixArray4D<real>::HostMirror locJ;
  if(hydro->haveCurrent) {
    locJ = Kokkos::create_mirror_view(this->hydro->J);
    Kokkos::deep_copy(locJ, this->hydro->J);
  }
#endif

  // Data format:
  // Similar to dump files, but no aggregation through cores
  // and include boundaries
  // Header [HEADER SIZE]
  // Field Name [NAMESIZE]
  // # of dimensions  (integet)
  // dimensions (n*integer)
  // raw data (real)
  // FieldName [NAMESIZE]...
  // FieldName [NAMESIZE] = "eof"

  std::string dot = std::string(".");
  std::string ext = std::string("idfx");
  std::string filename = filebase + dot + std::to_string(n)
                        + dot + std::to_string(idfx::prank) + dot + ext;
  n++;
  fileHdl = fopen(filename.c_str(),"wb");

  // Write Header
  char header[HEADERSIZE];
  std::snprintf(header, HEADERSIZE, "Idefix %s Debug DataBlock", IDEFIX_VERSION);
  fwrite (header, sizeof(char), HEADERSIZE, fileHdl);

  // Write Vc
  IdefixArray4D<real>::HostMirror locVc = Kokkos::create_mirror_view(this->hydro->Vc);
  Kokkos::deep_copy(locVc,this->hydro->Vc);
  dims[0] = this->np_tot[IDIR];
  dims[1] = this->np_tot[JDIR];
  dims[2] = this->np_tot[KDIR];
  dims[3] = NVAR;

  std::snprintf(fieldName,NAMESIZE,"Vc");

  WriteVariable(fileHdl, 4, dims, fieldName, locVc.data());

  if (this->gravity->haveSelfGravityPotential) {
    IdefixArray3D<real>::HostMirror locPot = Kokkos::create_mirror_view(this->gravity->phiP);
    Kokkos::deep_copy(locPot, this->gravity->phiP);

    dims[3] = 1;
    std::snprintf(fieldName,NAMESIZE,"Pot");

    WriteVariable(fileHdl, 4, dims, fieldName, locPot.data());
  }

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
  // Write Vs
  IdefixArray4D<real>::HostMirror locVs = Kokkos::create_mirror_view(this->hydro->Vs);
  Kokkos::deep_copy(locVs,this->hydro->Vs);
  dims[0] = this->np_tot[IDIR]+IOFFSET;
  dims[1] = this->np_tot[JDIR]+JOFFSET;
  dims[2] = this->np_tot[KDIR]+KOFFSET;
  dims[3] = DIMENSIONS;

  std::snprintf(fieldName,NAMESIZE,"Vs");

  WriteVariable(fileHdl, 4, dims, fieldName, locVs.data());

  // Write EMFs
  dims[0] = this->np_tot[IDIR];
  dims[1] = this->np_tot[JDIR];
  dims[2] = this->np_tot[KDIR];

  std::snprintf(fieldName,NAMESIZE,"Ex3");
  IdefixArray3D<real>::HostMirror locE = Kokkos::create_mirror_view(this->hydro->emf->ez);
  Kokkos::deep_copy(locE,this->hydro->emf->ez);
  WriteVariable(fileHdl, 3, dims, fieldName, locE.data());

  #if DIMENSIONS == 3
  std::snprintf(fieldName,NAMESIZE,"Ex1");
  Kokkos::deep_copy(locE,this->hydro->emf->ex);
  WriteVariable(fileHdl, 3, dims, fieldName, locE.data());
  std::snprintf(fieldName,NAMESIZE,"Ex2");
  Kokkos::deep_copy(locE,this->hydro->emf->ey);
  WriteVariable(fileHdl, 3, dims, fieldName, locE.data());
  #endif


  if(hydro->haveCurrent) {
    IdefixArray4D<real>::HostMirror locJ = Kokkos::create_mirror_view(this->hydro->J);
    Kokkos::deep_copy(locJ,this->hydro->J);
    dims[0] = this->np_tot[IDIR];
    dims[1] = this->np_tot[JDIR];
    dims[2] = this->np_tot[KDIR];
    dims[3] = 3;

    std::snprintf(fieldName,NAMESIZE,"J");

    WriteVariable(fileHdl, 4, dims, fieldName, locJ.data());
  }

#endif

  // Write end of file
  std::snprintf(fieldName,NAMESIZE,"eof");
  dims[0] = 1;
  WriteVariable(fileHdl, 1, dims, fieldName, locVc.data());

  fclose(fileHdl);
}
