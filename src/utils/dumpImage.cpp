// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2022 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include "dumpImage.hpp"
#include "output.hpp"
#include "idefix.hpp"


#define  HEADERSIZE 128

DumpImage::DumpImage(std::string filename, Output &output) {
  idfx::pushRegion("DumpImage::DumpImage");

  int nx[3];
  std::string fieldName;
  std::string eof ("eof");
  DataType type;
  int ndim;
  IdfxFileHandler fileHdl;
  Dump &dump = output.dump;

  idfx::cout << "DumpImage: loading restart file " << filename << "..." << std::flush;

  // open file
#ifdef WITH_MPI
  MPI_SAFE_CALL(MPI_File_open(MPI_COMM_WORLD, filename.c_str(),
                              MPI_MODE_RDONLY | MPI_MODE_UNIQUE_OPEN,
                              MPI_INFO_NULL, &fileHdl));
  dump.offset = 0;
#else
  fileHdl = fopen(filename.c_str(),"rb");
#endif

  // skip the header
#ifdef WITH_MPI
  dump.offset += HEADERSIZE;
#else
  fseek(fileHdl, HEADERSIZE, SEEK_SET);
#endif

  // First thing is to load the total domain size
  for(int dir=0 ; dir < 3; dir++) {
    dump.ReadNextFieldProperties(fileHdl, ndim, nx, type, fieldName);
    if(ndim>1) IDEFIX_ERROR("Wrong coordinate array dimensions while reading restart dump");
    // Store the size of the array
    this->np_int[dir] = nx[0];
    // Allocate arrays dynamically
    this->x[dir] = IdefixHostArray1D<real>("DumpImageX",np_int[dir]);
    this->xl[dir] = IdefixHostArray1D<real>("DumpImageXl",np_int[dir]);
    this->xr[dir] = IdefixHostArray1D<real>("DumpImageXr",np_int[dir]);

    // Read coordinates
    dump.ReadSerial(fileHdl, ndim, nx, type, reinterpret_cast<void*>( this->x[dir].data()) );
    dump.ReadNextFieldProperties(fileHdl, ndim, nx, type, fieldName);
    dump.ReadSerial(fileHdl, ndim, nx, type, reinterpret_cast<void*>( this->xl[dir].data()) );
    dump.ReadNextFieldProperties(fileHdl, ndim, nx, type, fieldName);
    dump.ReadSerial(fileHdl, ndim, nx, type, reinterpret_cast<void*>( this->xr[dir].data()) );
  }
    // Read the other fields

  while(true) {
    dump.ReadNextFieldProperties(fileHdl, ndim, nx, type, fieldName);
    if(fieldName.compare(eof) == 0) {
      break;
    } else if( ndim == 3) {
      // Load 3D field (raw data)
      // Make a new view of the right dimension for the raw data
      this->arrays[fieldName] = IdefixHostArray3D<real>("DumpImage"+fieldName,nx[2],nx[1],nx[0] );
      // Load it
      dump.ReadSerial(fileHdl,ndim,nx,type,
                      reinterpret_cast<void*>(this->arrays[fieldName].data()));
    } else if(fieldName.compare("time") == 0) {
      dump.ReadSerial(fileHdl, ndim, nx, type, &this->time);
    } else if(fieldName.compare("geometry")==0) {
      dump.ReadSerial(fileHdl, ndim, nx, type, &this->geometry);
    } else {
      int size=sizeof(double);
      for(int dim =0 ; dim<ndim ; dim++) {
        size = size*nx[dim];
      }
      double *scrch = new double[size];
      dump.ReadSerial(fileHdl,ndim, nx, type, scrch);
      delete[] scrch;
    }
  }
  // Close file

  #ifdef WITH_MPI
  MPI_SAFE_CALL(MPI_File_close(&fileHdl));
  #else
  fclose(fileHdl);
  #endif

  idfx::cout << "done." << std::endl;

  idfx::popRegion();
}
