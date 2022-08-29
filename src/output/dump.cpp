// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2022 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include "dump.hpp"
#include "gitversion.hpp"
#include "dataBlockHost.hpp"
#include "gridHost.hpp"
#include "output.hpp"

// Max size of array name
#define  NAMESIZE     16
#define  FILENAMESIZE   256
#define  HEADERSIZE 128

void Dump::Init(Input &input, DataBlock &data) {
  // Init the output period


  for (int dir=0; dir<3; dir++) {
    this->periodicity[dir] = (data.mygrid->lbound[dir] == periodic);
  }
  this->dumpFileNumber = 0;

  // Allocate scratch Array
  this->scrch = new real[ (data.np_int[IDIR]+IOFFSET)*
                          (data.np_int[JDIR]+JOFFSET)*
                          (data.np_int[KDIR]+KOFFSET)];

  #ifdef WITH_MPI
    Grid *grid = data.mygrid;

    int start[3];
    int size[3];
    int subsize[3];

    // Dimensions for cell-centered fields
    for(int dir = 0; dir < 3 ; dir++) {
      size[2-dir] = grid->np_int[dir];
      start[2-dir] = data.gbeg[dir]-data.nghost[dir];
      subsize[2-dir] = data.np_int[dir];
    }

    MPI_SAFE_CALL(MPI_Type_create_subarray(3, size, subsize, start,
                                           MPI_ORDER_C, realMPI, &this->descC));
    MPI_SAFE_CALL(MPI_Type_commit(&this->descC));

    // Dimensions for face-centered field
    for(int face = 0; face < 3 ; face++) {
      for(int dir = 0; dir < 3 ; dir++) {
        size[2-dir] = grid->np_int[dir];
        start[2-dir] = data.gbeg[dir]-data.nghost[dir];
        subsize[2-dir] = data.np_int[dir];
      }
      // Add the extra guy in the face direction
      size[2-face]++;
      subsize[2-face]++; // valid only for reading
                         //since it involves an overlap of data between procs

      MPI_SAFE_CALL(MPI_Type_create_subarray(3, size, subsize, start,
                                             MPI_ORDER_C, realMPI, &this->descSR[face]));
      MPI_SAFE_CALL(MPI_Type_commit(&this->descSR[face]));

      // Now for writing, it is only the last proc which keeps one additional cell
      if(grid->xproc[face] != grid->nproc[face] - 1  ) subsize[2-face]--;
      MPI_SAFE_CALL(MPI_Type_create_subarray(3, size, subsize, start,
                                             MPI_ORDER_C, realMPI, &this->descSW[face]));
      MPI_SAFE_CALL(MPI_Type_commit(&this->descSW[face]));
    }
    // Dimensions for edge-centered field
    #ifdef EVOLVE_VECTOR_POTENTIAL
      for(int nv = 0; nv <= AX3e ; nv++) {
        int edge; // Vector direction(=edge)

        // Map nv to a vector direction
        #if DIMENSIONS == 2
          if(nv==AX3e) {
            edge = KDIR;
          } else {
            IDEFIX_ERROR("Wrong direction for vector potential");
          }
        #elif DIMENSIONS == 3
          edge = nv;
        #else
          IDEFIX_ERROR("Cannot treat vector potential with that number of dimensions");
        #endif

        // load the array size
        for(int dir = 0; dir < 3 ; dir++) {
          size[2-dir] = grid->np_int[dir];
          start[2-dir] = data.gbeg[dir]-data.nghost[dir];
          subsize[2-dir] = data.np_int[dir];
        }

        // Extra cell in the dirs perp to field
        for(int i = 0 ; i < DIMENSIONS ; i++) {
          if(i!=edge) {
            size[2-i]++;
            subsize[2-i]++; // valid only for reading
                            //since it involves an overlap of data between procs
          }
        }
        MPI_SAFE_CALL(MPI_Type_create_subarray(3, size, subsize, start,
                                              MPI_ORDER_C, realMPI, &this->descER[nv]));
        MPI_SAFE_CALL(MPI_Type_commit(&this->descER[nv]));

        // Now for writing, it is only the last proc which keeps one additional cell,
        // so we remove what we added for reads
        for(int i = 0 ; i < DIMENSIONS ; i++) {
          if(i!=edge) {
            if(grid->xproc[i] != grid->nproc[i] - 1  ) {
              subsize[2-i]--;
            }
          }
        }
        MPI_SAFE_CALL(MPI_Type_create_subarray(3, size, subsize, start,
                                              MPI_ORDER_C, realMPI, &this->descEW[nv]));
        MPI_SAFE_CALL(MPI_Type_commit(&this->descEW[nv]));
      }
    #endif
  #endif
}

void Dump::WriteString(IdfxFileHandler fileHdl, char *str, int size) {
  #ifdef WITH_MPI
    MPI_Status status;
    MPI_SAFE_CALL(MPI_File_set_view(fileHdl, this->offset,
                                    MPI_BYTE, MPI_CHAR, "native", MPI_INFO_NULL ));
    if(idfx::prank==0) {
      MPI_SAFE_CALL(MPI_File_write(fileHdl, str, size, MPI_CHAR, &status));
    }
    offset=offset+size;
  #else
    fwrite (str, sizeof(char), size, fileHdl);
  #endif
}


void Dump::WriteSerial(IdfxFileHandler fileHdl, int ndim, int *dim,
                             DataType type, char* name, void* data ) {
  int ntot = 1;   // Number of elements to be written
  int size;

  if(type == DoubleType) size=sizeof(double);
  if(type == SingleType) size=sizeof(float);
  if(type == IntegerType) size=sizeof(int);

  // Write field name

  WriteString(fileHdl, name, NAMESIZE);

  #ifdef WITH_MPI
    MPI_Status status;
    MPI_Datatype MpiType;

    // Write data type
    MPI_SAFE_CALL(MPI_File_set_view(fileHdl, offset, MPI_BYTE,
                                    MPI_CHAR, "native", MPI_INFO_NULL ));
    if(idfx::prank==0) {
      MPI_SAFE_CALL(MPI_File_write(fileHdl, &type, 1, MPI_INT, &status));
    }
    offset=offset+sizeof(int);

    // Write dimensions
    MPI_SAFE_CALL(MPI_File_set_view(fileHdl, offset, MPI_BYTE,
                                    MPI_CHAR, "native", MPI_INFO_NULL ));
    if(idfx::prank==0) {
      MPI_SAFE_CALL(MPI_File_write(fileHdl, &ndim, 1, MPI_INT, &status));
    }
    offset=offset+sizeof(int);

    for(int n = 0 ; n < ndim ; n++) {
      MPI_SAFE_CALL(MPI_File_set_view(fileHdl, offset, MPI_BYTE,
                                      MPI_CHAR, "native", MPI_INFO_NULL ));
      if(idfx::prank==0) {
        MPI_SAFE_CALL(MPI_File_write(fileHdl, dim+n, 1, MPI_INT, &status));
      }
      offset=offset+sizeof(int);
      ntot = ntot * dim[n];
    }

    // Write raw data
    if(type == DoubleType) MpiType=MPI_DOUBLE;
    if(type == SingleType) MpiType=MPI_FLOAT;
    if(type == IntegerType) MpiType=MPI_INT;
    MPI_SAFE_CALL(MPI_File_set_view(fileHdl, offset, MPI_BYTE,
                                    MPI_CHAR, "native", MPI_INFO_NULL ));

    if(idfx::prank==0) {
      MPI_SAFE_CALL(MPI_File_write(fileHdl, data, ntot, MpiType, &status));
    }
    // increment offset accordingly
    offset += ntot*size;

  #else
    // Write type of data
    fwrite(&type, 1, sizeof(int), fileHdl);
    // Write dimensions of array
    fwrite(&ndim, 1, sizeof(int), fileHdl);
    for(int n = 0 ; n < ndim ; n++) {
      fwrite(dim+n, 1, sizeof(int), fileHdl);
      ntot = ntot * dim[n];
    }
    // Write raw data
    fwrite(data, ntot, size, fileHdl);
  #endif
}

void Dump::WriteDistributed(IdfxFileHandler fileHdl, int ndim, int *dim, int *gdim,
                                  char* name, IdfxDataDescriptor &descriptor, real* data ) {
    int64_t ntot = 1;   // Number of elements to be written

  // Define current datatype
  DataType type;
  #ifndef SINGLE_PRECISION
  type = DoubleType;
  #else
  type = SingleType;
  #endif

  // Write field name
  WriteString(fileHdl, name, NAMESIZE);

  #ifdef WITH_MPI
    MPI_Status status;
    MPI_Datatype MpiType;
    int64_t nglob = 1;

    // Write data type
    MPI_SAFE_CALL(MPI_File_set_view(fileHdl, offset, MPI_BYTE,
                                    MPI_CHAR, "native", MPI_INFO_NULL ));
    if(idfx::prank==0) {
      MPI_SAFE_CALL(MPI_File_write(fileHdl, &type, 1, MPI_INT, &status));
    }
    offset=offset+sizeof(int);

    // Write dimensions
    MPI_SAFE_CALL(MPI_File_set_view(fileHdl, offset, MPI_BYTE,
                                    MPI_CHAR, "native", MPI_INFO_NULL ));
    if(idfx::prank==0) {
      MPI_SAFE_CALL(MPI_File_write(fileHdl, &ndim, 1, MPI_INT, &status));
    }
    offset=offset+sizeof(int);

    for(int n = 0 ; n < ndim ; n++) {
      MPI_SAFE_CALL(MPI_File_set_view(fileHdl, offset, MPI_BYTE,
                                      MPI_CHAR, "native", MPI_INFO_NULL ));
      if(idfx::prank==0) {
        MPI_SAFE_CALL(MPI_File_write(fileHdl, gdim+n, 1, MPI_INT, &status));
      }
      offset=offset+sizeof(int);
      ntot = ntot * dim[n];
      nglob = nglob * gdim[n];
    }

    // Write raw data
    if(type == DoubleType) MpiType=MPI_DOUBLE;
    if(type == SingleType) MpiType=MPI_FLOAT;

    MPI_SAFE_CALL(MPI_File_set_view(fileHdl, offset, MpiType,
                                    descriptor, "native", MPI_INFO_NULL ));
    MPI_SAFE_CALL(MPI_File_write_all(fileHdl, data, ntot, MpiType, MPI_STATUS_IGNORE));

    offset=offset+nglob*sizeof(real);

  #else
    // Write type of data

    fwrite(&type, 1, sizeof(int), fileHdl);

    // Write dimensions of array
    // (in serial, dim and gdim are identical, so no need to differentiate)
    fwrite(&ndim, 1, sizeof(int), fileHdl);
    for(int n = 0 ; n < ndim ; n++) {
      fwrite(dim+n, 1, sizeof(int), fileHdl);
      ntot = ntot * dim[n];
    }

    // Write raw data
    fwrite(data, ntot, sizeof(real), fileHdl);
  #endif
}

void Dump::ReadNextFieldProperties(IdfxFileHandler fileHdl, int &ndim, int *dim,
                                         DataType &type, std::string &name) {
  char fieldName[NAMESIZE];
  #ifdef WITH_MPI
    // Read Name
    MPI_Status status;
    MPI_SAFE_CALL(MPI_File_set_view(fileHdl, this->offset, MPI_BYTE,
                                    MPI_CHAR, "native", MPI_INFO_NULL ));
    if(idfx::prank==0) {
      MPI_SAFE_CALL(MPI_File_read(fileHdl, fieldName, NAMESIZE, MPI_CHAR, &status));
    }
    offset=offset+NAMESIZE;
    // Broadcast
    MPI_SAFE_CALL(MPI_Bcast(fieldName, NAMESIZE, MPI_CHAR, 0, MPI_COMM_WORLD));
    name.assign(fieldName,strlen(fieldName));

    // Read Datatype
    MPI_SAFE_CALL(MPI_File_set_view(fileHdl, this->offset, MPI_BYTE,
                                    MPI_CHAR, "native", MPI_INFO_NULL ));
    if(idfx::prank==0) {
      MPI_SAFE_CALL(MPI_File_read(fileHdl, &type, 1, MPI_INT, &status));
    }
    offset=offset+sizeof(int);
    MPI_SAFE_CALL(MPI_Bcast(&type, 1, MPI_INT, 0, MPI_COMM_WORLD));

    // Read Dimensions
    MPI_SAFE_CALL(MPI_File_set_view(fileHdl, this->offset, MPI_BYTE,
                                    MPI_CHAR, "native", MPI_INFO_NULL ));
    if(idfx::prank==0) {
      MPI_SAFE_CALL(MPI_File_read(fileHdl, &ndim, 1, MPI_INT, &status));
    }
    offset=offset+sizeof(int);
    MPI_SAFE_CALL(MPI_Bcast(&ndim, 1, MPI_INT, 0, MPI_COMM_WORLD));

    MPI_SAFE_CALL(MPI_File_set_view(fileHdl, this->offset, MPI_BYTE,
                                    MPI_CHAR, "native", MPI_INFO_NULL ));
    if(idfx::prank==0) {
      MPI_SAFE_CALL(MPI_File_read(fileHdl, dim, ndim, MPI_INT, &status));
    }
    offset=offset+sizeof(int)*ndim;
    MPI_SAFE_CALL(MPI_Bcast(dim, ndim, MPI_INT, 0, MPI_COMM_WORLD));

  #else
    size_t numRead;

    // Read name
    numRead = fread(fieldName, sizeof(char), NAMESIZE, fileHdl);
    if(numRead<NAMESIZE) {
      IDEFIX_ERROR("Error: unexpected end of dump file");
    }
    name.assign(fieldName,strlen(fieldName));

    // Read datatype
    numRead = fread(&type, sizeof(int), 1, fileHdl);
    if(numRead<1) {
      IDEFIX_ERROR("Error: unexpected end of dump file");
    }

    // read dimensions
    numRead = fread(&ndim, sizeof(int), 1, fileHdl);
    if(numRead<1) {
      IDEFIX_ERROR("Error: unexpected end of dump file");
    }
    numRead = fread(dim, sizeof(int), ndim, fileHdl);
    if(numRead<ndim) {
      IDEFIX_ERROR("Error: unexpected end of dump file");
    }
  #endif
}

void Dump::ReadSerial(IdfxFileHandler fileHdl, int ndim, int *dim,
                            DataType type, void* data) {
  int size;
  int64_t ntot=1;
  // Get total size
  for(int i=0; i < ndim; i++) {
    ntot=ntot*dim[i];
  }
  if(type == DoubleType) size=sizeof(double);
  if(type == SingleType) size=sizeof(float);
  if(type == IntegerType) size=sizeof(int);
  if(type == BoolType) size=sizeof(bool);

  #ifdef WITH_MPI
    MPI_Status status;
    MPI_Datatype MpiType;

    if(type == DoubleType) MpiType=MPI_DOUBLE;
    if(type == SingleType) MpiType=MPI_FLOAT;
    if(type == IntegerType) MpiType=MPI_INT;
    if(type == BoolType) MpiType=MPI_CXX_BOOL;

    MPI_SAFE_CALL(MPI_File_set_view(fileHdl, this->offset, MPI_BYTE,
                                    MPI_CHAR, "native", MPI_INFO_NULL ));
    if(idfx::prank==0) {
      MPI_SAFE_CALL(MPI_File_read(fileHdl, data, ntot, MpiType, &status));
    }
    offset+= ntot*size;
    MPI_SAFE_CALL(MPI_Bcast(data, ntot, MpiType, 0, MPI_COMM_WORLD));

  #else
    size_t numRead;

    // Read raw data
    numRead = fread(data,size,ntot,fileHdl);
    if(numRead<ntot) {
      IDEFIX_ERROR("Error: unexpected end of dump file");
    }
  #endif
}

void Dump::ReadDistributed(IdfxFileHandler fileHdl, int ndim, int *dim, int *gdim,
                                 IdfxDataDescriptor &descriptor, void* data) {
  int64_t ntot=1;
  int64_t nglob=1;
  // Get total size
  for(int i=0; i < ndim; i++) {
    ntot=ntot*dim[i];
    nglob=nglob*gdim[i];
  }

  #ifdef WITH_MPI
    MPI_Datatype MpiType;

    #ifndef SINGLE_PRECISION
    MpiType = MPI_DOUBLE;
    #else
    MpiType = MPI_FLOAT;
    #endif

    MPI_SAFE_CALL(MPI_File_set_view(fileHdl, offset, MpiType,
                                    descriptor, "native", MPI_INFO_NULL ));
    MPI_SAFE_CALL(MPI_File_read_all(fileHdl, data, ntot, MpiType, MPI_STATUS_IGNORE));

    offset=offset+nglob*sizeof(real);
  #else
    size_t numRead;
    // Read raw data
    numRead = fread(data,sizeof(real),ntot,fileHdl);
    if(numRead<ntot) {
      IDEFIX_ERROR("Error: unexpected end of dump file");
    }
  #endif
}

int Dump::Read(DataBlock &data, Output& output, int readNumber ) {
  char filename[FILENAMESIZE];
  int nx[3];
  int nxglob[3];
  std::string fieldName;
  std::string eof ("eof");
  DataType type;
  int ndim;
  IdfxFileHandler fileHdl;

  idfx::pushRegion("Dump::Read");

  idfx::cout << "Dump: Reading restart file n " << readNumber << "..." << std::flush;

  // Reset timer
  timer.reset();

  // Set filename
  std::snprintf (filename, FILENAMESIZE, "dump.%04d.dmp", readNumber);

  // Make a local image of the datablock
  DataBlockHost dataHost(data);

  // open file
#ifdef WITH_MPI
  MPI_SAFE_CALL(MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY | MPI_MODE_UNIQUE_OPEN,
                              MPI_INFO_NULL, &fileHdl));
  this->offset = 0;
#else
  fileHdl = fopen(filename,"rb");
#endif
  // File is open

    // skip the header
#ifdef WITH_MPI
  this->offset += HEADERSIZE;
#else
  fseek(fileHdl, HEADERSIZE, SEEK_SET);
#endif

  // First thing is compare the total domain size
  for(int dir=0 ; dir < 3; dir++) {
    ReadNextFieldProperties(fileHdl, ndim, nx, type, fieldName);
    if(ndim>1) IDEFIX_ERROR("Wrong coordinate array dimensions while reading restart dump");
    if(nx[0] != data.mygrid->np_int[dir]) {
      idfx::cout << "dir " << dir << ", restart has " << nx[0] << " points " << std::endl;
      IDEFIX_ERROR("Domain size from the restart dump is different from the current one");
    }

    // Read coordinates
    ReadSerial(fileHdl, ndim, nx, type, scrch);

    // skip left and right edges arrays
    for (int iside=0; iside < 2; iside++) {
      ReadNextFieldProperties(fileHdl, ndim, nx, type, fieldName);
      ReadSerial(fileHdl, ndim, nx, type, scrch);
    }
    // Todo: check that coordinates are identical
  }
  // Coordinates are ok, load the bulk
  while(true) {
    ReadNextFieldProperties(fileHdl, ndim, nxglob, type, fieldName);
    //idfx::cout << "Next field is " << fieldName << " with " << ndim << " dimensions and (";
    //for(int i = 0 ; i < ndim ; i++) idfx::cout << nxglob[i] << " ";
    //idfx::cout << ") points." << std::endl;

    if(fieldName.compare(eof) == 0) {
      break;
    } else if(fieldName.compare(0,3,"Vc-") == 0) {
      // Next field is a cell-centered field

      // Matching Name is Vc-<<VcName>>
      int nv = -1;
      for(int n = 0 ; n < NVAR; n++) {
        if(fieldName.compare(3,3,data.hydro.VcName[n],0,3)==0) nv=n; // Found matching field
      }
      // Load it
      for(int dir = 0 ; dir < 3; dir++) {
        nx[dir] = dataHost.np_int[dir];
      }
      ReadDistributed(fileHdl, ndim, nx, nxglob, descC, scrch);

      if(nv<0) {
        IDEFIX_WARNING("Cannot find a field matching " + fieldName
                       + " in current running code. Skipping.");
      } else {
        // Load the scratch space in designated field
        for(int k = 0; k < nx[KDIR]; k++) {
          for(int j = 0 ; j < nx[JDIR]; j++) {
            for(int i = 0; i < nx[IDIR]; i++) {
              dataHost.Vc(nv,k+dataHost.beg[KDIR],j+dataHost.beg[JDIR],i+dataHost.beg[IDIR]) =
                                                      scrch[i + j*nx[IDIR] + k*nx[IDIR]*nx[JDIR]];
            }
          }
        }
      }
    } else if(fieldName.compare(0,3,"Vs-") == 0) {
      // Next field is a face-centered field

      // Matching Name is Vs-<<VcName>>
      #if MHD == YES
        int nv = -1;
        for(int n = 0 ; n < DIMENSIONS; n++) {
          if(fieldName.compare(3,4,data.hydro.VsName[n],0,4)==0) nv=n; // Found matching field
        }
        if(nv<0) {
          IDEFIX_ERROR("Cannot find a field matching " + fieldName
                              + " in current running code.");
        } else {
          // Load it
          for(int dir = 0 ; dir < 3; dir++) nx[dir] = dataHost.np_int[dir];
          nx[nv]++;   // Extra cell in the dir direction for cell-centered fields
          ReadDistributed(fileHdl, ndim, nx, nxglob, descSR[nv], scrch);

          for(int k = 0; k < nx[KDIR]; k++) {
            for(int j = 0 ; j < nx[JDIR]; j++) {
              for(int i = 0; i < nx[IDIR]; i++) {
                dataHost.Vs(nv,k+dataHost.beg[KDIR],j+dataHost.beg[JDIR],i+dataHost.beg[IDIR])
                            = scrch[i + j*nx[IDIR] + k*nx[IDIR]*nx[JDIR]];
              }
            }
          }
        }
      #else
        IDEFIX_WARNING("Code configured without MHD. Face-centered magnetic field components \
                        from the restart dump are skipped");
      #endif
    } else if(fieldName.compare(0,3,"Ve-") == 0) {
      #if MHD == YES && defined(EVOLVE_VECTOR_POTENTIAL)
        int nv = -1;
        for(int n = 0 ; n <= AX3e; n++) {
          if(fieldName.compare(3,4,data.hydro.VeName[n],0,4)==0) nv=n; // Found matching field
        }
        if(nv<0) {
          IDEFIX_ERROR("Cannot find a field matching " + fieldName
                              + " in current running code.");
        } else {
          int dir = 0;
          #if DIMENSIONS == 2
            if(nv==AX3e) {
              dir = KDIR;
            } else {
              IDEFIX_ERROR("Wrong direction for vector potential");
            }
          #elif DIMENSIONS == 3
            dir = nv;
          #else
            IDEFIX_ERROR("Cannot treat vector potential with that number of dimensions");
          #endif
          // Load it
          for(int i = 0 ; i < 3; i++) {
            nx[i] = dataHost.np_int[i];
          }
          // Extra cell in the dirs perp to field
          for(int i = 0 ; i < DIMENSIONS ; i++) {
            if(i!=dir) nx[i] ++;
          }

          ReadDistributed(fileHdl, ndim, nx, nxglob, descER[nv], scrch);

          for(int k = 0; k < nx[KDIR]; k++) {
            for(int j = 0 ; j < nx[JDIR]; j++) {
              for(int i = 0; i < nx[IDIR]; i++) {
                dataHost.Ve(nv,k+dataHost.beg[KDIR],j+dataHost.beg[JDIR],i+dataHost.beg[IDIR])
                            = scrch[i + j*nx[IDIR] + k*nx[IDIR]*nx[JDIR]];
              }
            }
          }
        }

      #else
      IDEFIX_WARNING("Code configured without vector potential support. Vector potentials \
                        from the restart dump are skipped");
      #endif
    } else if(fieldName.compare("time") == 0) {
      ReadSerial(fileHdl, ndim, nxglob, type, &data.t);
    } else if(fieldName.compare("dt") == 0) {
      ReadSerial(fileHdl, ndim, nxglob, type, &data.dt);
    } else if(fieldName.compare("vtkFileNumber")==0) {
      ReadSerial(fileHdl, ndim, nxglob, type, &output.vtk.vtkFileNumber);
    } else if(fieldName.compare("vtkLast")==0) {
      ReadSerial(fileHdl, ndim, nxglob, type, &output.vtkLast);
    } else if(fieldName.compare("dumpFileNumber")==0) {
      ReadSerial(fileHdl, ndim, nxglob, type, &this->dumpFileNumber);
    } else if(fieldName.compare("dumpLast")==0) {
      ReadSerial(fileHdl, ndim, nxglob, type, &output.dumpLast);
    } else if(fieldName.compare("analysisLast")==0) {
      ReadSerial(fileHdl, ndim, nxglob, type, &output.analysisLast);
    } else if(fieldName.compare("geometry")==0) {
      ReadSerial(fileHdl, ndim, nxglob, type, &this->geometry);
    } else if(fieldName.compare("periodicity")==0) {
      ReadSerial(fileHdl, ndim, nxglob, type, &this->periodicity);
    } else {
      ReadSerial(fileHdl,ndim, nxglob, type, scrch);
      IDEFIX_WARNING("Unknown field "+fieldName+" in restart dump. Skipping.");
    }
  }

  #ifdef WITH_MPI
  MPI_SAFE_CALL(MPI_File_close(&fileHdl));
  #else
  fclose(fileHdl);
  #endif

  // Send to device
  dataHost.SyncToDevice();

  idfx::cout << "done in " << timer.seconds() << " s." << std::endl;
  idfx::cout << "Restarting from t=" << data.t << "." << std::endl;

  idfx::popRegion();

  return(0);
}


int Dump::Write(DataBlock &data, Output& output) {
  char filename[FILENAMESIZE];
  char fieldName[NAMESIZE+1]; // +1 is just in case
  int nx[3];
  int nxtot[3];

  #ifndef SINGLE_PRECISION
  const DataType realType = DoubleType;
  #else
  const DataType realType = SingleType;
  #endif
  IdfxFileHandler fileHdl;

  idfx::pushRegion("Dump::Write");

  idfx::cout << "Dump: Write file n " << dumpFileNumber << "..." << std::flush;

  // Reset timer
  timer.reset();

  // Set filename
  std::snprintf(filename, FILENAMESIZE, "dump.%04d.dmp", dumpFileNumber);
  dumpFileNumber++;   // For next one

  // open file
#ifdef WITH_MPI
// Open file for creating, return error if file already exists.
  int err = MPI_File_open(MPI_COMM_WORLD, filename,
                              MPI_MODE_CREATE | MPI_MODE_RDWR
                              | MPI_MODE_EXCL | MPI_MODE_UNIQUE_OPEN,
                              MPI_INFO_NULL, &fileHdl);
  if (err != MPI_SUCCESS)  {
    // File exists, delete it before reopening
    if(idfx::prank == 0) {
      MPI_File_delete(filename,MPI_INFO_NULL);
    }
    MPI_SAFE_CALL(MPI_File_open(MPI_COMM_WORLD, filename,
                              MPI_MODE_CREATE | MPI_MODE_RDWR
                              | MPI_MODE_EXCL | MPI_MODE_UNIQUE_OPEN,
                              MPI_INFO_NULL, &fileHdl));
  }

  this->offset = 0;
#else
  fileHdl = fopen(filename,"wb");
#endif
  // File is open
  // First thing we need are coordinates: init a host mirror and sync it
  GridHost gridHost(*data.mygrid);
  gridHost.SyncFromDevice();

  // Test endianness
  std::string endian;
  int tmp1 = 1;
  unsigned char *tmp2 = (unsigned char *) &tmp1;
  if (*tmp2 != 0) {
    endian = "little";
  } else {
    endian = "big";
  }

  char header[HEADERSIZE];
  std::snprintf(header, HEADERSIZE, "Idefix %s Dump Data %s endian", GITVERSION, endian.c_str());
  WriteString(fileHdl, header, HEADERSIZE);

  for(int dir = 0; dir < 3 ; dir++) {
    // cell centers
    std::snprintf(fieldName, NAMESIZE, "x%d",dir+1);
    WriteSerial(fileHdl, 1, &gridHost.np_int[dir], realType, fieldName,
                reinterpret_cast<void*> (gridHost.x[dir].data()+gridHost.nghost[dir]));
    // cell left edges
    std::snprintf(fieldName, NAMESIZE, "xl%d",dir+1);
    WriteSerial(fileHdl, 1, &gridHost.np_int[dir], realType, fieldName,
                reinterpret_cast<void*> (gridHost.xl[dir].data()+gridHost.nghost[dir]));
    // cell right edges
    std::snprintf(fieldName, NAMESIZE, "xr%d",dir+1);
    WriteSerial(fileHdl, 1, &gridHost.np_int[dir], realType, fieldName,
                reinterpret_cast<void*> (gridHost.xr[dir].data()+gridHost.nghost[dir]));
  }

  // Then write raw data from Vc
  DataBlockHost dataHost(data);
  dataHost.SyncFromDevice();

  for(int nv = 0 ; nv <NVAR ; nv++) {
    std::snprintf(fieldName,NAMESIZE,"Vc-%s",data.hydro.VcName[nv].c_str());
    // Load the active domain in the scratch space
    for(int i = 0; i < 3 ; i++) {
      nx[i] = dataHost.np_int[i];
      nxtot[i] = gridHost.np_int[i];
    }

    for(int k = 0; k < nx[KDIR]; k++) {
      for(int j = 0 ; j < nx[JDIR]; j++) {
        for(int i = 0; i < nx[IDIR]; i++) {
          scrch[i + j*nx[IDIR] + k*nx[IDIR]*nx[JDIR]] = dataHost.Vc(nv,k+dataHost.beg[KDIR],
                                                                       j+dataHost.beg[JDIR],
                                                                       i+dataHost.beg[IDIR]);
        }
      }
    }
    WriteDistributed(fileHdl, 3, nx, nxtot, fieldName, this->descC, scrch);
  }

  #if MHD == YES
    // write staggered field components
    for(int nv = 0 ; nv <DIMENSIONS ; nv++) {
      std::snprintf(fieldName,NAMESIZE,"Vs-%s",data.hydro.VsName[nv].c_str());
      // Load the active domain in the scratch space
      for(int i = 0; i < 3 ; i++) {
        nx[i] = dataHost.np_int[i];
        nxtot[i] = gridHost.np_int[i];
      }
      // If it is the last datablock of the dimension, increase the size by one to get the last
      //active face of the staggered mesh.
      if(data.mygrid->xproc[nv] == data.mygrid->nproc[nv] - 1  ) nx[nv]++;
      nxtot[nv]++;

      for(int k = 0; k < nx[KDIR]; k++) {
        for(int j = 0 ; j < nx[JDIR]; j++) {
          for(int i = 0; i < nx[IDIR]; i++) {
            scrch[i + j*nx[IDIR] + k*nx[IDIR]*nx[JDIR] ] = dataHost.Vs(nv,k+dataHost.beg[KDIR],
                                                                          j+dataHost.beg[JDIR],
                                                                          i+dataHost.beg[IDIR]);
          }
        }
      }
      WriteDistributed(fileHdl, 3, nx, nxtot, fieldName, this->descSW[nv], scrch);
    }
    #ifdef EVOLVE_VECTOR_POTENTIAL
      // write edge field components
      for(int nv = 0 ; nv <= AX3e ; nv++) {
        std::snprintf(fieldName,NAMESIZE,"Ve-%s",data.hydro.VeName[nv].c_str());
        int edge = 0;
        #if DIMENSIONS == 2
          if(nv==AX3e) {
            edge = KDIR;
          } else {
            IDEFIX_ERROR("Wrong direction for vector potential");
          }
        #elif DIMENSIONS == 3
          edge = nv;
        #else
          IDEFIX_ERROR("Cannot treat vector potential with that number of dimensions");
        #endif
        // Load the active domain in the scratch space
        for(int i = 0; i < 3 ; i++) {
          nx[i] = dataHost.np_int[i];
          nxtot[i] = gridHost.np_int[i];
        }
        // If it is the last datablock of the dimension, increase the size by one in the direction
        // perpendicular to the vector.

        for(int i = 0 ; i < DIMENSIONS ; i++) {
          if(i != edge) {
            if(data.mygrid->xproc[i] == data.mygrid->nproc[i] - 1) nx[i]++;
            nxtot[i]++;
          }
        }
        for(int k = 0; k < nx[KDIR]; k++) {
          for(int j = 0 ; j < nx[JDIR]; j++) {
            for(int i = 0; i < nx[IDIR]; i++) {
              scrch[i + j*nx[IDIR] + k*nx[IDIR]*nx[JDIR] ] = dataHost.Ve(nv,k+dataHost.beg[KDIR],
                                                                            j+dataHost.beg[JDIR],
                                                                            i+dataHost.beg[IDIR]);
            }
          }
        }
        WriteDistributed(fileHdl, 3, nx, nxtot, fieldName, this->descEW[nv], scrch);
      }
    #endif
  #endif

  // Write some raw data
  nx[0] = 1;
  std::snprintf(fieldName,NAMESIZE, "time");
  WriteSerial(fileHdl, 1, nx, realType, fieldName, &data.t);
  std::snprintf(fieldName,NAMESIZE, "dt");
  WriteSerial(fileHdl, 1, nx, realType, fieldName, &data.dt);
  std::snprintf(fieldName,NAMESIZE, "vtkFileNumber");
  WriteSerial(fileHdl, 1, nx, IntegerType, fieldName, &output.vtk.vtkFileNumber);
  std::snprintf(fieldName,NAMESIZE, "vtkLast");
  WriteSerial(fileHdl, 1, nx, realType, fieldName, &output.vtkLast);
  std::snprintf(fieldName,NAMESIZE, "dumpFileNumber");
  WriteSerial(fileHdl, 1, nx, IntegerType, fieldName, &this->dumpFileNumber);
  std::snprintf(fieldName,NAMESIZE, "dumpLast");
  WriteSerial(fileHdl, 1, nx, realType, fieldName, &output.dumpLast);
  std::snprintf(fieldName,NAMESIZE, "analysisLast");
  WriteSerial(fileHdl, 1, nx, realType, fieldName, &output.analysisLast);
  std::snprintf(fieldName,NAMESIZE, "geometry");
  WriteSerial(fileHdl, 1, nx, IntegerType, fieldName, &this->geometry);

  nx[0] = 3;
  std::snprintf(fieldName,NAMESIZE, "periodicity");
  WriteSerial(fileHdl, 1, nx, IntegerType, fieldName, &this->periodicity);

  // Write end of file
  scrch[0] = 0.0;
  std::snprintf(fieldName,NAMESIZE,"eof");
  nx[0] = 1;
  WriteSerial(fileHdl, 1, nx, realType, fieldName, scrch);

#ifdef WITH_MPI
  MPI_SAFE_CALL(MPI_File_close(&fileHdl));
#else
  fclose(fileHdl);
#endif


  idfx::cout << "done in " << timer.seconds() << " s." << std::endl;
  idfx::popRegion();
  // One day, we will have a return code.

  return(0);
}
