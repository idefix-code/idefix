// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include "dumpImage.hpp"
#include "dataBlock.hpp"
#include "idefix.hpp"


#define  HEADERSIZE 128

DumpImage::DumpImage(std::string filename, DataBlock *data, bool enableDomainDecomposition) {
  idfx::pushRegion("DumpImage::DumpImage");

  int nx[3];
  std::string fieldName;
  std::string eof ("eof");
  DataType type;
  int ndim;
  IdfxFileHandler fileHdl;
  Dump dump(data);

  idfx::cout << "DumpImage: loading restart file " << filename << "..." << std::flush;

  // open file
#ifdef WITH_MPI
  MPI_SAFE_CALL(MPI_File_open(MPI_COMM_WORLD, filename.c_str(),
                              MPI_MODE_RDONLY | MPI_MODE_UNIQUE_OPEN,
                              MPI_INFO_NULL, &fileHdl));
  dump.offset = 0;
#else
  fileHdl = fopen(filename.c_str(),"rb");
  if(fileHdl == NULL) {
    std::stringstream msg;
    msg << "Failed to open dump file: " << filename << std::endl;
    IDEFIX_ERROR(msg);
  }
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
    this->np_glob[dir] = nx[0];
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

  if(enableDomainDecomposition) {
    #ifdef WITH_MPI
      GridBox gridBox = GetBox(data);
      // Create sub-x domains
      for(int dir = 0 ; dir < 3 ; dir ++) {
        IdefixHostArray1D<real> xLoc("DumpImageX",gridBox.size[dir]);
        IdefixHostArray1D<real> xlLoc("DumpImageXl",gridBox.size[dir]);
        IdefixHostArray1D<real> xrLoc("DumpImageXr",gridBox.size[dir]);
        for(int i = 0 ; i < gridBox.size[dir] ; i++) {
        xLoc(i) = x[dir](i+gridBox.start[dir]);
        xlLoc(i) = xl[dir](i+gridBox.start[dir]);
        xrLoc(i) = xr[dir](i+gridBox.start[dir]);
        }
        // Replace former arrays
        x[dir] = xLoc;
        xl[dir] = xlLoc;
        xr[dir] = xrLoc;
        np_int[dir] = gridBox.size[dir];
      }
      // Create dedicated reading datatype in the dump
      dump.CreateMPIDataType(gridBox, true);
    #else
      IDEFIX_WARNING("Can't create a DumpImage with domain decomposition without MPI enabled");
      enableDomainDecomposition = false;
    #endif
  }

  // Read the other fields
  while(true) {
    dump.ReadNextFieldProperties(fileHdl, ndim, nx, type, fieldName);
    if(fieldName.compare(eof) == 0) {
      break;
    } else if( ndim == 3) {
      // Load 3D field (raw data)
      // Make a new view of the right dimension for the raw data

      if(enableDomainDecomposition) {
        int nxloc[3];
        int nType = 0;  // = 0 for cell-centered; = 1 for face-centered, =2 for edge-centered
        int edgeDir = -1;
        int surfaceDir = -1;
        for(int dir = 0 ; dir < 3; dir++) {
          nxloc[dir] = this->np_int[dir];
          // Try to guess whether it's an edge or surface array depending on the extension
          // Surface: only one direction has +1 point, and it is the direction of the field
          // Edge: 2 directions with +1 point, the direction is that without +1
          if(nx[dir]==np_glob[dir]+1) {
            if(nType==0) {
              // This is either a surface or edge field
              // probably a surface
              surfaceDir = dir;
              // or an edge
              if(dir==JDIR) {
                edgeDir = IDIR;
              }
            } else {
              // surely an edge
              if(dir == JDIR) {
                edgeDir = KDIR;
              } else if(dir==KDIR && surfaceDir == -1) {
                edgeDir = JDIR;
              }
            }
            nType++;
            nxloc[dir]++;
          }
        }
        // Allocate an array of the right size
        this->arrays[fieldName] = IdefixHostArray3D<real>
                                    ("DumpImage"+fieldName,nxloc[2],nxloc[1],nxloc[0] );

        // load the data
        if(nType==0) {
          dump.ReadDistributed(fileHdl, ndim, nxloc, nx, dump.descCR,
                              reinterpret_cast<void*>(this->arrays[fieldName].data()) );
        } else if(nType==1) {
          // Edge type
          dump.ReadDistributed(fileHdl, ndim, nxloc, nx, dump.descSR[surfaceDir],
                              reinterpret_cast<void*>(this->arrays[fieldName].data()) );
        } else if(nType==2) {
          dump.ReadDistributed(fileHdl, ndim, nxloc, nx, dump.descER[edgeDir],
                              reinterpret_cast<void*>(this->arrays[fieldName].data()) );
        }
      } else {
        this->arrays[fieldName] = IdefixHostArray3D<real>("DumpImage"+fieldName,nx[2],nx[1],nx[0]);
        // Load it
        dump.ReadSerial(fileHdl,ndim,nx,type,
                        reinterpret_cast<void*>(this->arrays[fieldName].data()));
      }
    } else if(fieldName.compare("time") == 0) {
      dump.ReadSerial(fileHdl, ndim, nx, type, &this->time);
    } else if(fieldName.compare("geometry")==0) {
      dump.ReadSerial(fileHdl, ndim, nx, type, &this->geometry);
    } else if(fieldName.compare("centralMass")==0) {
      dump.ReadSerial(fileHdl, ndim, nx, type, &this->centralMass);
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

GridBox DumpImage::GetBox(DataBlock *data) {
  GridBox gridBox;
  for(int dir = 0 ; dir < 3 ; dir ++) {
    std::array<int,2> boundIndex;
    auto xl = this->xl[dir];
    auto xr = this->xr[dir];
    real xbeg = data->xbeg[dir];
    real xend = data->xend[dir];

    // left bound
    int i = 0;
    while(xl(i) <= xbeg  && i < xl.extent(0)) {
      boundIndex[0] = i;
      i++;
    }
    // right bound
    i = 0;
    while(xr(i) <= xend  && i < xr.extent(0)) {
      boundIndex[1] = i;
      i++;
    }

    // Bounds are stored in boundIndex
    // Now we initialize the gridBox in the desired dimension
    gridBox.start[dir] = boundIndex[0];
    gridBox.size[dir] = boundIndex[1] - boundIndex[0]+1;
    gridBox.sizeGlob[dir] = this->x[dir].extent(0);
  }
  return gridBox;
}
