// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2021 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include <string>
#include <sstream>
#include <iomanip>
#include "outputVtk.hpp"
#include "gitversion.hpp"
#include "../idefix.hpp"

#define VTK_RECTILINEAR_GRID    14
#define VTK_STRUCTURED_GRID     35

#ifndef VTK_FORMAT
  #if GEOMETRY == CARTESIAN || GEOMETRY == CYLINDRICAL
    #define VTK_FORMAT  VTK_RECTILINEAR_GRID
  #else
    #define VTK_FORMAT  VTK_STRUCTURED_GRID
  #endif
#endif

// Whether or not we wand the staggered field in vtk outputs
#define WRITE_STAGGERED_FIELD

// Wheter or not we want the electrical current in vtk outputs
// (this is only defined if hydro requires it)
#define WRITE_CURRENT

// Whether of not we want to write EMFs
#define WRITE_EMF

// Whether or not we want to write the cfl dt
#define WRITE_DT
/* ---------------------------------------------------------
    The following macros are specific to this file only
    and are used to ease up serial/parallel implementation
    for writing strings and real arrays
   --------------------------------------------------------- */



void OutputVTK::WriteHeaderString(const char* header, IdfxFileHandler fvtk) {
#ifdef WITH_MPI
  MPI_Status status;
  MPI_SAFE_CALL(MPI_File_set_view(fvtk, this->offset, MPI_BYTE,
                                  MPI_CHAR, "native", MPI_INFO_NULL ));
  if(idfx::prank==0) {
    MPI_SAFE_CALL(MPI_File_write(fvtk, header, strlen(header), MPI_CHAR, &status));
  }
  offset=offset+strlen(header);
#else
  fprintf (fvtk, "%s", header);
#endif
}

void OutputVTK::WriteHeaderFloat(float* buffer, int64_t nelem, IdfxFileHandler fvtk) {
#ifdef WITH_MPI
  MPI_Status status;
  MPI_SAFE_CALL(MPI_File_set_view(fvtk, this->offset, MPI_BYTE, MPI_CHAR,
                                  "native", MPI_INFO_NULL ));
  if(idfx::prank==0) {
    MPI_SAFE_CALL(MPI_File_write(fvtk, buffer, nelem, MPI_FLOAT, &status));
  }
  offset=offset+nelem*sizeof(float);
#else
  fwrite(buffer, sizeof(float), nelem, fvtk);
#endif
}

/* Main constructor */
OutputVTK::OutputVTK(Input &input, DataBlock &datain) {
  // Init the output period
  this->tperiod=input.GetReal("Output","vtk",0);
  this->tnext = datain.t;

  // Initialize the output structure
  // Create a local datablock as an image of gridin
  DataBlockHost data(datain);
  data.SyncFromDevice();

  // Pointer to global grid
  GridHost grid(*datain.mygrid);
  grid.SyncFromDevice();

  /* Note that there are two kinds of dimensions:
     - nx1, nx2, nx3, derived from the grid, which are the global dimensions
     - nx1loc,nx2loc,n3loc, which are the local dimensions of the current datablock
  */

  // Create the coordinate array required in VTK files
  this->nx1 = grid.np_int[IDIR];
  this->nx2 = grid.np_int[JDIR];
  this->nx3 = grid.np_int[KDIR];

  this->nx1loc = data.np_int[IDIR];
  this->nx2loc = data.np_int[JDIR];
  this->nx3loc = data.np_int[KDIR];

  // Vector array where we store the pencil before write
  this->Vwrite = new float[nx1loc+IOFFSET];

  // Temporary storage on host for 3D arrays
  this->vect3D = new float[nx1loc*nx2loc*nx3loc];

  // Essentially does nothing
  this->vtkFileNumber = 0;

  // Test endianness
  int tmp1 = 1;
  this->shouldSwapEndian = 0;
  unsigned char *tmp2 = (unsigned char *) &tmp1;
  if (*tmp2 != 0)
    this->shouldSwapEndian = 1;

  // Store coordinates for later use
  this->xnode = new float[nx1+IOFFSET];
  this->ynode = new float[nx2+JOFFSET];
  this->znode = new float[nx3+KOFFSET];

  for (int32_t i = 0; i < nx1 + IOFFSET; i++) {
    xnode[i] = BigEndian(grid.xl[IDIR](i + grid.nghost[IDIR]));
  }
  for (int32_t j = 0; j < nx2 + JOFFSET; j++)    {
    ynode[j] = BigEndian(grid.xl[JDIR](j + grid.nghost[JDIR]));
  }
  for (int32_t k = 0; k < nx3 + KOFFSET; k++) {
    if(DIMENSIONS==2)
      znode[k] = BigEndian(0.0);
    else
      znode[k] = BigEndian(grid.xl[KDIR](k + grid.nghost[KDIR]));
  }
#if VTK_FORMAT == VTK_STRUCTURED_GRID   // VTK_FORMAT
  /* -- Allocate memory for node_coord which is later used -- */
  node_coord = new float[(nx1+IOFFSET)*3];
#endif

  // Creat MPI view when using MPI I/O
#ifdef WITH_MPI
  int start[3];
  int size[3];
  int subsize[3];

  for(int dir = 0; dir < 3 ; dir++) {
    // VTK assumes Fortran array ordering, hence arrays dimensions are filled backwards
    start[2-dir] = datain.gbeg[dir]-grid.nghost[dir];
    size[2-dir] = grid.np_int[dir];
    subsize[2-dir] = datain.np_int[dir];
  }

  MPI_SAFE_CALL(MPI_Type_create_subarray(3, size, subsize, start, MPI_ORDER_C,
                                         MPI_FLOAT, &this->view));
  MPI_SAFE_CALL(MPI_Type_commit(&this->view));
#endif
}

int OutputVTK::CheckForWrite(DataBlock &datain) {
  // Do we need an output?
  if(datain.t < this->tnext) return(0);

  this->tnext+= this->tperiod;
  return(this->Write(datain));
}

int OutputVTK::Write(DataBlock &datain) {
  idfx::pushRegion("OutputVTK::Write");

  IdfxFileHandler fileHdl;
  std::string filename;

  idfx::cout << "OutputVTK::Write file n " << vtkFileNumber << "..." << std::flush;

  timer.reset();

  // Create a copy of the dataBlock on Host, and sync it.
  DataBlockHost data(datain);
  data.SyncFromDevice();

  std::stringstream ssfileName, ssvtkFileNum;
  ssvtkFileNum << std::setfill('0') << std::setw(4) << vtkFileNumber;
  ssfileName << "data." << ssvtkFileNum.str() << ".vtk";
  filename = ssfileName.str();

  // Open file and write header
#ifdef WITH_MPI
  MPI_SAFE_CALL(MPI_File_open(MPI_COMM_WORLD, filename.c_str(),
                              MPI_MODE_CREATE | MPI_MODE_RDWR | MPI_MODE_UNIQUE_OPEN,
                              MPI_INFO_NULL, &fileHdl));
  this->offset = 0;
#else
  fileHdl = fopen(filename.c_str(),"wb");
#endif

  WriteHeader(fileHdl);

  // Write field one by one
  for(int nv = 0 ; nv < NVAR ; nv++) {
    for(int k = data.beg[KDIR]; k < data.end[KDIR] ; k++ ) {
      for(int j = data.beg[JDIR]; j < data.end[JDIR] ; j++ ) {
        for(int i = data.beg[IDIR]; i < data.end[IDIR] ; i++ ) {
          vect3D[i-data.beg[IDIR] + (j-data.beg[JDIR])*nx1loc + (k-data.beg[KDIR])*nx1loc*nx2loc]
              = BigEndian(static_cast<float>(data.Vc(nv,k,j,i)));
        }
      }
    }
    WriteScalar(fileHdl, vect3D, datain.hydro.VcName[nv]);
  }

#if MHD == YES
  #ifdef WRITE_STAGGERED_FIELD
  for(int nv = 0 ; nv < DIMENSIONS ; nv++) {
    for(int k = data.beg[KDIR]; k < data.end[KDIR] ; k++ ) {
      for(int j = data.beg[JDIR]; j < data.end[JDIR] ; j++ ) {
        for(int i = data.beg[IDIR]; i < data.end[IDIR] ; i++ ) {
          vect3D[i-data.beg[IDIR] + (j-data.beg[JDIR])*nx1loc + (k-data.beg[KDIR])*nx1loc*nx2loc]
              = BigEndian(static_cast<float>(data.Vs(nv,k,j,i)));
        }
      }
    }
    WriteScalar(fileHdl, vect3D, datain.hydro.VsName[nv]);
  }
  #endif // WRITE_STAGGERED_FIELD

  #ifdef WRITE_CURRENT
  if(data.haveCurrent) {
    for(int nv = 0 ; nv < 3 ; nv++) {
      for(int k = data.beg[KDIR]; k < data.end[KDIR] ; k++ ) {
        for(int j = data.beg[JDIR]; j < data.end[JDIR] ; j++ ) {
          for(int i = data.beg[IDIR]; i < data.end[IDIR] ; i++ ) {
            vect3D[i-data.beg[IDIR] + (j-data.beg[JDIR])*nx1loc + (k-data.beg[KDIR])*nx1loc*nx2loc]
                = BigEndian(static_cast<float>(data.J(nv,k,j,i)));
          }
        }
      }
      std::string varname = "JX"+std::to_string(nv+1);
      WriteScalar(fileHdl, vect3D, varname);
    }
  }
  #endif // WRITE_CURRENT

  #ifdef WRITE_EMF
  {
    for(int k = data.beg[KDIR]; k < data.end[KDIR] ; k++ ) {
      for(int j = data.beg[JDIR]; j < data.end[JDIR] ; j++ ) {
        for(int i = data.beg[IDIR]; i < data.end[IDIR] ; i++ ) {
          vect3D[i-data.beg[IDIR] + (j-data.beg[JDIR])*nx1loc + (k-data.beg[KDIR])*nx1loc*nx2loc]
              = BigEndian(static_cast<float>(data.Ex3(k,j,i)));
        }
      }
    }
    std::string varname = "EX3";
    WriteScalar(fileHdl, vect3D, varname);
    #if DIMENSIONS == 3
    for(int k = data.beg[KDIR]; k < data.end[KDIR] ; k++ ) {
      for(int j = data.beg[JDIR]; j < data.end[JDIR] ; j++ ) {
        for(int i = data.beg[IDIR]; i < data.end[IDIR] ; i++ ) {
          vect3D[i-data.beg[IDIR] + (j-data.beg[JDIR])*nx1loc + (k-data.beg[KDIR])*nx1loc*nx2loc]
              = BigEndian(static_cast<float>(data.Ex1(k,j,i)));
        }
      }
    }
    varname = "EX1";
    WriteScalar(fileHdl, vect3D, varname);

    for(int k = data.beg[KDIR]; k < data.end[KDIR] ; k++ ) {
      for(int j = data.beg[JDIR]; j < data.end[JDIR] ; j++ ) {
        for(int i = data.beg[IDIR]; i < data.end[IDIR] ; i++ ) {
          vect3D[i-data.beg[IDIR] + (j-data.beg[JDIR])*nx1loc + (k-data.beg[KDIR])*nx1loc*nx2loc]
              = BigEndian(static_cast<float>(data.Ex2(k,j,i)));
        }
      }
    }
    varname = "EX2";
    WriteScalar(fileHdl, vect3D, varname);
    #endif // DIMENSIONS
  }
  #endif // WRITE_EMF
#endif// MHD

#ifdef WRITE_DT
  {
    for(int k = data.beg[KDIR]; k < data.end[KDIR] ; k++ ) {
      for(int j = data.beg[JDIR]; j < data.end[JDIR] ; j++ ) {
        for(int i = data.beg[IDIR]; i < data.end[IDIR] ; i++ ) {
          vect3D[i-data.beg[IDIR] + (j-data.beg[JDIR])*nx1loc + (k-data.beg[KDIR])*nx1loc*nx2loc]
              = BigEndian(static_cast<float>(1/data.InvDt(k,j,i)));
        }
      }
    }
    std::string varname = "Dt";
    WriteScalar(fileHdl, vect3D, varname);
  }
#endif

#ifdef WITH_MPI
  MPI_SAFE_CALL(MPI_File_close(&fileHdl));
#else
  fclose(fileHdl);
#endif

  vtkFileNumber++;
  // Make file number
  idfx::cout << "done in " << timer.seconds() << " s." << std::endl;

  idfx::popRegion();
  // One day, we will have a return code.
  return(0);
}


/* ********************************************************************* */
void OutputVTK::WriteHeader(IdfxFileHandler fvtk) {
/*!
* Write VTK header in parallel or serial mode.
*
* \param [in]  fvtk  pointer to file
* \param [in]  grid  pointer to an array of Grid structures
*
* \todo  Write the grid using several processors.
*********************************************************************** */

  std::string header;
  std::stringstream ssheader;
  float x1, x2, x3;

  /* -------------------------------------------
  1. File version and identifier
  ------------------------------------------- */

  ssheader << "# vtk DataFile Version 2.0" << std::endl;

  /* -------------------------------------------
  2. Header
  ------------------------------------------- */

  ssheader << "Idefix " << GITVERSION << " VTK Data" << std::endl;

  /* ------------------------------------------
  3. File format
  ------------------------------------------ */

  ssheader << "BINARY" << std::endl;

  /* ------------------------------------------
  4. Dataset structure
  ------------------------------------------ */

#if VTK_FORMAT == VTK_RECTILINEAR_GRID
  ssheader << "DATASET RECTILINEAR_GRID" << std::endl;
#elif VTK_FORMAT == VTK_STRUCTURED_GRID
  ssheader << "DATASET STRUCTURED_GRID" << std::endl;
#endif


  ssheader << "DIMENSIONS " << nx1 + IOFFSET << " " << nx2 + JOFFSET << " " << nx3 + KOFFSET
           << std::endl;

  header = ssheader.str();
  WriteHeaderString(header.c_str(), fvtk);

#if VTK_FORMAT == VTK_RECTILINEAR_GRID

  /* -- Write rectilinear grid information -- */
  std::stringstream coordx, coordy, coordz;

  coordx << "X_COORDINATES " << nx1 + IOFFSET << " float" << std::endl;
  header = coordx.str();
  WriteHeaderString(header.c_str(), fvtk);

  WriteHeaderFloat(xnode, nx1 + IOFFSET, fvtk);

  coordy << std::endl << "Y_COORDINATES " << nx2 + JOFFSET << " float" << std::endl;
  header = coordy.str();
  WriteHeaderString(header.c_str(), fvtk);

  WriteHeaderFloat(ynode, nx2 + JOFFSET, fvtk);

  coordz << std::endl << "Z_COORDINATES " << nx3 + KOFFSET << " float" << std::endl;
  header = coordz.str();
  WriteHeaderString(header.c_str(), fvtk);

  WriteHeaderFloat(znode, nx3 + KOFFSET, fvtk);

#elif VTK_FORMAT == VTK_STRUCTURED_GRID

  /* -- define node_coord -- */
  std::stringstream ssheader2;

  ssheader2 << "POINTS " << (nx1 + IOFFSET) * (nx2 + JOFFSET) * (nx3 + KOFFSET) << " float"
            << std::endl;
  header = ssheader2.str();
  WriteHeaderString(header.c_str(), fvtk);

  /* -- Write structured grid information -- */

  x1 = x2 = x3 = 0.0;
  for (int32_t k = 0; k < nx3 + KOFFSET; k++) {
    for (int32_t j = 0; j < nx2 + JOFFSET; j++) {
      for (int32_t i = 0; i < nx1 + IOFFSET; i++) {
        // BigEndian allows us to get back to little endian when needed
        D_EXPAND( x1 = BigEndian(xnode[i]);  ,
                  x2 = BigEndian(ynode[j]);  ,
                  x3 = BigEndian(znode[k]);  )

  #if (GEOMETRY == CARTESIAN) || (GEOMETRY == CYLINDRICAL)
        node_coord[3*i+IDIR] = BigEndian(x1);
        node_coord[3*i+JDIR] = BigEndian(x2);
        node_coord[3*i+KDIR] = BigEndian(x3);

  #elif GEOMETRY == POLAR
        node_coord[3*i+IDIR] = BigEndian(x1 * cos(x2));
        node_coord[3*i+JDIR] = BigEndian(x1 * sin(x2));
        node_coord[3*i+KDIR] = BigEndian(x3);

  #elif GEOMETRY == SPHERICAL
    #if DIMENSIONS == 2
        node_coord[3*i+IDIR] = BigEndian(x1 * sin(x2));
        node_coord[3*i+JDIR] = BigEndian(x1 * cos(x2));
        node_coord[3*i+KDIR] = BigEndian(0.0);

    #elif DIMENSIONS == 3
        node_coord[3*i+IDIR] = BigEndian(x1 * sin(x2) * cos(x3));
        node_coord[3*i+JDIR] = BigEndian(x1 * sin(x2) * sin(x3));
        node_coord[3*i+KDIR] = BigEndian(x1 * cos(x2));
    #endif // DIMENSIONS
  #endif // GEOMETRY
      }

      WriteHeaderFloat(node_coord, 3 * (nx1 + IOFFSET),fvtk);
    }
  }

#endif // VTK_FORMAT

  /* -----------------------------------------------------
  5. Dataset attributes [will continue by later calls
    to WriteVTK_Vector() or WriteVTK_Scalar()...]
  ----------------------------------------------------- */

  std::stringstream ssheader3;
  ssheader3 << std::endl << "CELL_DATA " << nx1 * nx2 * nx3 << std::endl;
  header = ssheader3.str();
  WriteHeaderString(header.c_str(), fvtk);
}
#undef VTK_STRUCTURED_GRID
#undef VTK_RECTILINEAR_GRID


/* ********************************************************************* */
void OutputVTK::WriteScalar(IdfxFileHandler fvtk, float* Vin,  std::string &var_name) {
/*!
* Write VTK scalar field.
*
* \param [in]   fvtk       pointer to file (handle)
* \param [in]   V          pointer to 3D data array
* \param [in] unit     the corresponding cgs unit (if specified, 1 otherwise)
* \param [in]   var_name   the variable name appearing in the VTK file
* \param [in]   grid       pointer to an array of Grid structures
*********************************************************************** */

  int i, j, k;
  std::stringstream ssheader;

  ssheader << std::endl << "SCALARS " << var_name.c_str() << " float" << std::endl;
  ssheader << "LOOKUP_TABLE default" << std::endl;
  std::string header(ssheader.str());

  WriteHeaderString(header.c_str(), fvtk);

#ifdef WITH_MPI
  MPI_SAFE_CALL(MPI_File_set_view(fvtk, this->offset, MPI_FLOAT, this->view,
                                  "native", MPI_INFO_NULL));
  MPI_SAFE_CALL(MPI_File_write_all(fvtk, Vin, nx1loc*nx2loc*nx3loc, MPI_FLOAT, MPI_STATUS_IGNORE));
  this->offset = this->offset + sizeof(float)*nx1*nx2*nx3;
#else
  fwrite(Vin,sizeof(float),nx1loc*nx2loc*nx3loc,fvtk);
#endif
}

/* ****************************************************************************/
/** Determines if the machine is little-endian.  If so,
  it will force the data to be big-endian.
@param in_number floating point number to be converted in big endian */
/* *************************************************************************** */

float OutputVTK::BigEndian(float in_number) {
  if (shouldSwapEndian) {
    unsigned char *bytes = (unsigned char*) &in_number;
    unsigned char tmp = bytes[0];
    bytes[0] = bytes[3];
    bytes[3] = tmp;
    tmp = bytes[1];
    bytes[1] = bytes[2];
    bytes[2] = tmp;
  }
  return(in_number);
}
