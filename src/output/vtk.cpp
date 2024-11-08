// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include "vtk.hpp"
#include <limits.h>
#include <string>
#include <cstdio>
#include <sstream>
#include <iomanip>
#if __has_include(<filesystem>)
  #include <filesystem> // NOLINT [build/c++17]
  namespace fs = std::filesystem;
#elif __has_include(<experimental/filesystem>)
  #include <experimental/filesystem>
  namespace fs = std::experimental::filesystem;
#else
  error "Missing the <filesystem> header."
#endif
#include "version.hpp"
#include "idefix.hpp"
#include "dataBlock.hpp"
#include "gridHost.hpp"
#include "output.hpp"
#include "fluid.hpp"

#define VTK_RECTILINEAR_GRID    14
#define VTK_STRUCTURED_GRID     35

#ifndef VTK_FORMAT
  #if GEOMETRY == CARTESIAN || GEOMETRY == CYLINDRICAL
    #define VTK_FORMAT  VTK_RECTILINEAR_GRID
  #else
    #define VTK_FORMAT  VTK_STRUCTURED_GRID
  #endif
#endif

void Vtk::WriteHeaderNodes(IdfxFileHandler fvtk) {
  int64_t size = node_coord.extent(0) *
             node_coord.extent(1) *
             node_coord.extent(2) *
             node_coord.extent(3);
  if(size> INT_MAX) {
    IDEFIX_WARNING("Possible overflow in I/O routine");
  }
#ifdef WITH_MPI
  int size_int = static_cast<int>(size);
  MPI_SAFE_CALL(MPI_File_set_view(fvtk, this->offset, MPI_FLOAT, this->nodeView,
                                  "native", MPI_INFO_NULL));
  MPI_SAFE_CALL(MPI_File_write_all(fvtk, node_coord.data(), size_int,
                                   MPI_FLOAT, MPI_STATUS_IGNORE));
  this->offset += sizeof(float)*(nx1+ioffset)*(nx2+joffset)*(nx3+koffset)*3;
#else
  if(fwrite(node_coord.data(),sizeof(float),size,fvtk) != size) {
    IDEFIX_ERROR("Unable to write to file. Check your filesystem permissions and disk quota.");
  }
#endif
}

/*init the object */
Vtk::Vtk(Input &input, DataBlock *datain, std::string filebase) {
  // Initialize the output structure
  // Create a local datablock as an image of gridin
  this->data = datain;

  this->filebase = filebase;

  // Pointer to global grid
  GridHost grid(*(datain->mygrid));
  grid.SyncFromDevice();

  // initialize output path
  if(input.CheckEntry("Output","vtk_dir")>=0) {
    outputDirectory = input.Get<std::string>("Output","vtk_dir",0);
  } else {
    outputDirectory = "./";
  }

  if(idfx::prank==0) {
    if(!fs::is_directory(outputDirectory)) {
      try {
        if(!fs::create_directory(outputDirectory)) {
          std::stringstream msg;
          msg << "Cannot create directory " << outputDirectory << std::endl;
          IDEFIX_ERROR(msg);
        }
      } catch(std::exception &e) {
        std::stringstream msg;
        msg << "Cannot create directory " << outputDirectory << std::endl;
        msg << e.what();
        IDEFIX_ERROR(msg);
      }
    }
  }

  /* Note that there are two kinds of dimensions:
     - nx1, nx2, nx3, derived from the grid, which are the global dimensions
     - nx1loc,nx2loc,n3loc, which are the local dimensions of the current datablock
  */

  for (int dir=0; dir<3; dir++) {
    this->periodicity[dir] = (datain->mygrid->lbound[dir] == periodic);
  }
  // Create the coordinate array required in VTK files
  this->nx1 = grid.np_int[IDIR];
  this->nx2 = grid.np_int[JDIR];
  this->nx3 = grid.np_int[KDIR];

  this->nx1loc = data->np_int[IDIR];
  this->nx2loc = data->np_int[JDIR];
  this->nx3loc = data->np_int[KDIR];

  this->ioffset = datain->mygrid->np_tot[IDIR] == 1 ? 0 : 1;
  this->joffset = datain->mygrid->np_tot[JDIR] == 1 ? 0 : 1;
  this->koffset = datain->mygrid->np_tot[KDIR] == 1 ? 0 : 1;

  // Temporary storage on host for 3D arrays
  this->vect3D = new float[nx1loc*nx2loc*nx3loc];

  // Store coordinates for later use
  this->xnode = new float[nx1+ioffset];
  this->ynode = new float[nx2+joffset];
  this->znode = new float[nx3+koffset];
  this->xcenter = new float[nx1];
  this->ycenter = new float[nx2];
  this->zcenter = new float[nx3];

  for (int32_t i = 0; i < nx1 + ioffset; i++) {
    if(grid.np_tot[IDIR] == 1) { // only one dimension in this direction
      xnode[i] = bigEndian(static_cast<float>(grid.x[IDIR](i)));
    } else {
      xnode[i] = bigEndian(static_cast<float>(grid.xl[IDIR](i + grid.nghost[IDIR])));
    }
  }
  for (int32_t j = 0; j < nx2 + joffset; j++)    {
    if(grid.np_tot[JDIR] == 1) { // only one dimension in this direction
      ynode[j] = bigEndian(static_cast<float>(grid.x[JDIR](j)));
    } else {
      ynode[j] = bigEndian(static_cast<float>(grid.xl[JDIR](j + grid.nghost[JDIR])));
    }
  }
  for (int32_t k = 0; k < nx3 + koffset; k++) {
    if(grid.np_tot[KDIR] == 1) {
      znode[k] = bigEndian(static_cast<float>(grid.x[KDIR](k)));
    } else {
      znode[k] = bigEndian(static_cast<float>(grid.xl[KDIR](k + grid.nghost[KDIR])));
    }
  }
  for (int32_t i = 0; i < nx1; i++) {
    if(grid.np_tot[IDIR] == 1) { // only one dimension in this direction
      xcenter[i] = xnode[i];
    } else {
      xcenter[i] = bigEndian(static_cast<float>(grid.x[IDIR](i + grid.nghost[IDIR])));
    }
  }
  for (int32_t j = 0; j < nx2; j++)    {
    if(grid.np_tot[JDIR] == 1) { // only one dimension in this direction
      ycenter[j] = ynode[j];
    } else {
      ycenter[j] = bigEndian(static_cast<float>(grid.x[JDIR](j + grid.nghost[JDIR])));
    }
  }
  for (int32_t k = 0; k < nx3; k++) {
    if(grid.np_tot[KDIR] == 1) {
      zcenter[k] = znode[k];
    } else {
      zcenter[k] = bigEndian(static_cast<float>(grid.x[KDIR](k + grid.nghost[KDIR])));
    }
  }
#if VTK_FORMAT == VTK_STRUCTURED_GRID   // VTK_FORMAT
  /* -- Allocate memory for node_coord which is later used -- */

  // initialize node array dimensions
  [[maybe_unused]] int nodestart[4];
  int nodesize[4];
  int nodesubsize[4];

  for(int dir = 0; dir < 3 ; dir++) {
    nodesize[2-dir] = datain->mygrid->np_int[dir];
    nodestart[2-dir] = datain->gbeg[dir]-datain->nghost[dir];
    nodesubsize[2-dir] = datain->np_int[dir];
  }

  // In the 4th dimension, we always have the 3 components
  nodesize[3] = 3;
  nodestart[3] = 0;
  nodesubsize[3] = 3;

  // Since we use cell-defined vtk variables,
  // we add one cell in each direction when we're looking at the last
  // sub-domain in each direction
  nodesize[2] += ioffset;
  nodesize[1] += joffset;
  nodesize[0] += koffset;

  if(datain->mygrid->xproc[0] == datain->mygrid->nproc[0]-1) nodesubsize[2] += ioffset;
  if(datain->mygrid->xproc[1] == datain->mygrid->nproc[1]-1) nodesubsize[1] += joffset;
  if(datain->mygrid->xproc[2] == datain->mygrid->nproc[2]-1) nodesubsize[0] += koffset;

  // Build an MPI view if needed
  #ifdef WITH_MPI
    // Keep communicator for later use
    MPI_SAFE_CALL(MPI_Type_create_subarray(4, nodesize, nodesubsize, nodestart,
                                          MPI_ORDER_C, MPI_FLOAT, &this->nodeView));
    MPI_SAFE_CALL(MPI_Type_commit(&this->nodeView));
  #endif

  // Allocate a node view on the host
  node_coord = IdefixHostArray4D<float>("VtkNodeCoord",nodesubsize[0],
                                                       nodesubsize[1],
                                                       nodesubsize[2],
                                                       nodesubsize[3]);

  // fill the node_coord array
  float x1 = 0.0;
  [[maybe_unused]] float x2 = 0.0;
  [[maybe_unused]] float x3 = 0.0;
  for (int32_t k = 0; k < nodesubsize[0]; k++) {
    for (int32_t j = 0; j < nodesubsize[1]; j++) {
      for (int32_t i = 0; i < nodesubsize[2]; i++) {
        // bigEndian allows us to get back to little endian when needed
          x1 = grid.xl[IDIR](i + data->gbeg[IDIR]);
          x2 = grid.xl[JDIR](j + data->gbeg[JDIR]);
          x3 = grid.xl[KDIR](k + data->gbeg[KDIR]);

  #if (GEOMETRY == CARTESIAN) || (GEOMETRY == CYLINDRICAL)
        node_coord(k,j,i,0) = bigEndian(x1);
        node_coord(k,j,i,1) = bigEndian(x2);
        node_coord(k,j,i,2) = bigEndian(x3);

  #elif GEOMETRY == POLAR
        node_coord(k,j,i,0) = bigEndian(x1 * std::cos(x2));
        node_coord(k,j,i,1) = bigEndian(x1 * std::sin(x2));
        node_coord(k,j,i,2) = bigEndian(x3);

  #elif GEOMETRY == SPHERICAL
    #if DIMENSIONS == 1
        node_coord(k,j,i,0) = bigEndian(x1);
        node_coord(k,j,i,1) = bigEndian(0.0f);
        node_coord(k,j,i,2) = bigEndian(0.0f);
    #elif DIMENSIONS == 2
        node_coord(k,j,i,0) = bigEndian(x1 * std::sin(x2));
        node_coord(k,j,i,1) = bigEndian(x1 * std::cos(x2));
        node_coord(k,j,i,2) = bigEndian(0.0f);

    #elif DIMENSIONS == 3
        node_coord(k,j,i,0) = bigEndian(x1 * std::sin(x2) * std::cos(x3));
        node_coord(k,j,i,1) = bigEndian(x1 * std::sin(x2) * std::sin(x3));
        node_coord(k,j,i,2) = bigEndian(x1 * std::cos(x2));
    #endif // DIMENSIONS
  #endif // GEOMETRY
      }
    }
  }

#endif // VTK_FORMAT == VTK_STRUCTURED_GRID

  // Create MPI view when using MPI I/O
#ifdef WITH_MPI
  int start[3];
  int size[3];
  int subsize[3];

  for(int dir = 0; dir < 3 ; dir++) {
    // VTK assumes Fortran array ordering, hence arrays dimensions are filled backwards
    start[2-dir] = data->gbeg[dir]-grid.nghost[dir];
    size[2-dir] = grid.np_int[dir];
    subsize[2-dir] = data->np_int[dir];
  }

  MPI_SAFE_CALL(MPI_Type_create_subarray(3, size, subsize, start, MPI_ORDER_C,
                                         MPI_FLOAT, &this->view));
  MPI_SAFE_CALL(MPI_Type_commit(&this->view));
  this->comm = datain->mygrid->CartComm;
  this->isRoot =   (data->mygrid->xproc[0] == 0)
                && (data->mygrid->xproc[1] == 0)
                && (data->mygrid->xproc[2] == 0);
#endif

  // Register variables that are required in restart dumps
  if(data->dump.get() != nullptr)
    data->dump->RegisterVariable(&vtkFileNumber, "vtkFileNumber");
}


int Vtk::Write() {
  idfx::pushRegion("Vtk::Write");

  IdfxFileHandler fileHdl;
  fs::path filename;

  timer.reset();

  std::stringstream ssfileName, ssvtkFileNum;
  ssvtkFileNum << std::setfill('0') << std::setw(4) << vtkFileNumber;
  ssfileName << filebase << "." << ssvtkFileNum.str() << ".vtk";
  filename = outputDirectory/ssfileName.str();

  idfx::cout << "Vtk: Write file " << ssfileName.str() << "..." << std::flush;

  // Check if file exists, if yes, delete it
  if(this->isRoot) {
    if(fs::exists(filename)) {
      fs::remove(filename);
    }
  }


  // Open file and write header
#ifdef WITH_MPI
  MPI_Barrier(this->comm);
  // Open file for creating, return error if file already exists.
  MPI_SAFE_CALL(MPI_File_open(this->comm, filename.c_str(),
                              MPI_MODE_CREATE | MPI_MODE_RDWR
                              | MPI_MODE_EXCL | MPI_MODE_UNIQUE_OPEN,
                              MPI_INFO_NULL, &fileHdl));
  this->offset = 0;
#else
  fileHdl = fopen(filename.c_str(),"wb");

  if(fileHdl == NULL) {
    std::stringstream msg;
    msg << "Unable to open file " << filename << std::endl;
    msg << "Check that you have write access and that you don't exceed your quota." << std::endl;
    IDEFIX_ERROR(msg);
  }
#endif

  WriteHeader(fileHdl, this->data->t);

  // Write field one by one
  for(auto const& [name, scalar] : vtkScalarMap) {
    auto Vcin = scalar.GetHostField();
    for(int k = data->beg[KDIR]; k < data->end[KDIR] ; k++ ) {
      for(int j = data->beg[JDIR]; j < data->end[JDIR] ; j++ ) {
        for(int i = data->beg[IDIR]; i < data->end[IDIR] ; i++ ) {
          vect3D[i-data->beg[IDIR] + (j-data->beg[JDIR])*nx1loc + (k-data->beg[KDIR])*nx1loc*nx2loc]
              = bigEndian(static_cast<float>(Vcin(k,j,i)));
        }
      }
    }
    WriteScalar(fileHdl, vect3D, name);
  }

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
void Vtk::WriteHeader(IdfxFileHandler fvtk, real time) {
/*!
* Write VTK header in parallel or serial mode.
*
* \param [in]  fvtk  pointer to file
* \param [in]  time  time of the simulation in code units
*
* \todo  Write the grid using several processors.
*********************************************************************** */

  std::string header;
  std::stringstream ssheader;

  /* -------------------------------------------
  1. File version and identifier
  ------------------------------------------- */

  ssheader << "# vtk DataFile Version 2.0" << std::endl;

  /* -------------------------------------------
  2. Header
  ------------------------------------------- */

  ssheader << "Idefix " << IDEFIX_VERSION << " VTK Data" << std::endl;

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
  // fields: geometry, periodicity, time, 6 NativeCoordinates (x1l, x2l, x3l, x1c, x2c, x3c)
  int nfields = 9;

  // Write grid geometry in the VTK file
  ssheader << "FIELD FieldData " << nfields << std::endl;
  ssheader << "GEOMETRY 1 1 int" << std::endl;

  // Flush the ascii header
  header = ssheader.str();
  WriteHeaderString(header.c_str(), fvtk);
  // reset the string stream
  ssheader.str(std::string());

  // convert time to single precision big endian
  int32_t geoBig = bigEndian(this->geometry);

  WriteHeaderBinary(&geoBig, 1, fvtk);
  // Done, add cariage return for next ascii write
  ssheader << std::endl;

  // write grid periodicity
  ssheader << "PERIODICITY 1 3 int" << std::endl;

  // Flush the ascii header
  header = ssheader.str();
  WriteHeaderString(header.c_str(), fvtk);
  // reset the string stream
  ssheader.str(std::string());

  int32_t perBig{-1};
  for (int dir=0; dir<3; dir++) {
    perBig = bigEndian(this->periodicity[dir]);
    WriteHeaderBinary(&perBig, 1, fvtk);
  }
  // Done, add cariage return for next ascii write
  ssheader << std::endl;

  ssheader << "TIME 1 1 float" << std::endl;
  // Flush the ascii header
  header = ssheader.str();
  WriteHeaderString(header.c_str(), fvtk);
  // reset the string stream
  ssheader.str(std::string());

  // convert time to single precision big endian
  float timeBE = bigEndian(static_cast<float>(time));

  WriteHeaderBinary(&timeBE, 1, fvtk);
  // Done, add cariage return for next ascii write
  ssheader << std::endl;

  // write x1l native coordinates
  ssheader << "X1L_NATIVE_COORDINATES 1 " << (nx1 + ioffset) << " float" << std::endl;
  // Flush the ascii header
  header = ssheader.str();
  WriteHeaderString(header.c_str(), fvtk);
  // reset the string stream
  ssheader.str(std::string());
  WriteHeaderBinary(this->xnode, nx1 + ioffset, fvtk);
  // Done, add cariage return for next ascii write
  ssheader << std::endl;

  // write x2l native coordinates
  ssheader << "X2L_NATIVE_COORDINATES 1 " << (nx2 + joffset) << " float" << std::endl;
  // Flush the ascii header
  header = ssheader.str();
  WriteHeaderString(header.c_str(), fvtk);
  // reset the string stream
  ssheader.str(std::string());
  WriteHeaderBinary(this->ynode, nx2 + joffset, fvtk);
  // Done, add cariage return for next ascii write
  ssheader << std::endl;

  // write x3l native coordinates
  ssheader << "X3L_NATIVE_COORDINATES 1 " << (nx3 + koffset) << " float" << std::endl;
  // Flush the ascii header
  header = ssheader.str();
  WriteHeaderString(header.c_str(), fvtk);
  // reset the string stream
  ssheader.str(std::string());
  WriteHeaderBinary(this->znode, nx3 + koffset, fvtk);
  // Done, add cariage return for next ascii write
  ssheader << std::endl;

  // write x1 native coordinates
  ssheader << "X1C_NATIVE_COORDINATES 1 " << nx1 << " float" << std::endl;
  // Flush the ascii header
  header = ssheader.str();
  WriteHeaderString(header.c_str(), fvtk);
  // reset the string stream
  ssheader.str(std::string());
  WriteHeaderBinary(this->xcenter, nx1, fvtk);
  // Done, add cariage return for next ascii write
  ssheader << std::endl;

  // write x2 native coordinates
  ssheader << "X2C_NATIVE_COORDINATES 1 " << nx2 << " float" << std::endl;
  // Flush the ascii header
  header = ssheader.str();
  WriteHeaderString(header.c_str(), fvtk);
  // reset the string stream
  ssheader.str(std::string());
  WriteHeaderBinary(this->ycenter, nx2, fvtk);
  // Done, add cariage return for next ascii write
  ssheader << std::endl;

  // write x3 native coordinates
  ssheader << "X3C_NATIVE_COORDINATES 1 " << nx3 << " float" << std::endl;
  // Flush the ascii header
  header = ssheader.str();
  WriteHeaderString(header.c_str(), fvtk);
  // reset the string stream
  ssheader.str(std::string());
  WriteHeaderBinary(this->zcenter, nx3, fvtk);
  // Done, add cariage return for next ascii write
  ssheader << std::endl;

  ssheader << "DIMENSIONS " << nx1 + ioffset << " " << nx2 + joffset << " " << nx3 + koffset
           << std::endl;

  header = ssheader.str();
  WriteHeaderString(header.c_str(), fvtk);

#if VTK_FORMAT == VTK_RECTILINEAR_GRID

  /* -- Write rectilinear grid information -- */
  std::stringstream coordx, coordy, coordz;

  coordx << "X_COORDINATES " << nx1 + ioffset << " float" << std::endl;
  header = coordx.str();
  WriteHeaderString(header.c_str(), fvtk);

  WriteHeaderBinary(xnode, nx1 + ioffset, fvtk);

  coordy << std::endl << "Y_COORDINATES " << nx2 + joffset << " float" << std::endl;
  header = coordy.str();
  WriteHeaderString(header.c_str(), fvtk);

  WriteHeaderBinary(ynode, nx2 + joffset, fvtk);

  coordz << std::endl << "Z_COORDINATES " << nx3 + koffset << " float" << std::endl;
  header = coordz.str();
  WriteHeaderString(header.c_str(), fvtk);

  WriteHeaderBinary(znode, nx3 + koffset, fvtk);

#elif VTK_FORMAT == VTK_STRUCTURED_GRID

  /* -- define node_coord -- */
  std::stringstream ssheader2;

  ssheader2 << "POINTS " << (nx1 + ioffset) * (nx2 + joffset) * (nx3 + koffset) << " float"
            << std::endl;
  header = ssheader2.str();
  WriteHeaderString(header.c_str(), fvtk);

  /* -- Write structured grid information -- */

  WriteHeaderNodes(fvtk);

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
void Vtk::WriteScalar(IdfxFileHandler fvtk, float* Vin,  const std::string &var_name) {
/*!
* Write VTK scalar field.
*
*********************************************************************** */

  std::stringstream ssheader;

  ssheader << std::endl << "SCALARS " << var_name.c_str() << " float" << std::endl;
  ssheader << "LOOKUP_TABLE default" << std::endl;
  std::string header(ssheader.str());

  WriteHeaderString(header.c_str(), fvtk);

#ifdef WITH_MPI
  MPI_SAFE_CALL(MPI_File_set_view(fvtk, this->offset, MPI_FLOAT, this->view,
                                  "native", MPI_INFO_NULL));

  int nwrite = nx1loc*nx2loc*nx3loc;
  //if(idfx::prank != 0) nwrite = 0;
  MPI_SAFE_CALL(MPI_File_write_all(fvtk, Vin, nwrite, MPI_FLOAT, MPI_STATUS_IGNORE));

  this->offset = this->offset + sizeof(float)*nx1*nx2*nx3;
#else
  if(fwrite(Vin,sizeof(float),nx1loc*nx2loc*nx3loc,fvtk) != nx1loc*nx2loc*nx3loc) {
    IDEFIX_ERROR("Unable to write to file. Check your filesystem permissions and disk quota.");
  }
#endif
}
