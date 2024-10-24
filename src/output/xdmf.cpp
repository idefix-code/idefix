// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <iomanip>
#if __has_include(<filesystem>)
  #include <filesystem> // NOLINT [build/c++17]
  namespace fs = std::filesystem;
#elif __has_include(<experimental/filesystem>)
  #include <experimental/filesystem>
  namespace fs = std::experimental::filesystem;
#else
  #error "Missing the <filesystem> header."
#endif

#include "xdmf.hpp"
#include "version.hpp"
#include "idefix.hpp"
#include "dataBlockHost.hpp"
#include "gridHost.hpp"
#include "output.hpp"

// Whether or not we write the time in the XDMF file
#define WRITE_TIME


Xdmf::Xdmf(Input &input, DataBlock *datain) {
  // Initialize the output structure
  // Create a local datablock as an image of gridin

  this->data = datain;

  // Pointer to global grid
  GridHost grid(*(data->mygrid));
  grid.SyncFromDevice();

  // initialize output path
  if(input.CheckEntry("Output","xdmf_dir")>=0) {
    outputDirectory = input.Get<std::string>("Output","xdmf_dir",0);
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
    this->periodicity[dir] = (data->mygrid->lbound[dir] == periodic);
  }
  // Create the coordinate array required in XDMF files
  this->nx1 = grid.np_int[IDIR];
  this->nx2 = grid.np_int[JDIR];
  this->nx3 = grid.np_int[KDIR];

  this->nx1loc = data->np_int[IDIR];
  this->nx2loc = data->np_int[JDIR];
  this->nx3loc = data->np_int[KDIR];

  this->nx1tot = grid.np_tot[IDIR];
  this->nx2tot = grid.np_tot[JDIR];
  this->nx3tot = grid.np_tot[KDIR];

  this->nx1loctot = data->np_tot[IDIR];
  this->nx2loctot = data->np_tot[JDIR];
  this->nx3loctot = data->np_tot[KDIR];


  this->ngx1 = grid.nghost[IDIR];
  this->ngx2 = grid.nghost[JDIR];
  this->ngx3 = grid.nghost[KDIR];

  // Store coordinates for later use
  this->xnode = new DUMP_DATATYPE[nx1+IOFFSET];
  this->ynode = new DUMP_DATATYPE[nx2+JOFFSET];
  this->znode = new DUMP_DATATYPE[nx3+KOFFSET];
  this->xcell = new DUMP_DATATYPE[nx1];
  this->ycell = new DUMP_DATATYPE[nx2];
  this->zcell = new DUMP_DATATYPE[nx3];

  for (int32_t i = 0; i < nx1 + IOFFSET; i++) {
    xnode[i] = static_cast<DUMP_DATATYPE>(grid.xl[IDIR](i + grid.nghost[IDIR]));
    if (i<nx1) xcell[i] = static_cast<DUMP_DATATYPE>(grid.x[IDIR](i + grid.nghost[IDIR]));
  }
  for (int32_t j = 0; j < nx2 + JOFFSET; j++) {
    if(DIMENSIONS==1) {
      ynode[j] = static_cast<DUMP_DATATYPE>(0.0);
      if (j<nx2) ycell[j] = static_cast<DUMP_DATATYPE>(0.0);
    } else {
      ynode[j] = static_cast<DUMP_DATATYPE>(grid.xl[JDIR](j + grid.nghost[JDIR]));
      if (j<nx2) ycell[j] = static_cast<DUMP_DATATYPE>(grid.x[JDIR](j + grid.nghost[JDIR]));
    }
  }
  for (int32_t k = 0; k < nx3 + KOFFSET; k++) {
    if(DIMENSIONS<3) {
      znode[k] = static_cast<DUMP_DATATYPE>(0.0);
      if (k<nx3) zcell[k] = static_cast<DUMP_DATATYPE>(0.0);
    } else {
      znode[k] = static_cast<DUMP_DATATYPE>(grid.xl[KDIR](k + grid.nghost[KDIR]));
      if (k<nx3) zcell[k] = static_cast<DUMP_DATATYPE>(grid.x[KDIR](k + grid.nghost[KDIR]));
    }
  }

  /* -- Allocate memory for node_coord which is later used -- */
  /* -- Data order that is saved is 3D/1D: Z-Y-X and 2D: Y-X-Z -- */
  for(int dir = 0; dir < 3 ; dir++) {
    this->nodesize[3-dir] = data->mygrid->np_int[dir];
    this->nodestart[3-dir] = data->gbeg[dir]-data->nghost[dir];
    this->nodesubsize[3-dir] = data->np_int[dir];

    this->cellsize[3-dir] = data->mygrid->np_int[dir];
    this->cellstart[3-dir] = data->gbeg[dir]-data->nghost[dir];
    this->cellsubsize[3-dir] = data->np_int[dir];
  }

  // In the 0th dimension, we always have the 3 components
  // as this is meshgrid on 3D arrays
  this->nodesize[0] = 3;
  this->nodestart[0] = 0;
  this->nodesubsize[0] = 3;

  this->cellsize[0] = 3;
  this->cellstart[0] = 0;
  this->cellsubsize[0] = 3;

  // Since we use cell-defined xdmf variables,
  // we add one cell in each direction when we're looking at the last
  // sub-domain in each direction
  this->nodesize[3] += IOFFSET;
  this->nodesize[2] += JOFFSET;
  this->nodesize[1] += KOFFSET;

  if(data->mygrid->xproc[0] == data->mygrid->nproc[0]-1) this->nodesubsize[3] += IOFFSET;
  if(data->mygrid->xproc[1] == data->mygrid->nproc[1]-1) this->nodesubsize[2] += JOFFSET;
  if(data->mygrid->xproc[2] == data->mygrid->nproc[2]-1) this->nodesubsize[1] += KOFFSET;

  // Allocate a node and cell views on the host
  node_coord = IdefixHostArray4D<DUMP_DATATYPE>("XdmfNodeCoord", nodesubsize[0],
                                                                 nodesubsize[1],
                                                                 nodesubsize[2],
                                                                 nodesubsize[3]);

  cell_coord = IdefixHostArray4D<DUMP_DATATYPE>("XdmfCellCoord", cellsubsize[0],
                                                                 cellsubsize[1],
                                                                 cellsubsize[2],
                                                                 cellsubsize[3]);
  /*
  field_data = IdefixHostArray3D<DUMP_DATATYPE>("XdmfFieldData", cellsubsize[1],
                                                                 cellsubsize[2],
                                                                 cellsubsize[3]);
  */
  // Temporary storage on host for 3D arrays
  this->vect3D = new DUMP_DATATYPE[nx1loc*nx2loc*nx3loc];

  // fill the node_coord array
  DUMP_DATATYPE x1 = 0.0;
  [[maybe_unused]] DUMP_DATATYPE x2 = 0.0;
  [[maybe_unused]] DUMP_DATATYPE x3 = 0.0;
  DUMP_DATATYPE x1_cell = 0.0;
  [[maybe_unused]] DUMP_DATATYPE x2_cell = 0.0;
  [[maybe_unused]] DUMP_DATATYPE x3_cell = 0.0;

  for (int32_t k = 0; k < nodesubsize[1]; k++) {
    for (int32_t j = 0; j < nodesubsize[2]; j++) {
      for (int32_t i = 0; i < nodesubsize[3]; i++) {
        D_EXPAND( x1 = grid.xl[IDIR](i + data->gbeg[IDIR]);  ,
                  x2 = grid.xl[JDIR](j + data->gbeg[JDIR]);  ,
                  x3 = grid.xl[KDIR](k + data->gbeg[KDIR]);  )
        if ( (k<(cellsubsize[1])) && (j<(cellsubsize[2])) && (i<(cellsubsize[3])) ) {
          D_EXPAND( x1_cell = grid.x[IDIR](i + data->gbeg[IDIR]);  ,
                    x2_cell = grid.x[JDIR](j + data->gbeg[JDIR]);  ,
                    x3_cell = grid.x[KDIR](k + data->gbeg[KDIR]);  )
        }
        #if (GEOMETRY == CARTESIAN) || (GEOMETRY == CYLINDRICAL)
        node_coord(0,k,j,i) = x1;
        node_coord(1,k,j,i) = x2;
        node_coord(2,k,j,i) = x3;
        if ( (k<(cellsubsize[1])) && (j<(cellsubsize[2])) && (i<(cellsubsize[3])) ) {
          cell_coord(0,k,j,i) = x1_cell;
          cell_coord(1,k,j,i) = x2_cell;
          cell_coord(2,k,j,i) = x3_cell;
        }

        #elif GEOMETRY == POLAR
        node_coord(0,k,j,i) = x1 * cos(x2);
        node_coord(1,k,j,i) = x1 * sin(x2);
        node_coord(2,k,j,i) = x3;
        if ( (k<(cellsubsize[1])) && (j<(cellsubsize[2])) && (i<(cellsubsize[3])) ) {
          cell_coord(0,k,j,i) = x1_cell* cos(x2_cell);
          cell_coord(1,k,j,i) = x1_cell * sin(x2_cell);
          cell_coord(2,k,j,i) = x3_cell;
        }

        #elif GEOMETRY == SPHERICAL
          #if DIMENSIONS == 1
        node_coord(0,k,j,i) = x1;
        node_coord(1,k,j,i) = 0.0;
        node_coord(2,k,j,i) = 0.0;
        if ( (k<(cellsubsize[1])) && (j<(cellsubsize[2])) && (i<(cellsubsize[3])) ) {
          cell_coord(0,k,j,i) = x1_cell;
          cell_coord(1,k,j,i) = 0.0;
          cell_coord(2,k,j,i) = 0.0;
        }
          #elif DIMENSIONS == 2
        node_coord(0,k,j,i) = x1 * sin(x2);
        node_coord(1,k,j,i) = x1 * cos(x2);
        node_coord(2,k,j,i) = 0.0;
        if ( (k<(cellsubsize[1])) && (j<(cellsubsize[2])) && (i<(cellsubsize[3])) ) {
          cell_coord(0,k,j,i) = x1_cell * sin(x2_cell);
          cell_coord(1,k,j,i) = x1_cell * cos(x2_cell);
          cell_coord(2,k,j,i) = 0.0;
        }
          #elif DIMENSIONS == 3
        node_coord(0,k,j,i) = x1 * sin(x2) * cos(x3);
        node_coord(1,k,j,i) = x1 * sin(x2) * sin(x3);
        node_coord(2,k,j,i) = x1 * cos(x2);
        if ( (k<(cellsubsize[1])) && (j<(cellsubsize[2])) && (i<(cellsubsize[3])) ) {
          cell_coord(0,k,j,i) = x1_cell * sin(x2_cell) * cos(x3_cell);
          cell_coord(1,k,j,i) = x1_cell * sin(x2_cell) * sin(x3_cell);
          cell_coord(2,k,j,i) = x1_cell * cos(x2_cell);
        }
          #endif // DIMENSIONS
        #endif // GEOMETRY
      }
    }
  }

  // Create MPI view when using MPI I/O
  #ifdef WITH_MPI
  #if (DIMENSIONS == 1) || (DIMENSIONS == 3)
  for(int dir = 0; dir < 3 ; dir++) {
    // XDMF assumes Fortran array ordering, hence arrays dimensions are filled backwards
    // So ordering is 3D/1D: Z-Y-X and 2D: Y-X-Z
    // offset in the destination array
    this->mpi_data_start[dir] = data->gbeg[2-dir]-grid.nghost[2-dir];
    this->mpi_data_size[dir] = grid.np_int[2-dir];
    this->mpi_data_subsize[dir] = data->np_int[2-dir];
  }
  #elif (DIMENSIONS == 2)
  for(int dir = 0; dir < DIMENSIONS ; dir++) {
    // XDMF assumes Fortran array ordering, hence arrays dimensions are filled backwards
    // So ordering is 3D/1D: Z-Y-X and 2D: Y-X-Z
    // offset in the destination array
    this->mpi_data_start[dir] = data->gbeg[DIMENSIONS-dir-1]-grid.nghost[DIMENSIONS-dir-1];
    this->mpi_data_size[dir] = grid.np_int[DIMENSIONS-dir-1];
    this->mpi_data_subsize[dir] = data->np_int[DIMENSIONS-dir-1];
  }
  for(int dir = DIMENSIONS; dir < 3 ; dir++) {
    // XDMF assumes Fortran array ordering, hence arrays dimensions are filled backwards
    // So ordering is 3D/1D: Z-Y-X and 2D: Y-X-Z
    // offset in the destination array
    this->mpi_data_start[dir] = data->gbeg[dir]-grid.nghost[dir];
    this->mpi_data_size[dir] = grid.np_int[dir];
    this->mpi_data_subsize[dir] = data->np_int[dir];
  }
  #endif
  #endif
}

int Xdmf::Write() {
  idfx::pushRegion("Xdmf::Write");
  fs::path filename;
  fs::path filename_xmf;
  hid_t err;
  idfx::cout << "Xdmf: Write file n " << xdmfFileNumber << "..." << std::flush;
  timer.reset();

  // Create a copy of the dataBlock on Host, and sync it.

  #if DIMENSIONS == 1
  int tot_dim = 1;
  #elif DIMENSIONS == 2
  int tot_dim = 2;
  #elif DIMENSIONS == 3
  int tot_dim = 3;
  #endif

  std::stringstream ssfileName, ssfileNameXmf, ssxdmfFileNum;
  #ifndef XDMF_DOUBLE
  std::string extension = ".flt";
  #else
  std::string extension = ".dbl";
  #endif
  ssxdmfFileNum << std::setfill('0') << std::setw(4) << xdmfFileNumber;
  ssfileName << "data." << ssxdmfFileNum.str() << extension << ".h5";
  ssfileNameXmf << "data." << ssxdmfFileNum.str() << extension << ".xmf";
  filename = outputDirectory/ssfileName.str();
  filename_xmf = outputDirectory/ssfileNameXmf.str();

  #ifdef WITH_MPI
  hid_t file_access = H5Pcreate(H5P_FILE_ACCESS);
  // #if MPI_POSIX == YES
  // H5Pset_fapl_mpiposix(file_access, MPI_COMM_WORLD, 1);
  // #else
  H5Pset_fapl_mpio(file_access,  MPI_COMM_WORLD, MPI_INFO_NULL);
  // #endif
  hid_t fileHdf = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, file_access);
  H5Pclose(file_access);
  #else
  hid_t fileHdf = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  #endif

  hid_t group_fields; // = static_cast<hid_t *>(malloc(sizeof(hid_t)));
  std::stringstream ssgroup_name;
  ssgroup_name << "/Timestep_" << xdmfFileNumber;
  hid_t timestep = H5Gcreate(fileHdf, ssgroup_name.str().c_str(), 0);

  WriteHeader(fileHdf, ssfileName.str(), filename_xmf, data->t, timestep, group_fields);

  /* ------------------------------------
      write cell-centered field data
   ------------------------------------ */

  // doesn't affect serial io
  hid_t plist_id_mpiio = 0; /* for collective MPI I/O */
  #ifdef WITH_MPI
  plist_id_mpiio = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id_mpiio,H5FD_MPIO_COLLECTIVE);
  #endif

  // Layout of the field data in memory
  hsize_t field_data_size[3], field_data_start[3], field_data_subsize[3], stride[3];
  #ifdef WITH_MPI
  for(int dir = 0; dir < 3 ; dir++) {
    field_data_size[dir] = static_cast<hsize_t>(this->mpi_data_size[dir]);
    field_data_subsize[dir] = static_cast<hsize_t>(this->mpi_data_subsize[dir]);
    field_data_start[dir] = static_cast<hsize_t>(this->mpi_data_start[dir]);
    stride[dir] = static_cast<hsize_t>(1);
  }
  #else
  #if (DIMENSIONS == 1) || (DIMENSIONS == 3)
  field_data_size[0] = static_cast<hsize_t>(this->nx3);
  field_data_size[1] = static_cast<hsize_t>(this->nx2);
  field_data_size[2] = static_cast<hsize_t>(this->nx1);
  field_data_subsize[0] = static_cast<hsize_t>(this->nx3loc);
  field_data_subsize[1] = static_cast<hsize_t>(this->nx2loc);
  field_data_subsize[2] = static_cast<hsize_t>(this->nx1loc);
  #elif DIMENSIONS == 2
  field_data_size[0] = static_cast<hsize_t>(this->nx2);
  field_data_size[1] = static_cast<hsize_t>(this->nx1);
  field_data_size[2] = static_cast<hsize_t>(this->nx3);
  field_data_subsize[0] = static_cast<hsize_t>(this->nx2loc);
  field_data_subsize[1] = static_cast<hsize_t>(this->nx1loc);
  field_data_subsize[2] = static_cast<hsize_t>(this->nx3loc);
  #endif
  for(int dir = 0; dir < 3 ; dir++) {
    field_data_start[dir] = static_cast<hsize_t>(0);
    stride[dir] = static_cast<hsize_t>(1);
  }
  #endif

  int rank = DIMENSIONS;
  hsize_t dimens[3], offset[3];
  hid_t dataspace = H5Screate_simple(rank, field_data_size, NULL);
  #ifdef WITH_MPI
  err = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
                            field_data_start, stride,
                            field_data_subsize, NULL);
  #endif

  #if (DIMENSIONS == 1) || (DIMENSIONS == 3)
  dimens[0] = nx3loc; dimens[1] = nx2loc; dimens[2] = nx1loc;
  #elif DIMENSIONS == 2
  dimens[0] = nx2loc; dimens[1] = nx1loc; dimens[2] = nx3loc;
  #endif
  hid_t memspace = H5Screate_simple(rank, dimens, NULL);

  offset[0] = 0; offset[1] = 0; offset[2] = 0;
  err = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset, stride, field_data_subsize, NULL);

  // Write field one by one
  for(auto const& [name, scalar] : xdmfScalarMap) {
    auto Vcin = scalar.GetHostField();
    for(int k = data->beg[KDIR]; k < data->end[KDIR] ; k++ ) {
      for(int j = data->beg[JDIR]; j < data->end[JDIR] ; j++ ) {
        for(int i = data->beg[IDIR]; i < data->end[IDIR] ; i++ ) {
          vect3D[i-data->beg[IDIR] + (j-data->beg[JDIR])*nx1loc + (k-data->beg[KDIR])*nx1loc*nx2loc]
              = static_cast<DUMP_DATATYPE>(Vcin(k,j,i));
        }
      }
    }
    WriteScalar(vect3D, name, field_data_size, ssfileName.str(), filename_xmf,
                memspace, dataspace, plist_id_mpiio, static_cast<hid_t&>(group_fields));
  }
  WriteFooter(ssfileName.str(), filename_xmf);

  #ifdef WITH_MPI
  H5Pclose(plist_id_mpiio);
  #endif
  H5Sclose(memspace);
  H5Sclose(dataspace);
  H5Gclose(group_fields); // Close group "vars"
  H5Gclose(timestep);
  H5Fclose(fileHdf);

  xdmfFileNumber++;
  // Make file number
  idfx::cout << "done in " << timer.seconds() << " s." << std::endl;
  idfx::popRegion();
  // One day, we will have a return code.
  return(0);
}


/* ********************************************************************* */
void Xdmf::WriteHeader(
                       hid_t fileHdf,
                       const std::string filename,
                       const std::string filename_xmf,
                       real time,
                       hid_t &timestep,
                       hid_t &group_fields) {
/*!
* Write XDMF header in parallel or serial mode.
*
* \param [in]  fileHdf  pointer to file
* \param [in]  filename  hdf5 filename
* \param [in]  filename_xmf  xmf filename
* \param [in]  time  simulation time in code units
*
*********************************************************************** */

  std::stringstream ssheader;
  int rank;
  std::vector<std::string> directions{"X", "Y", "Z"};

  hid_t dataspace, memspace, dataset;
  hid_t strspace, stratt, string_type;
  hid_t tspace, tattr;
  hid_t unit_info, unit_attr;
  hid_t group;
  hid_t file_access = 0;
  #ifdef WITH_MPI
  hid_t plist_id_mpiio = 0; /* for collective MPI I/O */
  #endif
  hid_t err;

  hsize_t dimstr;

  #ifdef WRITE_TIME
  tspace  = H5Screate(H5S_SCALAR);
  tattr   = H5Acreate(timestep, "time", H5T_NATIVE_DOUBLE, tspace, H5P_DEFAULT);
  err = H5Awrite(tattr, H5T_NATIVE_DOUBLE, &time);
  H5Aclose(tattr);
  H5Sclose(tspace);
  #endif

  unit_info = H5Screate(H5S_SCALAR);
  double unit;
  #if defined(UNIT_DENSITY) || defined(UNIT_MASS)
  #ifdef UNIT_DENSITY
  unit_attr   = H5Acreate(timestep, "density_unit", H5T_NATIVE_DOUBLE, unit_info, H5P_DEFAULT);
  unit = UNIT_DENSITY;
  err = H5Awrite(unit_attr, H5T_NATIVE_DOUBLE, &unit);
  #endif
  #ifdef UNIT_MASS
  unit_attr   = H5Acreate(timestep, "mass_unit", H5T_NATIVE_DOUBLE, unit_info, H5P_DEFAULT);
  unit = UNIT_MASS;
  err = H5Awrite(unit_attr, H5T_NATIVE_DOUBLE, &unit);
  #endif
  #else
  unit_attr   = H5Acreate(timestep, "density_unit", H5T_NATIVE_DOUBLE, unit_info, H5P_DEFAULT);
  unit = 1.0;
  err = H5Awrite(unit_attr, H5T_NATIVE_DOUBLE, &unit);
  #endif
  H5Aclose(unit_attr);
  H5Sclose(unit_info);

  unit_info = H5Screate(H5S_SCALAR);
  #if defined(UNIT_VELOCITY) || defined(UNIT_TIME)
  #ifdef UNIT_VELOCITY
  unit_attr   = H5Acreate(timestep, "velocity_unit", H5T_NATIVE_DOUBLE, unit_info, H5P_DEFAULT);
  unit = UNIT_VELOCITY;
  err = H5Awrite(unit_attr, H5T_NATIVE_DOUBLE, &unit);
  #endif
  #ifdef UNIT_TIME
  unit_attr   = H5Acreate(timestep, "time_unit", H5T_NATIVE_DOUBLE, unit_info, H5P_DEFAULT);
  unit = UNIT_TIME;
  err = H5Awrite(unit_attr, H5T_NATIVE_DOUBLE, &unit);
  #endif
  #else
  unit_attr   = H5Acreate(timestep, "velocity_unit", H5T_NATIVE_DOUBLE, unit_info, H5P_DEFAULT);
  unit = 1.0;
  err = H5Awrite(unit_attr, H5T_NATIVE_DOUBLE, &unit);
  #endif
  H5Aclose(unit_attr);
  H5Sclose(unit_info);

  unit_info = H5Screate(H5S_SCALAR);
  #ifdef UNIT_LENGTH
  unit_attr   = H5Acreate(timestep, "length_unit", H5T_NATIVE_DOUBLE, unit_info, H5P_DEFAULT);
  unit = UNIT_LENGTH;
  err = H5Awrite(unit_attr, H5T_NATIVE_DOUBLE, &unit);
  #else
  unit_attr   = H5Acreate(timestep, "length_unit", H5T_NATIVE_DOUBLE, unit_info, H5P_DEFAULT);
  unit = 1.0;
  err = H5Awrite(unit_attr, H5T_NATIVE_DOUBLE, &unit);
  #endif
  H5Aclose(unit_attr);
  H5Sclose(unit_info);

  dimstr = 1;

  ssheader << "Idefix " << IDEFIX_VERSION << " XDMF Data";
  strspace = H5Screate_simple(1, &dimstr, NULL);
  string_type = H5Tcopy(H5T_C_S1);
  H5Tset_size(string_type, strlen( ssheader.str().c_str() ));
  H5Tset_strpad(string_type, H5T_STR_SPACEPAD);
  stratt = H5Acreate(timestep, "version", string_type, strspace, H5P_DEFAULT);
  err = H5Awrite(stratt, string_type, ssheader.str().c_str());
  H5Aclose(stratt);
  H5Sclose(strspace);

  #ifndef XDMF_DOUBLE
  std::string datatype = "float";
  #else
  std::string datatype = "double";
  #endif
  strspace = H5Screate_simple(1, &dimstr, NULL);
  string_type = H5Tcopy(H5T_C_S1);
  H5Tset_size(string_type, strlen( datatype.c_str() ));
  H5Tset_strpad(string_type, H5T_STR_SPACEPAD);
  stratt = H5Acreate(timestep, "dump_datatype", string_type, strspace, H5P_DEFAULT);
  err = H5Awrite(stratt, string_type, datatype.c_str());
  H5Aclose(stratt);
  H5Sclose(strspace);

  #if GEOMETRY == CARTESIAN
  std::string geometry = "cartesian";
  #elif GEOMETRY == CYLINDRICAL
  std::string geometry = "cylindrical";
  #elif GEOMETRY == POLAR
  std::string geometry = "polar";
  #elif GEOMETRY == SPHERICAL
  std::string geometry = "spherical";
  #endif
  strspace = H5Screate_simple(1, &dimstr, NULL);
  string_type = H5Tcopy(H5T_C_S1);
  H5Tset_size(string_type, strlen( geometry.c_str() ));
  H5Tset_strpad(string_type, H5T_STR_SPACEPAD);
  stratt = H5Acreate(timestep, "geometry", string_type, strspace, H5P_DEFAULT);
  err = H5Awrite(stratt, string_type, geometry.c_str());
  H5Aclose(stratt);
  H5Sclose(strspace);

  /* Create group_fields "vars" (cell-centered vars) */
  group_fields = H5Gcreate(timestep, "vars", 0);

  /* Define "coords" attribute of group_fields "vars" */
  dimstr = 3;
  std::string coords_label = "/cell_coords/X /cell_coords/Y /cell_coords/Z ";
  strspace = H5Screate_simple(1, &dimstr, NULL);
  string_type = H5Tcopy(H5T_C_S1);
  H5Tset_size( string_type, strlen("/cell_coords/X ") );
  H5Tset_strpad(string_type, H5T_STR_SPACEPAD);
  stratt = H5Acreate(group_fields, "coords", string_type, strspace, H5P_DEFAULT);
  err = H5Awrite(stratt, string_type, coords_label.c_str());
  H5Aclose(stratt);
  H5Sclose(strspace);

  hsize_t dimens[DIMENSIONS];
  hsize_t start[DIMENSIONS];
  hsize_t stride[DIMENSIONS];
  hsize_t count[DIMENSIONS];

  /* Create group "cell_coords" (centered mesh) */
  group = H5Gcreate(fileHdf, "cell_coords", 0);

  #if (DIMENSIONS == 1) || (DIMENSIONS == 3)
  for (int dir = 0; dir < DIMENSIONS; dir++) {
    dimens[dir] = this->cellsize[dir+1];
  }
  #elif DIMENSIONS == 2
  for (int dir = 0; dir < DIMENSIONS; dir++) {
    dimens[dir] = this->cellsize[DIMENSIONS+dir];
  }
  #endif
  rank   = DIMENSIONS;
  dataspace = H5Screate_simple(rank, dimens, NULL);

  #if (DIMENSIONS == 1) || (DIMENSIONS == 3)
  for (int dir = 0; dir < DIMENSIONS; dir++) {
    start[dir]  = this->cellstart[dir+1];
    stride[dir] = 1;
    count[dir]  = this->cellsubsize[dir+1];
  }
  #elif DIMENSIONS == 2
  for (int dir = 0; dir < DIMENSIONS; dir++) {
    start[dir]  = this->cellstart[DIMENSIONS+dir];
    stride[dir] = 1;
    count[dir]  = this->cellsubsize[DIMENSIONS+dir];
  }
  #endif
  err = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
                            start, stride, count, NULL);

  memspace = H5Screate_simple(rank,count,NULL);

  /* ------------------------------------
       write cell centered mesh
   ------------------------------------ */
  #ifdef WITH_MPI
  plist_id_mpiio = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id_mpiio, H5FD_MPIO_COLLECTIVE);
  #endif

  DUMP_DATATYPE *cell_mesh;
  for (int dir = 0; dir < 3; dir++) {
    dataset = H5Dcreate(group, directions[dir].c_str(), H5_DUMP_DATATYPE, dataspace, H5P_DEFAULT);
    cell_mesh = Kokkos::subview (this->cell_coord,
                                            dir,
                                            Kokkos::ALL(),
                                            Kokkos::ALL(),
                                            Kokkos::ALL()).data();
    #ifdef WITH_MPI
    err = H5Dwrite(dataset, H5_DUMP_DATATYPE, memspace,
                   dataspace, plist_id_mpiio, cell_mesh);
    #else
    err = H5Dwrite(dataset, H5_DUMP_DATATYPE, memspace,
                   dataspace, H5P_DEFAULT, cell_mesh);
    #endif
    H5Dclose(dataset);
  }

  #ifdef WITH_MPI
  H5Pclose(plist_id_mpiio);
  #endif
  H5Sclose(memspace);
  H5Sclose(dataspace);
  H5Gclose(group); /* Close group "cell_coords" */

  /* Create group "node_coords" (node mesh) */
  group = H5Gcreate(fileHdf, "node_coords", 0);

  #if (DIMENSIONS == 1) || (DIMENSIONS == 3)
  for (int dir = 0; dir < DIMENSIONS; dir++) {
    dimens[dir] = this->nodesize[dir+1];
  }
  #elif DIMENSIONS == 2
  for (int dir = 0; dir < DIMENSIONS; dir++) {
    dimens[dir] = this->nodesize[DIMENSIONS+dir];
  }
  #endif

  dataspace = H5Screate_simple(rank, dimens, NULL);

  #if (DIMENSIONS == 1) || (DIMENSIONS == 3)
  for (int dir = 0; dir < DIMENSIONS; dir++) {
    start[dir]  = this->nodestart[dir+1];
    stride[dir] = 1;
    count[dir]  = this->nodesubsize[dir+1];
    // if (grid->rbound[dir] != 0) count[nd] += 1;
  }
  #elif DIMENSIONS == 2
  for (int dir = 0; dir < DIMENSIONS; dir++) {
    start[dir]  = this->nodestart[DIMENSIONS+dir];
    stride[dir] = 1;
    count[dir]  = this->nodesubsize[DIMENSIONS+dir];
    // if (grid->rbound[dir] != 0) count[nd] += 1;
  }
  #endif

  err = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
                            start, stride, count, NULL);

  // for (int dir = 0; dir < DIMENSIONS; dir++) {
  //  dimens[dir] = this->nodesubsize[dir+1];
    // if (grid->rbound[nr] != 0) dimens[nd] += 1;
  // }
  memspace = H5Screate_simple(rank,count,NULL);

/* ------------------------------------
          write node centered mesh
   ------------------------------------ */

  #ifdef WITH_MPI
  plist_id_mpiio = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id_mpiio, H5FD_MPIO_COLLECTIVE );
  #endif

  for (int dir = 0; dir < 3; dir++) {
    dataset = H5Dcreate(group, directions[dir].c_str(), H5_DUMP_DATATYPE, dataspace, H5P_DEFAULT);

    DUMP_DATATYPE *node_mesh = Kokkos::subview (this->node_coord,
                                              dir,
                                              Kokkos::ALL(),
                                              Kokkos::ALL(),
                                              Kokkos::ALL()).data();
    #ifdef WITH_MPI
    err = H5Dwrite(dataset, H5_DUMP_DATATYPE, memspace,
                   dataspace, plist_id_mpiio, node_mesh);
    #else
    err = H5Dwrite(dataset, H5_DUMP_DATATYPE, memspace,
                   dataspace, H5P_DEFAULT, node_mesh);
    #endif
    H5Dclose(dataset);
  }

  #ifdef WITH_MPI
  H5Pclose(plist_id_mpiio);
  #endif
  H5Sclose(memspace);
  H5Sclose(dataspace);
  H5Gclose(group); /* Close group "node_coords" */

  if (idfx::prank==0) {
    std::vector<std::string> directions {"X", "Y", "Z"};
    std::stringstream ssxmfcontent;
    /*      XMF file generation      */
    /* -------------------------------------------
    1. File version and identifier
    ------------------------------------------- */

    ssxmfcontent << "<?xml version=\"1.0\" ?>" << std::endl;
    ssxmfcontent << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>" << std::endl;
    ssxmfcontent << "<Xdmf Version=\"2.0\">" << std::endl;
    ssxmfcontent << " <Domain>" << std::endl;
    ssxmfcontent << "   <Grid Name=\"node_mesh\" GridType=\"Uniform\">" << std::endl;
    ssxmfcontent << "    <Time Value=\"" << time << "\"/>" << std::endl;
    #if DIMENSIONS == 1
    ssxmfcontent << "     <Topology TopologyType=\"1DSMesh\" NumberOfElements=\"";
    int tot_dim = 1;
    #elif DIMENSIONS == 2
    ssxmfcontent << "     <Topology TopologyType=\"2DSMesh\" NumberOfElements=\"";
    int tot_dim = 2;
    #elif DIMENSIONS == 3
    ssxmfcontent << "     <Topology TopologyType=\"3DSMesh\" NumberOfElements=\"";
    int tot_dim = 3;
    #endif
    for (int dir = 3-tot_dim; dir < 3; dir++) {
      ssxmfcontent << this->nodesize[dir+1];
      if (dir<2)
        ssxmfcontent << " ";
    }
    ssxmfcontent << "\"/>" << std::endl;
    ssxmfcontent << "     <Geometry GeometryType=\"";
    for (int dir = 0; dir < tot_dim; dir++) {
      ssxmfcontent << directions[dir];
      if (dir<(tot_dim-1))
        ssxmfcontent << "_";
    }
    ssxmfcontent << "\">" << std::endl;
    for (int dir=0; dir<tot_dim; dir++) {
      ssxmfcontent << "       <DataItem Dimensions=\"";
      for (int tmp_dir = 3-tot_dim; tmp_dir < 3; tmp_dir++) {
        ssxmfcontent << this->nodesize[tmp_dir+1];
        if (tmp_dir<2)
          ssxmfcontent << " ";
      }
      #ifndef XDMF_DOUBLE
      ssxmfcontent << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">" << std::endl;
      #else
      ssxmfcontent << "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
      #endif
      ssxmfcontent << "        " << filename << ":/node_coords/" << directions[dir] << std::endl;
      ssxmfcontent << "       </DataItem>" << std::endl;
    }
    ssxmfcontent << "     </Geometry>" << std::endl;

    std::ofstream xmfFile(filename_xmf, std::ios::trunc);
    xmfFile << ssxmfcontent.str();
    xmfFile.close();
  }
}



/* ********************************************************************* */
void Xdmf::WriteScalar(
                       DUMP_DATATYPE* Vin,
                       const std::string &var_name,
                       const hsize_t *dims,
                       const std::string filename,
                       const std::string filename_xmf,
                       hid_t &memspace,
                       hid_t &dataspace,
                       hid_t &plist_id_mpiio,
                       hid_t &group_fields) {
/*!
* Write HDF5 scalar field.
*
*********************************************************************** */

  std::stringstream ssdataset_name;
  std::string dataset_name;
  std::stringstream ssgroup_name;
  ssgroup_name << "/Timestep_" << xdmfFileNumber << "/vars/";

  dataset_name = var_name.c_str();
  std::string dataset_label = dataset_name.c_str();
  std::transform(dataset_label.begin(), dataset_label.end(), dataset_label.begin(), ::tolower);
  hid_t err, dataset;

  // We define the dataset that contain the fields.

  dataset = H5Dcreate(group_fields, var_name.c_str(), H5_DUMP_DATATYPE,
                        dataspace, H5P_DEFAULT);
  #ifdef WITH_MPI
  err = H5Dwrite(dataset, H5_DUMP_DATATYPE, memspace, dataspace,
                 plist_id_mpiio, Vin);
  #else
  err = H5Dwrite(dataset, H5_DUMP_DATATYPE, memspace, dataspace,
                 H5P_DEFAULT, Vin);
  #endif
  H5Dclose(dataset);

  if (idfx::prank == 0) {
    std::stringstream ssxmfcontent;
    ssxmfcontent << "     <Attribute Name=\"" << dataset_label;
    ssxmfcontent << "\" AttributeType=\"Scalar\" Center=\"Cell\">" << std::endl;

    ssxmfcontent << "       <DataItem Dimensions=\"";
    for (int dir = 0; dir < DIMENSIONS; dir++) {
     ssxmfcontent << dims[dir];
     if (dir<(DIMENSIONS-1))
       ssxmfcontent << " ";
    }

    #ifndef XDMF_DOUBLE
    ssxmfcontent << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">" << std::endl;
    #else
    ssxmfcontent << "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
    #endif
    ssxmfcontent << "        " << filename << ":";
    ssxmfcontent << ssgroup_name.str() << dataset_name << std::endl;
    ssxmfcontent << "       </DataItem>" << std::endl;
    ssxmfcontent << "     </Attribute>" << std::endl;

    std::ofstream xmfFile(filename_xmf, std::ios::app);
    xmfFile << ssxmfcontent.str();
    xmfFile.close();
  }
}

/* ********************************************************************* */
void Xdmf::WriteFooter(
                       const std::string filename,
                       const std::string filename_xmf
                      ) {
/*!
* Write XDMF footer.
*
* \param [in]  filename  hdf5 filename
* \param [in]  filename_xmf  xmf filename
*
*********************************************************************** */
  // We define the dataset.
  std::vector<size_t> dims_cells(3), offset_cells(3), count_cells(3);
  std::vector<std::string> directions {"X", "Y", "Z"};

  for(int dir = 0; dir < 3 ; dir++) {
    dims_cells[dir] = std::size_t(this->cellsize[dir+1]);
    offset_cells[dir] = std::size_t(this->cellstart[dir+1]);
    count_cells[dir] = std::size_t(this->cellsubsize[dir+1]);
  }

  if (idfx::prank==0) {
    std::stringstream ssxmfcontent;
    #if DIMENSIONS == 1
    int tot_dim = 1;
    #elif DIMENSIONS == 2
    int tot_dim = 2;
    #elif DIMENSIONS == 3
    int tot_dim = 3;
    #endif

    for(int dir = 0; dir < tot_dim ; dir++) {
      ssxmfcontent << "     <Attribute Name=\"" << directions[dir];
      ssxmfcontent << "\" AttributeType=\"Scalar\" Center=\"Cell\">" << std::endl;
      ssxmfcontent << "       <DataItem Dimensions=\"";
      for (int tmp_dir = 3-tot_dim; tmp_dir < 3; tmp_dir++) {
        ssxmfcontent << dims_cells[tmp_dir];
        if (tmp_dir < 2)
          ssxmfcontent << " ";
      }
      #ifndef XDMF_DOUBLE
      ssxmfcontent << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">" << std::endl;
      #else
      ssxmfcontent << "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
      #endif
      ssxmfcontent << "        " << filename << ":/cell_coords/" << directions[dir] << std::endl;
      ssxmfcontent << "       </DataItem>" << std::endl;
      ssxmfcontent << "     </Attribute>" << std::endl;
    }
    ssxmfcontent << "   </Grid>" << std::endl;
    ssxmfcontent << " </Domain>" << std::endl;
    ssxmfcontent << "</Xdmf>" << std::endl;

    std::ofstream xmfFile(filename_xmf, std::ios::app);
    xmfFile << ssxmfcontent.str();
    xmfFile.close();
  }
}
