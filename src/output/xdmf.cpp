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
#include "xdmf.hpp"
#include "gitversion.hpp"
#include "idefix.hpp"
#include "dataBlockHost.hpp"
#include "gridHost.hpp"
#include "output.hpp"

// using namespace HighFive;

// Whether or not we write the time in the XDMF file
#define WRITE_TIME

/*init the object */
void Xdmf::Init(Input &input, DataBlock &datain) {
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

  for (int dir=0; dir<3; dir++) {
    this->periodicity[dir] = (datain.mygrid->lbound[dir] == periodic);
  }
  // Create the coordinate array required in XDMF files
  this->nx1 = grid.np_int[IDIR];
  this->nx2 = grid.np_int[JDIR];
  this->nx3 = grid.np_int[KDIR];

  this->nx1loc = data.np_int[IDIR];
  this->nx2loc = data.np_int[JDIR];
  this->nx3loc = data.np_int[KDIR];

  // Temporary storage on host for 3D arrays
  this->vect3D = new DUMP_DATATYPE[nx1loc*nx2loc*nx3loc];

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
  /* -- Data order that is saved is Z-Y-X -- */
  for(int dir = 0; dir < 3 ; dir++) {
    this->nodesize[3-dir] = datain.mygrid->np_int[dir];
    this->nodestart[3-dir] = datain.gbeg[dir]-datain.nghost[dir];
    this->nodesubsize[3-dir] = datain.np_int[dir];

    this->cellsize[3-dir] = datain.mygrid->np_int[dir];
    this->cellstart[3-dir] = datain.gbeg[dir]-datain.nghost[dir];
    this->cellsubsize[3-dir] = datain.np_int[dir];
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

  if(datain.mygrid->xproc[0] == datain.mygrid->nproc[0]-1) this->nodesubsize[3] += IOFFSET;
  if(datain.mygrid->xproc[1] == datain.mygrid->nproc[1]-1) this->nodesubsize[2] += JOFFSET;
  if(datain.mygrid->xproc[2] == datain.mygrid->nproc[2]-1) this->nodesubsize[1] += KOFFSET;

  // Allocate a node and cell views on the host
  node_coord = IdefixHostArray4D<DUMP_DATATYPE>("XdmfNodeCoord", nodesubsize[0],
                                                                 nodesubsize[1],
                                                                 nodesubsize[2],
                                                                 nodesubsize[3]);

  cell_coord = IdefixHostArray4D<DUMP_DATATYPE>("XdmfCellCoord", cellsubsize[0],
                                                                 cellsubsize[1],
                                                                 cellsubsize[2],
                                                                 cellsubsize[3]);
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
        D_EXPAND( x1 = data.xl[IDIR](i + grid.nghost[IDIR] );  ,
                  x2 = data.xl[JDIR](j + grid.nghost[JDIR]);  ,
                  x3 = data.xl[KDIR](k + grid.nghost[KDIR]);  )
        if ( (k<(cellsubsize[1])) && (j<(cellsubsize[2])) && (i<(cellsubsize[3])) ) {
          D_EXPAND( x1_cell = data.x[IDIR](i + grid.nghost[IDIR] );  ,
                    x2_cell = data.x[JDIR](j + grid.nghost[JDIR]);  ,
                    x3_cell = data.x[KDIR](k + grid.nghost[KDIR]);  )
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
  for(int dir = 0; dir < 3 ; dir++) {
    // XDMF assumes Fortran array ordering, hence arrays dimensions are filled backwards
    // So ordering is Z-Y-X
    // offset in the destination array
    this->mpi_data_start[2-dir] = datain.gbeg[dir]-grid.nghost[dir];
    this->mpi_data_size[2-dir] = grid.np_int[dir];
    this->mpi_data_subsize[2-dir] = datain.np_int[dir];
  }
  #endif
}

int Xdmf::Write(DataBlock &datain, Output &output) {
  idfx::pushRegion("Xdmf::Write");
  std::string filename;
  std::string filename_xmf;

  idfx::cout << "Xdmf: Write file n " << xdmfFileNumber << "..." << std::flush;

  timer.reset();

  // Create a copy of the dataBlock on Host, and sync it.
  DataBlockHost data(datain);
  data.SyncFromDevice();

  std::stringstream ssfileName, ssfileNameXmf, ssxdmfFileNum;
  #ifndef XDMF_DOUBLE
  std::string extension = ".flt";
  #else
  std::string extension = ".dbl";
  #endif
  ssxdmfFileNum << std::setfill('0') << std::setw(4) << xdmfFileNumber;
  ssfileName << "data." << ssxdmfFileNum.str() << extension << ".h5";
  ssfileNameXmf << "data." << ssxdmfFileNum.str() << extension << ".xmf";
  filename = ssfileName.str();
  filename_xmf = ssfileNameXmf.str();

  // Open file and write header
#ifdef WITH_MPI
  // MPI-IO requires informing HDF5 that we want something other than
  // the default behaviour. This is done through property lists. We
  // need a file access property list.
  HighFive::FileAccessProps fapl;
  // We tell HDF5 to use MPI-IO
  fapl.add(HighFive::MPIOFileAccess{MPI_COMM_WORLD, MPI_INFO_NULL});
  // We also specify that we want all meta-data related operations
  // to use MPI collective operations. This implies that all MPI ranks
  // in the communicator must participate in any HDF5 operation that
  // reads or writes metadata. Essentially, this is safe if all MPI ranks
  // participate in all HDF5 operations.
  fapl.add(HighFive::MPIOCollectiveMetadata{});

  // Now we can create the file as usual.
  HighFive::File fileHdf(filename, HighFive::File::ReadWrite | HighFive::File::Truncate, fapl);

#else
  // Open a file
  File fileHdf(filename, HighFive::File::ReadWrite | HighFive::File::Truncate);
#endif

  WriteHeader(fileHdf, filename, filename_xmf, datain.t);

  std::stringstream ssgroup_name;
  ssgroup_name << "/Timestep_" << xdmfFileNumber << "/vars";

  // We can create a group as usual.
  HighFive::Group group_fields = fileHdf.createGroup(ssgroup_name.str());

  // Write field one by one
  for(int nv = 0 ; nv < NVAR ; nv++) {
    for(int k = data.beg[KDIR]; k < data.end[KDIR] ; k++ ) {
      for(int j = data.beg[JDIR]; j < data.end[JDIR] ; j++ ) {
        for(int i = data.beg[IDIR]; i < data.end[IDIR] ; i++ ) {
          vect3D[i-data.beg[IDIR] + (j-data.beg[JDIR])*nx1loc + (k-data.beg[KDIR])*nx1loc*nx2loc]
              = static_cast<DUMP_DATATYPE>(data.Vc(nv,k,j,i));
        }
      }
    }
    WriteScalar(vect3D, datain.hydro.VcName[nv], data, filename, filename_xmf, group_fields);
  }
  // Write user-defined variables (when required by output)
  if(output.userDefVariablesEnabled) {
    // Walk the map and make an output for each key of the map
    // (and we thank c++11 for its cute way of doing this)
    for(auto const &variable : output.userDefVariables) {
      for(int k = data.beg[KDIR]; k < data.end[KDIR] ; k++ ) {
        for(int j = data.beg[JDIR]; j < data.end[JDIR] ; j++ ) {
          for(int i = data.beg[IDIR]; i < data.end[IDIR] ; i++ ) {
            vect3D[i-data.beg[IDIR] + (j-data.beg[JDIR])*nx1loc + (k-data.beg[KDIR])*nx1loc*nx2loc]
                = static_cast<DUMP_DATATYPE>(variable.second(k,j,i));
          }
        }
      }
      WriteScalar(vect3D, variable.first, data, filename, filename_xmf, group_fields);
    }
  }

  // Write vector potential if we're using this
  #ifdef EVOLVE_VECTOR_POTENTIAL
  for(int nv = 0 ; nv <= AX3e ; nv++) {
    for(int k = data.beg[KDIR]; k < data.end[KDIR] ; k++ ) {
      for(int j = data.beg[JDIR]; j < data.end[JDIR] ; j++ ) {
        for(int i = data.beg[IDIR]; i < data.end[IDIR] ; i++ ) {
          vect3D[i-data.beg[IDIR] + (j-data.beg[JDIR])*nx1loc + (k-data.beg[KDIR])*nx1loc*nx2loc]
              = static_cast<DUMP_DATATYPE>(data.Ve(nv,k,j,i));
        }
      }
    }
    WriteScalar(vect3D, datain.hydro.VeName[nv], data, filename, filename_xmf, group_fields);
  }
  #endif
  WriteFooter(fileHdf, filename, filename_xmf);

  fileHdf.flush();
  xdmfFileNumber++;
  // Make file number
  idfx::cout << "done in " << timer.seconds() << " s." << std::endl;

  idfx::popRegion();
  // One day, we will have a return code.
  return(0);
}


/* ********************************************************************* */
void Xdmf::WriteHeader(
                       HighFive::File fileHdf,
                       const std::string filename,
                       const std::string filename_xmf,
                       real time) {
/*!
* Write XDMF header in parallel or serial mode.
*
* \param [in]  fileHdf  pointer to file
* \param [in]  filename  hdf5 filename
* \param [in]  filename_xmf  xmf filename
* \param [in]  time  simulation time in code units
*
*********************************************************************** */

  std::string header;
  std::stringstream ssheader;
  //std::stringstream ssgroup_name;

  /* -------------------------------------------
  1. Header
  ------------------------------------------- */
  // Create a dummy dataset of one single integer that will hold some metadata as attributes
  auto info_dataset = fileHdf.createDataSet("info",
                                            HighFive::DataSpace(1),
                                            HighFive::create_datatype<int>());
  ssheader << "Idefix " << GITVERSION << " XDMF Data";
  info_dataset.createAttribute<std::string>("version", ssheader.str());

  /* ------------------------------------------
  2. File format
  ------------------------------------------ */

  #ifndef XDMF_DOUBLE
  info_dataset.createAttribute<std::string>("dtype", "float");
  #else
  info_dataset.createAttribute<std::string>("dtype", "double");
  #endif

  #ifdef WRITE_TIME
  info_dataset.createAttribute("time", time);
  #endif

  // We can create groups as usual.
  HighFive::Group group_cells  = fileHdf.createGroup("cell_coords");
  HighFive::Group group_nodes  = fileHdf.createGroup("node_coords");

  // We define the dataset.
  std::vector<size_t> dims_nodes(3), offset_nodes(3), count_nodes(3);
  std::vector<size_t> dims_cells(3), offset_cells(3), count_cells(3);
  std::vector<std::string> directions{"X", "Y", "Z"};

  for(int dir = 0; dir < 3 ; dir++) {
    dims_nodes[dir] = std::size_t(this->nodesize[dir+1]);
    offset_nodes[dir] = std::size_t(this->nodestart[dir+1]);
    count_nodes[dir] = std::size_t(this->nodesubsize[dir+1]);

    dims_cells[dir] = std::size_t(this->cellsize[dir+1]);
    offset_cells[dir] = std::size_t(this->cellstart[dir+1]);
    count_cells[dir] = std::size_t(this->cellsubsize[dir+1]);
  }

  // We create the dataset for the nodes and cells
  HighFive::DataSet *dataset_nodes = static_cast<HighFive::DataSet *>
                                     ( malloc(3*sizeof(HighFive::DataSet)) );
  HighFive::DataSet *dataset_cells = static_cast<HighFive::DataSet *>
                                     ( malloc(3*sizeof(HighFive::DataSet)) );

  for(int dir = 0; dir < 3 ; dir++) {
    dataset_nodes[dir] = group_nodes.createDataSet<DUMP_DATATYPE>
                                     (directions[dir], HighFive::DataSpace(dims_nodes));
    dataset_cells[dir] = group_cells.createDataSet<DUMP_DATATYPE>
                                     (directions[dir], HighFive::DataSpace(dims_cells));
  }

#ifdef WITH_MPI
  auto xfer_props = HighFive::DataTransferProps{};
  xfer_props.add(HighFive::UseCollectiveIO{});

  for(int dir = 0; dir < 3 ; dir++) {
    dataset_nodes[dir].select(offset_nodes, count_nodes).write_raw(
                       Kokkos::subview (this->node_coord,
                                        dir,
                                        Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL()).data(),
                                        xfer_props);
    dataset_cells[dir].select(offset_cells, count_cells).write_raw(
                       Kokkos::subview (this->cell_coord,
                                        dir,
                                        Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL()).data(),
                                        xfer_props);
  }
#else
  for(int dir = 0; dir < 3 ; dir++) {
    dataset_nodes[dir].write_raw(
                       Kokkos::subview (this->node_coord,
                                        dir,
                                        Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL()).data() );
    dataset_cells[dir].write_raw(
                       Kokkos::subview (this->cell_coord,
                                        dir,
                                        Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL()).data() );
  }
#endif
  if (idfx::prank==0) {
    std::stringstream ssxmfcontent;
    /*      XMF file generation      */
    /* -------------------------------------------
    1. File version and identifier
    ------------------------------------------- */

    ssxmfcontent << "<?xml version=\"1.0\" ?>" << std::endl;
    ssxmfcontent << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>" << std::endl;
    ssxmfcontent << "<Xdmf Version=\"2.0\">" << std::endl;
    ssxmfcontent << " <Domain>" << std::endl;
    ssxmfcontent << "  <Grid Name=\"node_mesh\" GridType=\"Uniform\">" << std::endl;
    ssxmfcontent << "    <Time Value=\"" << time << "\"/>" << std::endl;
    #if DIMENSIONS == 1
    ssxmfcontent << "     <Topology TopologyType=\"1DSMesh\" NumberOfElements=\"";
    ssxmfcontent << this->nodesize[1] << "\"/>" << std::endl;
    ssxmfcontent << "     <Geometry GeometryType=\"X\">" << std::endl;
    ssxmfcontent << "       <DataItem Dimensions=\"" << this->nodesize[1];
    #ifndef XDMF_DOUBLE
    ssxmfcontent << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">" << std::endl;
    #else
    ssxmfcontent << "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
    #endif
    ssxmfcontent << "        ./" << filename << ":/node_coords/" << directions[0] << std::endl;
    ssxmfcontent << "       </DataItem>" << std::endl;
    #elif DIMENSIONS == 2
    ssxmfcontent << "     <Topology TopologyType=\"2DSMesh\" NumberOfElements=\"";
    for (int dir = 0; dir < 2; dir++) {
      ssxmfcontent << this->nodesize[dir+1];
      if (dir<1)
        ssxmfcontent << " ";
    }
    ssxmfcontent << "\"/>" << std::endl;
    ssxmfcontent << "     <Geometry GeometryType=\"X_Y\">" << std::endl;
    for (int dir = 0; dir < 2; dir++) {
      ssxmfcontent << "       <DataItem Dimensions=\"";
      for (int tmp_dir = 0; tmp_dir < 2; tmp_dir++) {
        ssxmfcontent << this->nodesize[tmp_dir+1];
        if (tmp_dir<1)
          ssxmfcontent << " ";
      }
      #ifndef XDMF_DOUBLE
      ssxmfcontent << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">" << std::endl;
      #else
      ssxmfcontent << "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
      #endif
      ssxmfcontent << "        ./" << filename << ":/node_coords/" << directions[dir] << std::endl;
      ssxmfcontent << "       </DataItem>" << std::endl;
    }
    #elif DIMENSIONS == 3
    ssxmfcontent << "     <Topology TopologyType=\"3DSMesh\" NumberOfElements=\"";
    for (int dir = 0; dir < 3; dir++) {
      ssxmfcontent << this->nodesize[dir+1];
      if (dir<2)
        ssxmfcontent << " ";
    }
    ssxmfcontent << "\"/>" << std::endl;
    ssxmfcontent << "     <Geometry GeometryType=\"X_Y_Z\">" << std::endl;
    for (int dir = 0; dir < 3; dir++) {
      ssxmfcontent << "       <DataItem Dimensions=\"";
      for (int tmp_dir = 0; tmp_dir < 3; tmp_dir++) {
        ssxmfcontent << this->nodesize[tmp_dir+1];
        if (tmp_dir<2)
          ssxmfcontent << " ";
      }
      #ifndef XDMF_DOUBLE
      ssxmfcontent << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">" << std::endl;
      #else
      ssxmfcontent << "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
      #endif
      ssxmfcontent << "        ./" << filename << ":/node_coords/" << directions[dir] << std::endl;
      ssxmfcontent << "       </DataItem>" << std::endl;
    }
    #endif
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
                       DataBlockHost data,
                       const std::string filename,
                       const std::string filename_xmf,
                       HighFive::Group group_fields) {
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

  // We define the dataset.
  std::vector<size_t> dims(3);
#ifdef WITH_MPI
  for(int dir = 0; dir < 3 ; dir++)
    dims[dir] = std::size_t(this->mpi_data_size[dir]);
#else
  dims[0] = std::size_t(this->nx3loc);
  dims[1] = std::size_t(this->nx2loc);
  dims[2] = std::size_t(this->nx1loc);
#endif

  // We create the dataset for the scalar
  HighFive::DataSet dataset = group_fields.createDataSet<DUMP_DATATYPE>
                              ( dataset_name, HighFive::DataSpace(dims) );

#ifdef WITH_MPI
  auto xfer_props = HighFive::DataTransferProps{};
  xfer_props.add(HighFive::UseCollectiveIO{});

  // Each MPI rank writes a non-overlapping part of the array.
  std::vector<std::size_t> offset(3), count(3);

  for(int dir = 0; dir < 3 ; dir++) {
    offset[dir] = std::size_t(this->mpi_data_start[dir]);
    count[dir]  = std::size_t(this->mpi_data_subsize[dir]);
  }
  dataset.select(offset, count).write_raw(Vin, xfer_props);
#else
  dataset.write_raw(Vin);
#endif
  if (idfx::prank == 0) {
    std::stringstream ssxmfcontent;
    ssxmfcontent << "     <Attribute Name=\"" << dataset_label;
    ssxmfcontent << "\" AttributeType=\"Scalar\" Center=\"Cell\">" << std::endl;
    #if DIMENSIONS == 1
    #ifndef XDMF_DOUBLE
    ssxmfcontent << "       <DataItem Dimensions=\"" << dims[0];
    ssxmfcontent << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">" << std::endl;
    #else
    ssxmfcontent << "       <DataItem Dimensions=\"" << dims[0];
    ssxmfcontent << "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
    #endif
    #elif DIMENSIONS == 2
    ssxmfcontent << "       <DataItem Dimensions=\"";
    for (int dir = 0; dir < 2; dir++) {
     ssxmfcontent << dims[dir];
     if (dir<1)
       ssxmfcontent << " ";
    }
    #ifndef XDMF_DOUBLE
    ssxmfcontent << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">" << std::endl;
    #else
    ssxmfcontent << "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
    #endif
    #elif DIMENSIONS == 3
    ssxmfcontent << "       <DataItem Dimensions=\"";
    for (int dir = 0; dir < 3; dir++) {
     ssxmfcontent << dims[dir];
     if (dir<2)
       ssxmfcontent << " ";
    }
    #ifndef XDMF_DOUBLE
    ssxmfcontent << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">" << std::endl;
    #else
    ssxmfcontent << "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
    #endif
    #endif
    ssxmfcontent << "        ./" << filename << ":";
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
                       HighFive::File fileHdf,
                       const std::string filename,
                       const std::string filename_xmf
                      ) {
/*!
* Write XDMF header in parallel or serial mode.
*
* \param [in]  fileHdf  pointer to file
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

    for(int dir = 0; dir < 3 ; dir++) {
      ssxmfcontent << "     <Attribute Name=\"" << directions[dir];
      ssxmfcontent << "\" AttributeType=\"Scalar\" Center=\"Cell\">" << std::endl;
      ssxmfcontent << "       <DataItem Dimensions=\"";
      for (int tmp_dir = 0; tmp_dir < 3; tmp_dir++) {
        ssxmfcontent << dims_cells[dir];
        if (tmp_dir < (tot_dim-1))
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
