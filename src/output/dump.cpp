// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include <algorithm>
#include <unordered_set>
#if __has_include(<filesystem>)
  #include <filesystem>
  namespace fs = std::filesystem;
#elif __has_include(<experimental/filesystem>)
  #include <experimental/filesystem>
  namespace fs = std::experimental::filesystem;
#else
  #error "Missing the <filesystem> header."
#endif
#include <iomanip>
#include "dump.hpp"
#include "version.hpp"
#include "dataBlockHost.hpp"
#include "gridHost.hpp"
#include "output.hpp"
#include "fluid.hpp"

// Max size of array name
#define  NAMESIZE     16
#define  FILENAMESIZE   256
#define  HEADERSIZE 128

// Register a variable to be dumped (and read)

void Dump::RegisterVariable(IdefixArray3D<real>& in,
                        std::string name,
                        int dir,
                        DumpField::ArrayLocation loc) {
    dumpFieldMap.emplace(name, DumpField(in, loc, dir));
}

void  Dump::RegisterVariable(IdefixHostArray3D<real>& in,
                        std::string name,
                        int dir,
                        DumpField::ArrayLocation loc) {
    dumpFieldMap.emplace(name, DumpField(in, loc, dir));
}

void  Dump::RegisterVariable(IdefixArray4D<real>& in,
                        std::string name,
                        int varnum,
                        int dir,
                        DumpField::ArrayLocation loc) {
    dumpFieldMap.emplace(name, DumpField(in, varnum, loc, dir));
}

void  Dump::RegisterVariable(IdefixHostArray4D<real>& in,
                        std::string name,
                        int varnum,
                        int dir,
                        DumpField::ArrayLocation loc) {
    dumpFieldMap.emplace(name, DumpField(in, varnum, loc, dir));
}



void Dump::CreateMPIDataType(GridBox gb, bool read) {
  #ifdef WITH_MPI
    int start[3];
    int size[3];
    int subsize[3];

    // the grid is required to now the current MPÏ domain decomposition
    Grid *grid = data->mygrid;

    // Dimensions for cell-centered fields
    for(int dir = 0; dir < 3 ; dir++) {
      size[2-dir] = gb.sizeGlob[dir];
      start[2-dir] = gb.start[dir];
      subsize[2-dir] = gb.size[dir];
    }
    if(read) {
      MPI_SAFE_CALL(MPI_Type_create_subarray(3, size, subsize, start,
                                            MPI_ORDER_C, realMPI, &this->descCR));
      MPI_SAFE_CALL(MPI_Type_commit(&this->descCR));
    } else {
      MPI_SAFE_CALL(MPI_Type_create_subarray(3, size, subsize, start,
                                            MPI_ORDER_C, realMPI, &this->descCW));
      MPI_SAFE_CALL(MPI_Type_commit(&this->descCW));
    }

    // Dimensions for face-centered field
    for(int face = 0; face < 3 ; face++) {
      for(int dir = 0; dir < 3 ; dir++) {
        size[2-dir] = gb.sizeGlob[dir];
        start[2-dir] = gb.start[dir];
        subsize[2-dir] = gb.size[dir];
      }
      if(read) {
        // Add the extra guy in the face direction
        size[2-face]++;
        subsize[2-face]++; // valid only for reading
                          //since it involves an overlap of data between procs

        MPI_SAFE_CALL(MPI_Type_create_subarray(3, size, subsize, start,
                                              MPI_ORDER_C, realMPI, &this->descSR[face]));
        MPI_SAFE_CALL(MPI_Type_commit(&this->descSR[face]));
      } else {
        // Now for writing, it is only the last proc which keeps one additional cell
        size[2-face]++;
        if(grid->xproc[face] == grid->nproc[face] - 1  ) subsize[2-face]++;
        MPI_SAFE_CALL(MPI_Type_create_subarray(3, size, subsize, start,
                                              MPI_ORDER_C, realMPI, &this->descSW[face]));
        MPI_SAFE_CALL(MPI_Type_commit(&this->descSW[face]));
      }
    }
    // Dimensions for edge-centered field
    for(int nv = 0; nv < 3 ; nv++) {
      // load the array size
      for(int dir = 0; dir < 3 ; dir++) {
        size[2-dir] = gb.sizeGlob[dir];
        start[2-dir] = gb.start[dir];
        subsize[2-dir] = gb.size[dir];
      }

      if(read) {
        // Extra cell in the dirs perp to field
        for(int i = 0 ; i < DIMENSIONS ; i++) {
          if(i!=nv) {
            size[2-i]++;
            subsize[2-i]++; // valid only for reading
                            //since it involves an overlap of data between procs
          }
        }
        MPI_SAFE_CALL(MPI_Type_create_subarray(3, size, subsize, start,
                                              MPI_ORDER_C, realMPI, &this->descER[nv]));
        MPI_SAFE_CALL(MPI_Type_commit(&this->descER[nv]));
      } else {
        // Now for writing, it is only the last proc which keeps one additional cell,
        // so we remove what we added for reads
        for(int i = 0 ; i < DIMENSIONS ; i++) {
          if(i!=nv) {
            size[2-i]++;
            if(grid->xproc[i] == grid->nproc[i] - 1  ) {
              subsize[2-i]++;
            }
          }
        }
        MPI_SAFE_CALL(MPI_Type_create_subarray(3, size, subsize, start,
                                              MPI_ORDER_C, realMPI, &this->descEW[nv]));
        MPI_SAFE_CALL(MPI_Type_commit(&this->descEW[nv]));
      }
    }
  #endif
}

void Dump::Init(DataBlock *datain) {
  idfx::pushRegion("Dump::Init");
  this->data = datain;

  for (int dir=0; dir<3; dir++) {
    this->periodicity[dir] = (data->mygrid->lbound[dir] == periodic);
  }
  this->dumpFileNumber = 0;

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

  // Allocate scratch Array
  // It must be able to handle either a full domain 1D Array and
  // a sub-3D somain
  int64_t nmax = (data->np_int[IDIR]+IOFFSET)*
                  (data->np_int[JDIR]+JOFFSET)*
                  (data->np_int[KDIR]+KOFFSET);
  nmax = std::max(nmax,static_cast<int64_t>(data->mygrid->np_tot[IDIR]));
  nmax = std::max(nmax,static_cast<int64_t>(data->mygrid->np_tot[JDIR]));
  nmax = std::max(nmax,static_cast<int64_t>(data->mygrid->np_tot[KDIR]));

  this->scrch = new real[nmax];

  #ifdef WITH_MPI
    Grid *grid = data->mygrid;
    GridBox gb;
    for(int dir = 0; dir < 3 ; dir++) {
      gb.start[dir] = data->gbeg[dir]-data->nghost[dir];
      gb.size[dir] = data->np_int[dir];
      gb.sizeGlob[dir] = grid->np_int[dir];
    }
    // Create MPI datatypes for read/write
    CreateMPIDataType(gb, false);
    CreateMPIDataType(gb, true);
  #endif

  // Register variables that are needed in restart dumps
  this->RegisterVariable(&dumpFileNumber, "dumpFileNumber");
  this->RegisterVariable(&geometry, "geometry");
  this->RegisterVariable(periodicity, "periodicity", 3);

  idfx::popRegion();
}


Dump::Dump(Input &input, DataBlock *datain) {
  // Constructor with an input object, in which case,
  // We use the outputdirectory provided in the input

    // initialize output path
  if(input.CheckEntry("Output","dmp_dir")>=0) {
    outputDirectory = input.Get<std::string>("Output","dmp_dir",0);
  } else {
    outputDirectory = "./";
  }
  Init(datain);
}

Dump::Dump(DataBlock *datain) {
  // Constructor without an input object, in which case,
  // We use the default output directory

  outputDirectory = "./";
  Init(datain);
}

Dump::~Dump() {
  delete scrch;
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
  if(type == BoolType) size=sizeof(bool);

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

void Dump::Skip(IdfxFileHandler fileHdl, int ndim, int *dim,
                            DataType type) {
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
    offset+= ntot*size;
  #else
    fseek(fileHdl, ntot*size, SEEK_CUR);
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

// Helper function to convert filesystem::file_time into std::time_t
// see https://stackoverflow.com/questions/56788745/
// This conversion "hack" is required in C++17 as no proper conversion bewteen
// fs::last_write_time and std::time_t
// exists in the standard library until C++20
template <typename TP>
std::time_t to_time_t(TP tp) {
    auto sctp = std::chrono::time_point_cast<std::chrono::system_clock::duration>
                (tp - TP::clock::now() + std::chrono::system_clock::now());
    return std::chrono::system_clock::to_time_t(sctp);
}

int Dump::GetLastDumpInDirectory(fs::path &directory) {
  int num = -1;

  std::time_t youngFileTime;
  bool first = true;
  for (const auto & entry : fs::directory_iterator(directory)) {
      // Check file extension
      if(entry.path().extension().string().compare(".dmp")==0) {
        auto fileTime = to_time_t(fs::last_write_time(entry.path()));
        // Check which one is the most recent
        if(first || fileTime>youngFileTime) {
          // std::tm *gmt = std::gmtime(&fileTime);
          // idfx::cout << "file " << entry.path() << "is the most recent with "
          //            << std::put_time(gmt, "%d %B %Y %H:%M:%S") << std::endl;

          // Ours is more recent, extract the dump file number
          try {
            num = std::stoi(entry.path().filename().string().substr(5,4));
            first = false;
            youngFileTime = fileTime;
          } catch (...) {
            // the file name does not follow the convention "filebase.xxxx.dmp"
            // ->We skip it
          }
        }
      }
  }
  return(num);
}
bool Dump::Read(Output& output, int readNumber ) {
  fs::path filename;
  int nx[3];
  int nxglob[3];
  std::string fieldName;
  std::string eof ("eof");
  DataType type;
  int ndim;
  IdfxFileHandler fileHdl;

  idfx::pushRegion("Dump::Read");

  fs::path readDir = this->outputDirectory;

  if(readNumber<0) {
    // We actually don't know which file we're supposed to read, so let's guess
    readNumber = GetLastDumpInDirectory(readDir);
    if(readNumber < 0) {
      idfx::cout << "Dump: cannot find a valid dump in " << this->outputDirectory << std::endl;
      if(outputDirectory.compare("./")!=0) {
        idfx::cout << "Dump: reverting to the current directory." << std::endl;
        readDir = ".";
        readNumber = GetLastDumpInDirectory(readDir);
      }
      if(readNumber<0) {
        IDEFIX_WARNING("cannot find a valid restart dump.");
        return(false);
      }
    }
  }

  // Reset timer
  timer.reset();

  // Set filename
  std::stringstream ssdumpFileNum,ssFileName;
  ssdumpFileNum << std::setfill('0') << std::setw(4) << readNumber;
  ssFileName << "dump." << ssdumpFileNum.str() << ".dmp";
  filename = readDir/ssFileName.str();

  idfx::cout << "Dump: Reading " << filename << "..." << std::flush;
  // open file
#ifdef WITH_MPI
  MPI_SAFE_CALL(MPI_File_open(MPI_COMM_WORLD, filename.c_str(),
                              MPI_MODE_RDONLY | MPI_MODE_UNIQUE_OPEN,
                              MPI_INFO_NULL, &fileHdl));
  this->offset = 0;
#else
  fileHdl = fopen(filename.c_str(),"rb");
  if(fileHdl == NULL) {
    std::stringstream msg;
    msg << "Failed to open dump file: " << std::string(filename) << std::endl;
    IDEFIX_ERROR(msg);
  }
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
    if(nx[0] != data->mygrid->np_int[dir]) {
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

  std::unordered_set<std::string> notFound {};
  for(auto it = dumpFieldMap.begin(); it != dumpFieldMap.end(); it++) {
    notFound.insert(it->first);
  }

  // Coordinates are ok, load the bulk
  while(true) {
    ReadNextFieldProperties(fileHdl, ndim, nxglob, type, fieldName);

    /*idfx::cout << "Next field is " << fieldName << " with " << ndim << " dimensions and (";
    for(int i = 0 ; i < ndim ; i++) idfx::cout << nxglob[i] << " ";
    idfx::cout << ") points." << std::endl;*/

    if(fieldName.compare(eof) == 0) {
      // We have reached end of dump file
      break;
    } else {
      if(auto it = dumpFieldMap.find(fieldName) ; it != dumpFieldMap.end()) {
        // This key has been registered
        notFound.erase(fieldName);
        DumpField &scalar = it->second;
        if(scalar.GetType() == DumpField::Type::IdefixArray) {
          // Distributed idefix array
          int direction = scalar.GetDirection();

          // Load it
          for(int dir = 0 ; dir < 3; dir++) {
            nx[dir] = data->np_int[dir];
          }

          if(scalar.GetLocation() == DumpField::ArrayLocation::Face) {
            nx[direction]++;   // Extra cell in the dir direction for face-centered fields
          }
          if(scalar.GetLocation() == DumpField::ArrayLocation::Edge) {
            // Extra cell in the dirs perp to field
            for(int i = 0 ; i < DIMENSIONS ; i++) {
              if(i!=direction) nx[i] ++;
            }
          }
          if(scalar.GetLocation() == DumpField::ArrayLocation::Center) {
            ReadDistributed(fileHdl, ndim, nx, nxglob, descCR, scrch);
          } else if(scalar.GetLocation() == DumpField::ArrayLocation::Face) {
            ReadDistributed(fileHdl, ndim, nx, nxglob, descSR[direction], scrch);
          } else if(scalar.GetLocation() == DumpField::ArrayLocation::Edge) {
            ReadDistributed(fileHdl, ndim, nx, nxglob, descER[direction], scrch);
          }
          auto toRead = scalar.GetHostField<IdefixHostArray3D<real>>();
          // Load the scratch space in designated field
          for(int k = 0; k < nx[KDIR]; k++) {
            for(int j = 0 ; j < nx[JDIR]; j++) {
              for(int i = 0; i < nx[IDIR]; i++) {
                toRead(k+data->beg[KDIR],j+data->beg[JDIR],i+data->beg[IDIR]) =
                                                        scrch[i + j*nx[IDIR] + k*nx[IDIR]*nx[JDIR]];
              }
            }
          }
          scalar.SyncFrom(toRead);
        } else {
          // Fundamental Type
          // Check that size matches
          if(nxglob[0] != scalar.GetSize()) {
            idfx::cout << "nxglob=" << nxglob[0] << " scalar=" << scalar.GetSize() << std::endl;
            IDEFIX_ERROR("Size of field "+fieldName+" do not match");
          }
          // Todo: check that type matches
          void *ptr = scalar.GetHostField<void *>();

          ReadSerial(fileHdl, ndim, nxglob, type, ptr);
        }
      } else {
        Skip(fileHdl, ndim, nxglob, type);
        // Key has not been registered, throw a warning
        IDEFIX_WARNING("Cannot find a field matching " + fieldName
                       + " in current running code. Skipping.");
      }
    }
  }
  if (notFound.size() > 0) {
    std::stringstream msg {};
    msg << "The following fields were not found in " << filename << ": ";
    for(auto it = notFound.begin(); it != notFound.end(); it++) {
      msg << *it << ' ';
    }
    IDEFIX_WARNING(msg);
  }
  #ifdef WITH_MPI
  MPI_SAFE_CALL(MPI_File_close(&fileHdl));
  #else
  fclose(fileHdl);
  #endif

  idfx::cout << "done in " << timer.seconds() << " s." << std::endl;
  idfx::cout << "Restarting from t=" << data->t << "." << std::endl;

  idfx::popRegion();

  return(true);
}


int Dump::Write(Output& output) {
  fs::path filename;
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


  // Set filenames
  std::stringstream ssdumpFileNum,ssFileName;
  ssdumpFileNum << std::setfill('0') << std::setw(4) << dumpFileNumber;
  ssFileName << "dump." << ssdumpFileNum.str() << ".dmp";
  filename = outputDirectory/ssFileName.str();

  dumpFileNumber++;   // For next one

  // Check if file exists, if yes, delete it
  if(idfx::prank==0) {
    if(fs::exists(filename)) {
      fs::remove(filename);
    }
  }

  // open file
#ifdef WITH_MPI
  MPI_Barrier(MPI_COMM_WORLD);
  // Open file for creating, return error if file already exists.
  MPI_SAFE_CALL(MPI_File_open(MPI_COMM_WORLD, filename.c_str(),
                              MPI_MODE_CREATE | MPI_MODE_RDWR
                              | MPI_MODE_EXCL | MPI_MODE_UNIQUE_OPEN,
                              MPI_INFO_NULL, &fileHdl));
  this->offset = 0;
#else
  fileHdl = fopen(filename.c_str(),"wb");
#endif
  // File is open
  // First thing we need are coordinates: init a host mirror and sync it
  GridHost gridHost(*data->mygrid);
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
  std::snprintf(header, HEADERSIZE, "Idefix %s Dump Data %s endian",
                IDEFIX_VERSION, endian.c_str());
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

  for(auto const& [name, scalar] : dumpFieldMap) {
    // Todo: replace these C char by std::string
    std::snprintf(fieldName,NAMESIZE,"%s",name.c_str());
    if(scalar.GetType() == DumpField::Type::IdefixArray) {
      auto toWrite = scalar.GetHostField<IdefixHostArray3D<real>>();
      int dir = scalar.GetDirection();
      for(int i = 0; i < 3 ; i++) {
        nx[i] = data->np_int[i];
        nxtot[i] = gridHost.np_int[i];
      }

      if(scalar.GetLocation() == DumpField::ArrayLocation::Face) {
        // If it is the last datablock of the dimension, increase the size by one to get the last
        //active face of the staggered mesh.
        if(data->mygrid->xproc[dir] == data->mygrid->nproc[dir] - 1  ) nx[dir]++;
        nxtot[dir]++;
      }

      if(scalar.GetLocation() == DumpField::ArrayLocation::Edge) {
        // If it is the last datablock of the dimension, increase the size by one in the direction
        // perpendicular to the vector.
        for(int i = 0 ; i < DIMENSIONS ; i++) {
          if(i != dir) {
            if(data->mygrid->xproc[i] == data->mygrid->nproc[i] - 1) nx[i]++;
            nxtot[i]++;
          }
        }
      }

      // Load the dataset in the scratch array
      for(int k = 0; k < nx[KDIR]; k++) {
        for(int j = 0 ; j < nx[JDIR]; j++) {
          for(int i = 0; i < nx[IDIR]; i++) {
            scrch[i + j*nx[IDIR] + k*nx[IDIR]*nx[JDIR]] = toWrite(k+data->beg[KDIR],
                                                                  j+data->beg[JDIR],
                                                                  i+data->beg[IDIR]);
          }
        }
      }

      if(scalar.GetLocation() == DumpField::ArrayLocation::Center) {
        WriteDistributed(fileHdl, 3, nx, nxtot, fieldName, this->descCW, scrch);
      } else if(scalar.GetLocation() == DumpField::ArrayLocation::Face) {
        WriteDistributed(fileHdl, 3, nx, nxtot, fieldName, this->descSW[dir], scrch);
      } else if(scalar.GetLocation() == DumpField::ArrayLocation::Edge) {
         WriteDistributed(fileHdl, 3, nx, nxtot, fieldName, this->descEW[dir], scrch);
      } else {
        IDEFIX_ERROR("Unknown scalar type for dump write");
      }

    } else {
      // Scalar type if a fundamental type, not distributed
      DataType thisType;
      if(scalar.GetType()==DumpField::Type::Int) thisType = DataType::IntegerType;
      if(scalar.GetType()==DumpField::Type::Single) thisType = DataType::SingleType;
      if(scalar.GetType()==DumpField::Type::Double) thisType = DataType::DoubleType;
      if(scalar.GetType()==DumpField::Type::Bool) thisType = DataType::BoolType;

      nx[0] = scalar.GetSize();

      WriteSerial(fileHdl, 1, nx, thisType, fieldName, scalar.GetHostField<void*>());
    }
  }

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
