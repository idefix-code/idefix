// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#define PYBIND11_DETAILED_ERROR_MESSAGES

#include "pydefix.hpp"
#include <pybind11/embed.h> // everything needed for embedding
#include <pybind11/numpy.h> // for numpy arrays
#include <pybind11/stl.h>   // For STL vectors and containers
#include <string>
#include <vector>
#include "idefix.hpp"
#include "dataBlock.hpp"
#include "dataBlockHost.hpp"


namespace py = pybind11;

int Pydefix::ninstance = 0;

/************************************
 * DataBlockHost Python binding
 * **********************************/

IdefixHostArray4D<real> GatherIdefixArray(IdefixHostArray4D<real> in, DataBlockHost dataHost) {
  idfx::cout << "I am being called" << std::endl;
  idfx::cout << "My dims=" << in.extent(0) << " " << in.extent(1) << " " << in.extent(2) << " " << in.extent(3) << std::endl;
  idfx::cout << "First element=" << in(0,0,0,0) << std::endl; 
  Grid *grid = dataHost.data->mygrid;
  IdefixHostArray4D<real> out;
  #ifdef WITH_MPI
    if(idfx::psize > 1) {
      const int nvar = in.extent(0);
      out = IdefixHostArray4D<real>("pydefix::GatheredArray",nvar,grid->np_tot[KDIR],grid->np_tot[JDIR],grid->np_tot[IDIR]);
      if(idfx::prank == 0) {
        for(int rank = 0 ; rank < idfx::psize ; rank++) {
          // np_tot: total size of the incoming array
          // np_int: size that should be copied into global
          // beg: offset in the incoming array where copy should begin
          // gbeg: offset in the global array where copy should be begin
          MPI_Status status;
          std::array<int,3> np_int,np_tot, beg, gbeg;
          IdefixHostArray4D<real> buf;

          if(rank==0) {
            np_int = dataHost.np_int;
            np_tot = dataHost.np_tot;
            gbeg = dataHost.gbeg;
            beg = dataHost.beg;

            // Add back boundaries
            for(int dir = 0 ; dir < DIMENSIONS ; dir++) {
              if(dataHost.lbound[dir] != internal) {
                np_int[dir] += dataHost.nghost[dir];
                gbeg[dir] -= dataHost.nghost[dir];
                beg[dir] -= dataHost.nghost[dir];
              }
              if(dataHost.rbound[dir] != internal) {
                np_int[dir] += dataHost.nghost[dir];
              }
            }
            buf = in;
          } else { // target rank >0
            // Fetch the local array size
            MPI_Recv(np_int.data(), 3, MPI_INT, rank, 010, MPI_COMM_WORLD, &status);
            MPI_Recv(np_tot.data(), 3, MPI_INT, rank, 011, MPI_COMM_WORLD, &status);
            MPI_Recv(beg.data(), 3, MPI_INT, rank, 012, MPI_COMM_WORLD, &status);
            MPI_Recv(gbeg.data(), 3, MPI_INT, rank, 013, MPI_COMM_WORLD, &status);

            buf = IdefixHostArray4D<real>("pydefix::tempArray",nvar,np_tot[KDIR],np_tot[JDIR],np_tot[IDIR]);
            // Fetch data
            MPI_Recv(buf.data(), nvar*np_tot[IDIR]*np_tot[JDIR]*np_tot[KDIR], realMPI, rank, 014, MPI_COMM_WORLD,&status);
          } // target rank
          // Copy data from the buffer
          for(int n = 0 ; n < nvar ; n++) {
            for(int k = 0 ; k < np_int[KDIR] ; k++) {
              const int kt = k+gbeg[KDIR];
              for(int j = 0 ; j < np_int[JDIR] ; j++) {
                const int jt = j+gbeg[JDIR];
                for(int i = 0 ; i < np_int[IDIR] ; i++) {
                  const int it = i+gbeg[IDIR];
                  out(n,kt,jt,it) = buf(n,k+beg[KDIR],j+beg[JDIR],i+beg[IDIR]);
                }
              }
            }
          }// End for
        }// End loop on target rank for root process
      } else { // MPI prank >0
        std::array<int,3> np_int = dataHost.np_int;
        std::array<int,3> np_tot = dataHost.np_tot;
        std::array<int,3> gbeg = dataHost.gbeg;
        std::array<int,3> beg = dataHost.beg;

        // Add back boundaries
        for(int dir = 0 ; dir < DIMENSIONS ; dir++) {
          if(dataHost.lbound[dir] != internal) {
            np_int[dir] += dataHost.nghost[dir];
            gbeg[dir] -= dataHost.nghost[dir];
            beg[dir] -= dataHost.nghost[dir];
          }
          if(dataHost.rbound[dir] != internal) {
            np_int[dir] += dataHost.nghost[dir];
          }
        }

        // send the local array size
        MPI_Send(np_int.data(), 3, MPI_INT, 0, 010, MPI_COMM_WORLD);
        MPI_Send(np_tot.data(), 3, MPI_INT, 0, 011, MPI_COMM_WORLD);
        MPI_Send(beg.data(), 3, MPI_INT, 0, 012, MPI_COMM_WORLD);
        MPI_Send(gbeg.data(), 3, MPI_INT, 0, 013, MPI_COMM_WORLD);
        MPI_Send(in.data(), nvar*np_tot[IDIR]*np_tot[JDIR]*np_tot[KDIR], realMPI, 0, 014, MPI_COMM_WORLD);        // Allocate array
      }
      // All is transfered
      MPI_Bcast(out.data(), nvar*grid->np_tot[KDIR]*grid->np_tot[JDIR]*grid->np_tot[IDIR], realMPI, 0, MPI_COMM_WORLD);
    } else { // MPI with a single proc
      out = in;
    }
  #else // No MPI
    out = in;
  #endif
  idfx::cout << "final dims=" << out.extent(0) << " " << out.extent(1) << " " << out.extent(2) << " " << out.extent(3) << std::endl;
  return out;
} 

PYBIND11_EMBEDDED_MODULE(pydefix, m) {
  py::class_<DataBlockHost>(m, "DataBlockHost")
    .def(py::init<>())
    .def_readwrite("x", &DataBlockHost::x)
    .def_readwrite("xr", &DataBlockHost::xr)
    .def_readwrite("xl", &DataBlockHost::xl)
    .def_readwrite("dx", &DataBlockHost::dx)
    .def_readwrite("dV", &DataBlockHost::dV)
    .def_readwrite("A", &DataBlockHost::A)
    .def_readwrite("Vc", &DataBlockHost::Vc)
    .def_readwrite("dustVc", &DataBlockHost::dustVc)
    #if MHD == YES
      .def_readwrite("Vs", &DataBlockHost::Vs)
      .def_readwrite("Ve", &DataBlockHost::Ve)
      .def_readwrite("J", &DataBlockHost::J)
      .def_readwrite("Ex1", &DataBlockHost::Ex1)
      .def_readwrite("Ex2", &DataBlockHost::Ex2)
      .def_readwrite("Ex3", &DataBlockHost::Ex3)
    #endif
    .def_readwrite("InvDt", &DataBlockHost::InvDt)
    .def_readwrite("t",&DataBlockHost::t)
    .def_readwrite("dt",&DataBlockHost::dt);

    m.attr("RHO") = RHO;
    m.attr("VX1") = VX1;
    m.attr("VX2") = VX2;
    m.attr("VX3") = VX3;
    m.attr("PRS") = PRS;
    #if MHD == YES
      m.attr("BX1") = BX1;
      m.attr("BX2") = BX2;
      m.attr("BX3") = BX3;
      D_EXPAND(
        m.attr("BX1s") = BX1s; ,
        m.attr("BX2s") = BX2s; ,
        m.attr("BX3s") = BX3s; )
    #endif
    m.attr("IDIR") = IDIR;
    m.attr("JDIR") = JDIR;
    m.attr("KDIR") = KDIR;

    m.attr("prank") = idfx::prank;
    m.attr("psize") = idfx::psize;

    m.def("GatherIdefixArray",&GatherIdefixArray, "Gather arrays from domain decomposition");
}


template<typename... Ts>
void Pydefix::CallScript(std::string scriptName, std::string funcName, Ts... args) {
  idfx::pushRegion("Pydefix::CallScript");
  try {
    py::module_ script = py::module_::import(scriptName.c_str());
    py::object result = script.attr(funcName.c_str())(&args...);
  } catch(std::exception &e) {
    std::stringstream message;
    message << "An exception occured while calling the Python interpreter" << std::endl
                << "in file \"" << scriptName << ".py\" function \"" << funcName << "\":"
                << std::endl
                << e.what() << std::endl;
    IDEFIX_ERROR(message);
  }
  idfx::popRegion();
}


Pydefix::Pydefix(Input &input) {
  // Check that the input has a [python] block
  if(input.CheckBlock("Python")) {
    this->isActive = true;
    ninstance++;
    // Check whether we need to start an interpreter
    if(ninstance==1) {
      idfx::cout << "Pydefix: start Python interpreter." << std::endl;

      py::initialize_interpreter();
    }
    this->scriptFilename = input.Get<std::string>("Python","script",0);
    if(scriptFilename.substr(scriptFilename.length() - 3, 3).compare(".py")==0) {
      IDEFIX_ERROR("You should not include the python script .py extension in your input file");
    }
    if(input.CheckEntry("Python","output_function")>0) {
      this->haveOutput = true;
      this->outputFunctionName = input.Get<std::string>("Python","output_function",0);
    }
    if(input.CheckEntry("Python","initflow_function")>0) {
      this->haveInitflow = true;
      this->initflowFunctionName = input.Get<std::string>("Python","initflow_function",0);
    }
  }
}

void Pydefix::Output(DataBlock &data, int n) {
  idfx::pushRegion("Pydefix::Output");
  if(!this->isActive) {
    IDEFIX_ERROR("Python Outputs requires the [python] block to be defined in your input file.");
  }
  if(!this->haveOutput) {
    IDEFIX_ERROR("No python output function has been defined "
                  "in your input file [python]:output_function");
  }
  DataBlockHost dataHost(data);
  dataHost.SyncFromDevice();
  this->CallScript(this->scriptFilename,this->outputFunctionName,dataHost, n);
  idfx::popRegion();
}

void Pydefix::InitFlow(DataBlock &data) {
  idfx::pushRegion("Pydefix::InitFlow");
  if(!this->isActive) {
    IDEFIX_ERROR("Python Initflow requires the [python] block to be defined in your input file.");
  }
  if(!this->haveInitflow) {
    IDEFIX_ERROR("No python initflow function has been defined "
                  "in your input file [python]:initflow_function");
  }
  DataBlockHost dataHost(data);
  dataHost.SyncFromDevice();
  this->CallScript(this->scriptFilename,this->initflowFunctionName,dataHost);
  dataHost.SyncToDevice();
  idfx::popRegion();
}

void Pydefix::ShowConfig() {
  if(isActive == false) {
    idfx::cout << "Pydefix: DISABLED." << std::endl;
  } else {
    idfx::cout << "Pydefix: ENABLED." << std::endl;
    if(haveOutput) {
      idfx::cout << "Pydefix: output function ENABLED." << std::endl;
    } else {
      idfx::cout << "Pydefix: output function DISABLED." << std::endl;
    }
    if(haveInitflow) {
      idfx::cout << "Pydefix: initflow function ENABLED." << std::endl;
    } else {
      idfx::cout << "Pydefix: initflow function DISABLED." << std::endl;
    }
  }
}

Pydefix::~Pydefix() {
  if(isActive) {
    if(ninstance == 1) {
      py::finalize_interpreter();
      idfx::cout << "Pydefix: shutdown Python interpreter." << std::endl;
    }
    ninstance--;
    isActive = false;
  }
}


/*
py::array_t<real> Pydefix::toNumpyArray(const IdefixHostArray3D<real>& in) {
  py::array_t<real, py::array::c_style> array({in.extent(0),in.extent(1),in.extent(2)},in.data());
  return array;
}

py::array_t<real> Pydefix::toNumpyArray(const IdefixHostArray4D<real>& in) {
  py::array_t<real, py::array::c_style> array({in.extent(0),in.extent(1),in.extent(2),in.extent(3)},in.data());
  return array;
}
*/
