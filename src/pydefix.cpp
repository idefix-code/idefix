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
  if(!this->haveOutput) {
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
