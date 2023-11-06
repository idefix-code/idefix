// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

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
    .def_readwrite("xr", &DataBlockHost::x)
    .def_readwrite("xl", &DataBlockHost::x)
    .def_readwrite("dx", &DataBlockHost::x)
    .def_readwrite("dV", &DataBlockHost::x)
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
    .def_readwrite("InvDt", &DataBlockHost::InvDt);
}


Pydefix::Pydefix() {
  ninstance++;
  if(ninstance==1) {
    idfx::cout << "Pydefix: start Python interpreter." << std::endl;

    py::initialize_interpreter();
  }
}

Pydefix::~Pydefix() {
  if(ninstance == 1) {
    py::finalize_interpreter();
    idfx::cout << "Pydefix: shutdown Python interpreter." << std::endl;
  }
  ninstance--;
}

void Pydefix::CallScript(DataBlock *data, std::string scriptName, std::string funcName) {
  idfx::pushRegion("Pydefix::CallScript");
  DataBlockHost d(*data);
  d.SyncFromDevice();
  try {
    //auto Vc = pydefix.toNumpyArray(d.Vc);
    py::module_ script = py::module_::import(scriptName.c_str());

    //py::module_ embeded = py::module_::import("embeded");
    //py::object myV = py::cast(Vc);
    py::object result = script.attr(funcName.c_str())(nCalls, d);
  } catch(std::exception &e) {
    std::stringstream message;
    message << "An exception occured while calling the Python interpreter" << std::endl
                << "in file \"" << scriptName << ".py\" function \"" << funcName << "\":"
                << std::endl
                << e.what() << std::endl;
    IDEFIX_ERROR(message);
  }
  nCalls++;
  idfx::popRegion();
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
