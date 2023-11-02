// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include <pybind11/embed.h> // everything needed for embedding
#include <pybind11/numpy.h> // for numpy arrays
#include "idefix.hpp"
#include "pydefix.hpp"

namespace py = pybind11;

int Pydefix::ninstance = 0;

Pydefix::Pydefix() {
  ninstance++;
  if(ninstance==1) {
    idfx::cout << "Pydefix: start interpreter." << std::endl;

    py::initialize_interpreter();
  }
  py::print("Hello,from python"); // use the Python API
}

Pydefix::~Pydefix() {
  if(ninstance == 1) {
    py::finalize_interpreter();
    idfx::cout << "Pydefix: shutdown interpreter." << std::endl;
  }
  ninstance--;
  idfx::cout << "Pydefix: destroyed." << std::endl;
}

py::array_t<real> Pydefix::toNumpyArray(const IdefixHostArray3D<real>& in) {
  py::array_t<real, py::array::c_style> array({in.extent(0),in.extent(1),in.extent(2)},in.data());
  return array;
}

py::array_t<real> Pydefix::toNumpyArray(const IdefixHostArray4D<real>& in) {
  py::array_t<real, py::array::c_style> array({in.extent(0),in.extent(1),in.extent(2),in.extent(3)},in.data());
  return array;
}
