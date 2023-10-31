// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include <pybind11/embed.h> // everything needed for embedding
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
