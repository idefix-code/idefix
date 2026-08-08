// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include <memory>
#include <string>

#include "selfGravity.hpp"
#include "selfGravityIterative.hpp"
#ifdef WITH_FFT
  #include "selfGravityFFT.hpp"
#endif

std::unique_ptr<SelfGravity> SelfGravity::Create(Input &input, DataBlock *data) {
  std::string strSolver = input.GetOrSet<std::string>("SelfGravity","solver",0,"BICGSTAB");

  std::unique_ptr<SelfGravity> ptr;
  if(strSolver == "FFT" || strSolver == "fft") {
    #ifdef WITH_FFT
      ptr = std::make_unique<SelfGravityFFT>();
    #else
      IDEFIX_ERROR("[SelfGravity]: FFT solver requested but Idefix was not compiled with FFT support.");
    #endif
  } else {
    ptr = std::make_unique<SelfGravityIterative>();
  }

  ptr->Init(input, data);
  return ptr;
}
