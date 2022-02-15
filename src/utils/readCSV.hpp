// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2022 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef UTILS_READCSV_HPP_
#define UTILS_READCSV_HPP_

#include <string>
#include "idefix.hpp"
#include "readCSV.hpp"

class ReadCSV {
 public:
  ReadCSV(std::string filename, char delimiter);
  IdefixArray1D<real> x;
  IdefixArray1D<real> y;
  IdefixArray2D<real> data;
};


KOKKOS_FUNCTION real K_ReadCSVFetch(const real x, const real y,
                                    const IdefixArray1D<real> &xin,
                                    const IdefixArray1D<real> &yin,
                                    const IdefixArray2D<real> &datain);
#endif //UTILS_READCSV_HPP_
