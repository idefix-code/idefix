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

class ReadCSVDeviceContext {
 public:
  IdefixArray1D<real> xin;
  IdefixArray1D<real> yin;
  IdefixArray2D<real> data;

  // Fetch function that should be called inside idefix_loop
  KOKKOS_INLINE_FUNCTION
  real Fetch(const real x, const real y) const {
    const int nx = xin.extent(0);
    const int ny = yin.extent(0);

    int endx = nx-1;
    int endy = ny-1;
    if(x < xin(0) || x >xin(endx)) {
      return(NAN);
    }
    if(y < yin(0) || y >yin(endy)) {
      return(NAN);
    }

    int i = static_cast<int> ( (x - xin(0)) / (xin(endx) - xin(0)) * nx );
    int j = static_cast<int> ( (y - yin(0)) / (yin(endx) - yin(0)) * ny );

    if(xin(i) > x || xin(i+1) < x) { // x-points are not evenly distributed
      i = 0;
      while(xin(i) > x) {
        i++;
      }
    }
    if(yin(j) > y || yin(j+1) < y) { // y-points are not evenly distributed
      j = 0;
      while(yin(j) > y) {
        j++;
      }
    }

    real delta_x = (x - xin(i)) / (xin(i+1) - xin(i));
    real delta_y = (y - yin(j)) / (yin(j+1) - yin(j));

    // Rectilinear interpolation, in principle only valid for regularly-spaced datapoints.
    real value =   (1.0-delta_x) * (1.0-delta_y) * data(i  ,j  )
                        + (1.0-delta_x) * (delta_y)     * data(i  ,j+1)
                        + (delta_x)     * (1.0-delta_y) * data(i+1,j  )
                        + (delta_x)     * (delta_y)     * data(i+1,j+1);
    return(value);
  }
};

class ReadCSV {
 public:
  ReadCSV(std::string filename, char delimiter);
  ReadCSVDeviceContext GetDeviceContext() const;
  IdefixArray1D<real> x;
  IdefixArray1D<real> y;
  IdefixArray2D<real> data;
};


#endif //UTILS_READCSV_HPP_
