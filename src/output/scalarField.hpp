// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef OUTPUT_SCALARFIELD_HPP_
#define OUTPUT_SCALARFIELD_HPP_
#include "idefix.hpp"


// Forward class declaration

class ScalarField {
 public:
  enum Type {Device3D, Device4D, Host3D, Host4D};

  explicit ScalarField(IdefixArray4D<real>& in, const int varnum):
    d4Darray{in}, var{varnum}, type{Device4D} {};
  explicit ScalarField(IdefixHostArray4D<real>& in, const int varnum):
    h4Darray{in}, var{varnum}, type{Host4D} {};
  explicit ScalarField(IdefixArray3D<real>& in):
    d3Darray{in}, type{Device3D} {};
  explicit ScalarField(IdefixHostArray3D<real>& in):
    h3Darray{in}, type{Host3D} {};

  IdefixHostArray3D<real> GetHostField() const {
    if(type==Host3D) {
      return(h3Darray);
    } else if(type==Host4D) {
      IdefixHostArray3D<real> arr3D = Kokkos::subview(
                                      h4Darray, var, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
      return(arr3D);
    } else if(type==Device3D) {
      IdefixHostArray3D<real> arr3D = Kokkos::create_mirror(d3Darray);
      Kokkos::deep_copy(arr3D,d3Darray);
      return(arr3D);
    } else if(type==Device4D) {
      IdefixArray3D<real> arrDev3D = Kokkos::subview(
                                      d4Darray, var, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
      IdefixHostArray3D<real> arr3D = Kokkos::create_mirror(arrDev3D);
      Kokkos::deep_copy(arr3D,arrDev3D);
      return(arr3D);
    } else {
      IDEFIX_ERROR("unknown field");
      return(h3Darray);
    }
  }

 private:
  IdefixArray4D<real> d4Darray;
  IdefixArray3D<real> d3Darray;
  IdefixHostArray4D<real> h4Darray;
  IdefixHostArray3D<real> h3Darray;
  int var;
  Type type;
};

#endif // OUTPUT_SCALARFIELD_HPP_
