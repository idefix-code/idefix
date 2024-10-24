// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef UTILS_BIGENDIAN_HPP_
#define UTILS_BIGENDIAN_HPP_

#include <cstddef>
#include <type_traits>
#include "idefix.hpp"

/* ****************************************************************************/
/** Determines if the machine is little-endian.  If so,
  it will force the data to be big-endian.
@param in_number floating point number to be converted in big endian */
/* *************************************************************************** */

class BigEndian {
 public:
  BigEndian() {
    // Test endianness
    union {
        uint32_t i;
        char c[4];
    } bint = {0x01020304};

    if(bint.c[0] == 1)
      this->shouldSwapEndian = false;
    else
      this->shouldSwapEndian = true;
  }

  // Swap when needed
  template <class T>
  T operator() (T in_number) {
    static_assert(std::is_arithmetic_v<T> == true);
    T out_number;
    if (this->shouldSwapEndian) {
      constexpr int size = sizeof(T);
      union {
        T u;
        unsigned char byte[size];
      } in, out;
      in.u = in_number;
      for(int n = 0 ; n < size ; n++) {
        out.byte[size-n-1] = in.byte[n];
      }
      out_number = out.u;
    } else {
      out_number = in_number;
    }
    return(out_number);
  }

 private:
  // Endianness swaping flag
  bool shouldSwapEndian;
};

#endif // UTILS_BIGENDIAN_HPP_
