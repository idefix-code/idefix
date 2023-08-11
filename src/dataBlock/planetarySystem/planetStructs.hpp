// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef DATABLOCK_PLANETARYSYSTEM_PLANETSTRUCTS_HPP_
#define DATABLOCK_PLANETARYSYSTEM_PLANETSTRUCTS_HPP_

#include "idefix.hpp"
#include "input.hpp"


// Structs used in Planetary system class

struct Force {
    real f_inner[3];
    real f_ex_inner[3];
    real f_outer[3];
    real f_ex_outer[3];
    KOKKOS_FUNCTION void operator+=(Force const volatile& f) volatile {
        for (int i = 0; i < 3; ++i) {
            f_inner[i] += f.f_inner[i];
            f_ex_inner[i] += f.f_ex_inner[i];
            f_outer[i] += f.f_outer[i];
            f_ex_outer[i] += f.f_ex_outer[i];
        }
    }
};
struct Point {
    real x;
    real y;
    real z;
};

struct PointSpeed {
  real x;
  real y;
  real z;
  real vx;
  real vy;
  real vz;
};

// arithmetics of pointspeed
// Addition
static PointSpeed operator+(PointSpeed a, PointSpeed const &b ) {
    a.x += b.x;
    a.y += b.y;
    a.z += b.z;
    a.vx += b.vx;
    a.vy += b.vy;
    a.vz += b.vz;
    return(a);
  }

// Substraction
static PointSpeed operator-(PointSpeed a, PointSpeed const &b ) {
    a.x -= b.x;
    a.y -= b.y;
    a.z -= b.z;
    a.vx -= b.vx;
    a.vy -= b.vy;
    a.vz -= b.vz;
    return(a);
  }

// product
template<typename T>
static PointSpeed operator*(PointSpeed a, T const &b ) {
  a.x *= b;
  a.y *= b;
  a.z *= b;
  a.vx *= b;
  a.vy *= b;
  a.vz *= b;
  return(a);
}
template<typename T>
static PointSpeed operator*(T const &b, PointSpeed a) {
  a.x *= b;
  a.y *= b;
  a.z *= b;
  a.vx *= b;
  a.vy *= b;
  a.vz *= b;
  return(a);
}

// division
template<typename T>
static PointSpeed operator/(PointSpeed a, T const &b ) {
  a.x /= b;
  a.y /= b;
  a.z /= b;
  a.vx /= b;
  a.vy /= b;
  a.vz /= b;
  return(a);
}

#endif //DATABLOCK_PLANETARYSYSTEM_PLANETSTRUCTS_HPP_
