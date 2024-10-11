// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

// Most of these macros are identical to the one defined in the Pluto code (c) A. Mignone & col.

#ifndef MACROS_HPP_
#define MACROS_HPP_

#include <cstdio>

#if COMPONENTS == 1
  #define EXPAND(a,b,c) a
  #define SELECT(a,b,c) a
#endif

#if COMPONENTS == 2
  #define EXPAND(a,b,c) a b
  #define SELECT(a,b,c) b
#endif

#if COMPONENTS == 3
  #define EXPAND(a,b,c) a b c
  #define SELECT(a,b,c) c

#endif

#if DIMENSIONS == 1
  #define D_EXPAND(a,b,c)  a
  #define D_SELECT(a,b,c)  a
  #define IOFFSET 1
  #define JOFFSET 0
  #define KOFFSET 0
#endif

#if DIMENSIONS == 2
  #define D_EXPAND(a,b,c) a b
  #define D_SELECT(a,b,c) b
  #define IOFFSET 1
  #define JOFFSET 1
  #define KOFFSET 0
#endif

#if DIMENSIONS == 3
  #define D_EXPAND(a,b,c) a b c
  #define D_SELECT(a,b,c) c
  #define IOFFSET 1
  #define JOFFSET 1
  #define KOFFSET 1
#endif

// Spatial averages macros.
//    The following set of macros provide a compact way to perform multi-D
//    averages from cell centered values to interfaces.
//    For instance, \C AVERAGE_X(q,k,j,i) will simply take the
//    arithmetic average betwen q(i-1) and q(i) at the i-1/2 interface.
//    Likewise, AVERAGE_YZ(q,k,j,i) will produce an average at the
//    j-1/2 and k-1/2 edge.


#define AVERAGE_3D_X(q,k,j,i)   (0.5*(q(k,j,i) + q(k,j,i-1)))

#if DIMENSIONS == 1

  #define AVERAGE_3D_Y(q,k,j,i)    (q(0,0,i))
  #define AVERAGE_3D_Z(q,k,j,i)    (q(0,0,i))

  #define AVERAGE_3D_XY(q,k,j,i)   AVERAGE_3D_X(q,0,0,i)
  #define AVERAGE_3D_XZ(q,k,j,i)   AVERAGE_3D_X(q,0,0,i)
  #define AVERAGE_3D_YZ(q,k,j,i)   (q(0,0,i))

  #define AVERAGE_3D_XYZ(q,k,j,i)  0.5*(q(0,0,i) + q(0,0,i-1))

#elif DIMENSIONS == 2

  #define AVERAGE_3D_Y(q,k,j,i)   (0.5*(q(k,j,i) + q(k,j-1,i)))
  #define AVERAGE_3D_Z(q,k,j,i)    (q(0,j,i))

  #define AVERAGE_3D_XY(q,k,j,i)   ( 0.25*(  q(k,j,i)   + q(k,j,i-1)      \
                                        + q(k,j-1,i) + q(k,j-1,i-1)) )
  #define AVERAGE_3D_XZ(q,k,j,i)   (0.5*(q(0,j,i) + q(0,j,i-1)))
  #define AVERAGE_3D_YZ(q,k,j,i)   (0.5*(q(0,j,i) + q(0,j-1,i)))

  #define AVERAGE_3D_XYZ(q,k,j,i)  (0.25*(  q(0,j,i)   + q(0,j,i-1)       \
                                       + q(0,j-1,i) + q(0,j-1,i-1)))

#elif DIMENSIONS == 3

  #define AVERAGE_3D_Y(q,k,j,i)   (0.5*(q(k,j,i) + q(k,j-1,i)))
  #define AVERAGE_3D_Z(q,k,j,i)   (0.5*(q(k,j,i) + q(k-1,j,i)))

  #define AVERAGE_3D_XY(q,k,j,i)  (0.25*(  q(k,j,i)   + q(k,j,i-1)        \
                                       + q(k,j-1,i) + q(k,j-1,i-1)) )
  #define AVERAGE_3D_XZ(q,k,j,i)  (0.25*(  q(k,j,i)   + q(k,j,i-1)        \
                                       + q(k-1,j,i) + q(k-1,j,i-1)) )
  #define AVERAGE_3D_YZ(q,k,j,i)  (0.25*(  q(k,j,i)   + q(k,j-1,i)        \
                                       + q(k-1,j,i) + q(k-1,j-1,i)) )
  #define AVERAGE_3D_XYZ(q,k,j,i) (0.125*(  q(k,j,i)   + q(k,j,i-1)       \
                                       + q(k,j-1,i) + q(k,j-1,i-1)        \
                                       + q(k-1,j,i)   + q(k-1,j,i-1)      \
                                       + q(k-1,j-1,i) + q(k-1,j-1,i-1)))
#endif

// Same but for 4D views
#define AVERAGE_4D_X(q,n,k,j,i)   (0.5*(q(n,k,j,i) + q(n,k,j,i-1)))

#if DIMENSIONS == 1

  #define AVERAGE_4D_Y(q,n,k,j,i)    (q(n,0,0,i))
  #define AVERAGE_4D_Z(q,n,k,j,i)    (q(n,0,0,i))

  #define AVERAGE_4D_XY(q,n,k,j,i)   AVERAGE_4D_X(q,n,0,0,i)
  #define AVERAGE_4D_XZ(q,n,k,j,i)   AVERAGE_4D_X(q,n,0,0,i)
  #define AVERAGE_4D_YZ(q,n,k,j,i)   (q(n,0,0,i))

  #define AVERAGE_4D_XYZ(q,n,k,j,i)  0.5*(q(n,0,0,i) + q(n,0,0,i-1))

#elif DIMENSIONS == 2


  #define AVERAGE_4D_Y(q,n,k,j,i)   (0.5*(q(n,k,j,i) + q(n,k,j-1,i)))
  #define AVERAGE_4D_Z(q,n,k,j,i)    (q(n,0,j,i))

  #define AVERAGE_4D_XY(q,n,k,j,i)   ( 0.25*(  q(n,k,j,i)   + q(n,k,j,i-1)    \
                                        + q(n,k,j-1,i) + q(n,k,j-1,i-1)) )
  #define AVERAGE_4D_XZ(q,n,k,j,i)   (0.5*(q(n,0,j,i) + q(n,0,j,i-1)))
  #define AVERAGE_4D_YZ(q,n,k,j,i)   (0.5*(q(n,0,j,i) + q(n,0,j-1,i)))

  #define AVERAGE_4D_XYZ(q,n,k,j,i)  (0.25*(  q(n,0,j,i)   + q(n,0,j,i-1)    \
                                       + q(n,0,j-1,i) + q(n,0,j-1,i-1)))

#elif DIMENSIONS == 3

  #define AVERAGE_4D_Y(q,n,k,j,i)   (0.5*(q(n,k,j,i) + q(n,k,j-1,i)))
  #define AVERAGE_4D_Z(q,n,k,j,i)   (0.5*(q(n,k,j,i) + q(n,k-1,j,i)))

  #define AVERAGE_4D_XY(q,n,k,j,i)  (0.25*(  q(n,k,j,i)   + q(n,k,j,i-1)     \
                                       + q(n,k,j-1,i) + q(n,k,j-1,i-1)) )
  #define AVERAGE_4D_XZ(q,n,k,j,i)  (0.25*(  q(n,k,j,i)   + q(n,k,j,i-1)     \
                                       + q(n,k-1,j,i) + q(n,k-1,j,i-1)) )
  #define AVERAGE_4D_YZ(q,n,k,j,i)  (0.25*(  q(n,k,j,i)   + q(n,k,j-1,i)     \
                                       + q(n,k-1,j,i) + q(n,k-1,j-1,i)) )
  #define AVERAGE_4D_XYZ(q,n,k,j,i) (0.125*(  q(n,k,j,i)   + q(n,k,j,i-1)    \
                                       + q(n,k,j-1,i) + q(n,k,j-1,i-1)       \
                                       + q(n,k-1,j,i)   + q(n,k-1,j,i-1)     \
                                       + q(n,k-1,j-1,i) + q(n,k-1,j-1,i-1)))
#endif



#define MPI_SAFE_CALL(cmd) {                                          \
  int mpiErrNo = (cmd);                                               \
  if (MPI_SUCCESS != mpiErrNo) {                                      \
    char msg[MPI_MAX_ERROR_STRING];                                   \
    int len;                                                          \
    std::stringstream stream;                                         \
    MPI_Error_string(mpiErrNo, msg, &len);                            \
    stream << "MPI failed with error code :" << mpiErrNo              \
                        << " " << msg << std::endl;                   \
    IDEFIX_ERROR(stream);                                             \
  }                                                                   \
}


#ifdef RUNTIME_CHECKS
  #define RUNTIME_CHECK_HOST(condition, ...) { \
    if(!(condition)) { \
      char msg[255]; \
      snprintf(msg, sizeof(msg), __VA_ARGS__); \
      IDEFIX_ERROR(msg); \
    } \
  }

  #if defined(KOKKOS_ENABLE_HIP) || defined(KOKKOS_ENABLE_CUDA)
    // string formatting functions can't be accessed from GPU kernels,
    // so we fallback to a simple error message
    #define RUNTIME_CHECK_KERNEL(condition, fallback_msg, ...) { \
      if(!(condition)) Kokkos::abort(fallback_msg); \
    }
  #else
    #define RUNTIME_CHECK_KERNEL(condition, fallback_msg, ...) { \
      if(!(condition)) { \
        char msg[255]; \
        snprintf(msg, sizeof(msg), __VA_ARGS__); \
        Kokkos::abort(msg); \
      } \
    }
  #endif // if defined(KOKKOS_ENABLE_HIP) || defined(KOKKOS_ENABLE_CUDA)

#else
  #define RUNTIME_CHECK_HOST(condition, msg, ...) {}
  #define RUNTIME_CHECK_KERNEL(condition, fallback_msg, ...) {}
#endif // ifdef RUNTIME_CHECKS


#endif // MACROS_HPP_
