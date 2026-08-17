// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef MPI_COMMARRAY_HPP_
#define MPI_COMMARRAY_HPP_

//To enable to force transfers and validate that it is working in any conditions.
//#define WITH_GPU_FORCE_COPY

#ifdef WITH_MPI
  #include <mpi.h>
#endif //WITH_MPI

#include <cassert>
#include <vector>
#include <array>

#include "idefix.hpp"

namespace idefix {

/**
 * Provide a specific array used for MPI communications. It embeds the standard
 * array as usual, and also optionally a communication copy on the host, depending
 * if we have access to WITH_MPI_GPU_DIRECT or not.
 */
template <class T>
class IdefixCommArrayGpuDirect : public T {
 public:
  //inherit constructors
  using T::T;

  /**
   * Accès it as a simple device view for usiage in kokkos kernels
   * (because need to not carry the host part).
   */
  T & deviceView(void) {
    return *this;
  }

  /**
   * If needed, transfers the data from the device to the host to be ready to make
   * a communication.
   */
  void syncCommData(void) {}

  /**
   * If needed, transfers the data from the device to the host to be ready to make
   * a communication.
   */
  void syncCommDataAsync(void) {}

  /**
   * If needed, transerts the data to the device after a communication to be ready
   * to use it.
   */
  void syncDeviceData(void) {}

  /**
   * If needed, transerts the data to the device after a communication to be ready
   * to use it.
   */
  void syncDeviceDataAsync(void) {}

  /**
   * Get access to the communication array. It will be either on host or device
   * depending on the WITH_MPI_GPU_DIRECT configuration.
   */
  void * commData(void) {
      return this->data();
  }
};

/**
 * Provide a specific array used for MPI communications. It embeds the standard
 * array as usual, and also optionally a communication copy on the host, depending
 * if we have access to WITH_MPI_GPU_DIRECT or not.
 */
template <class T>
class IdefixCommArrayNoGpuDirect : public T {
 public:
  //inherit constructors
  using T::T;

  /**
   * Accès it as a simple device view for usiage in kokkos kernels
   * (because need to not carry the host part).
   */
  T & deviceView(void) {
    return *this;
  }

  /**
   * If needed, transfers the data from the device to the host to be ready to make
   * a communication.
   */
  void syncCommData(void) {
    Kokkos::deep_copy(this->commArray, *this);
  }

  /**
   * If needed, transfers the data from the device to the host to be ready to make
   * a communication.
   */
  void syncCommDataAsync(void) {
    Kokkos::deep_copy(Kokkos::DefaultExecutionSpace(), this->commArray, *this);
  }

  /**
   * If needed, transerts the data to the device after a communication to be ready
   * to use it.
   */
  void syncDeviceData(void) {
    Kokkos::deep_copy(*this, this->commArray);
  }

  /**
   * If needed, transerts the data to the device after a communication to be ready
   * to use it.
   */
  void syncDeviceDataAsync(void) {
    Kokkos::deep_copy(Kokkos::DefaultExecutionSpace(), *this, this->commArray);
  }

  /**
   * Get access to the communication array. It will be either on host or device
   * depending on the WITH_MPI_GPU_DIRECT configuration.
   */
  void * commData(void) {
    return this->commArray.data();
  }

 private:
  static typename T::host_mirror_type initCommArray(T & deviceArray) {
    #ifdef WITH_GPU_FORCE_COPY
      //for validation purpose, force the transfers by using a copy anytime.
      return Kokkos::create_mirror(deviceArray);
    #else
      //if no GPU DIRECT, it creates a copy only if deviceArray is on GPU,
      //if on host, it just reference it
      return Kokkos::create_mirror_view(deviceArray);
    #endif
  }

 private:
  /**
   * Buffer to contain a copy of the data on host if required. If WITH_MPI_GPU_DIRECT
   * is enabled it directly points the device buffer as no copy is required.
   */
  typename T::host_mirror_type commArray{initCommArray(*this)};
};

#ifdef WITH_MPI_GPU_DIRECT
  template <class T> using IdefixCommArray = IdefixCommArrayGpuDirect<T>;
#else
  template <class T> using IdefixCommArray = IdefixCommArrayNoGpuDirect<T>;
#endif

/** Wrap the idefix 1D array to use it for communication purpose. */
template <class T> using IdefixCommArray1D = IdefixCommArray< IdefixArray1D<T> >;
/** Wrap the idefix 2D array to use it for communication purpose. */
template <class T> using IdefixCommArray2D = IdefixCommArray< IdefixArray2D<T> >;
/** Wrap the idefix 3D array to use it for communication purpose. */
template <class T> using IdefixCommArray3D = IdefixCommArray< IdefixArray3D<T> >;
/** Wrap the idefix 4D array to use it for communication purpose. */
template <class T> using IdefixCommArray4D = IdefixCommArray< IdefixArray4D<T> >;

template <class T>
struct Idefix_MPI_Request {
  IdefixCommArray<T> commArray;
  #ifdef WITH_MPI
    MPI_Request request;
  #endif //WITH_MPI
};

/** Wrap the idefix 1D array to use it for communication purpose. */
template <class T> using MPI_Request_1D = Idefix_MPI_Request< IdefixArray1D<T> >;
/** Wrap the idefix 2D array to use it for communication purpose. */
template <class T> using MPI_Request_2D = Idefix_MPI_Request< IdefixArray2D<T> >;
/** Wrap the idefix 3D array to use it for communication purpose. */
template <class T> using MPI_Request_3D = Idefix_MPI_Request< IdefixArray3D<T> >;
/** Wrap the idefix 4D array to use it for communication purpose. */
template <class T> using MPI_Request_4D = Idefix_MPI_Request< IdefixArray4D<T> >;

} // namespace idefix

#endif // MPI_COMMARRAY_HPP_
