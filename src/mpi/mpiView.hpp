// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef MPI_MPIVIEW_HPP_
#define MPI_MPIVIEW_HPP_

#include <mpi.h>

#include <cassert>

#include "idefix.hpp"
#include "buffer.hpp"

/**
 * Provide a specific array used for MPI communications. It embeds the standard
 * array as usual, and also carry optionally a communication copy on host depending
 * if we have access to WITH_MPI_GPU_DIRECT or not.
 */
template <class T>
class IdefixCommArray : public T {
 public:
  //inherit constructors
  using T::T;

  /**
   * If needed, transfers the data from the device to the host to be ready to make
   * a communication.
   */
  void syncCommData(void) {
    Kokkos::deep_copy(this->commArray, *this);
  }

  /**
   * If needed, transerts the data to the device after a communication to be ready
   * to use it.
   */
  void syncDeviceData(void) {
    Kokkos::deep_copy(*this, this->commArray);
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
    #ifdef WITH_MPI_GPU_DIRECT
      //in gpu direct, we use directly the device buffer
      //if running on host, it will also just use directly the host storage.
      return deviceArray;
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

/** Wrap the idefix 1D array to use it for communication purpose. */
template <class T> using IdefixCommArray1D = IdefixCommArray< IdefixArray1D<T> >;
/** Wrap the idefix 2D array to use it for communication purpose. */
template <class T> using IdefixCommArray2D = IdefixCommArray< IdefixArray2D<T> >;
/** Wrap the idefix 3D array to use it for communication purpose. */
template <class T> using IdefixCommArray3D = IdefixCommArray< IdefixArray3D<T> >;
/** Wrap the idefix 4D array to use it for communication purpose. */
template <class T> using IdefixCommArray4D = IdefixCommArray< IdefixArray4D<T> >;

/*template <class T>
struct IdefixMpiViewRequest
{
  Mpi
  MPI_Request request;
};*/

/**
 * Provide a wrapper of the standard MPI_Sendrecv by handling a kokkkos
 * view. Use it with caution as it will allocate/deallocate a temporary
 * host buffer if required (when WITH_MPI_GPU_DIRECT is disabled.).
 */
template <class T, class U, class V>
int idefix_MPI_View_Sendrecv(Kokkos::View<T, U, V> sendbuf, int sendcount, MPI_Datatype sendtype,
  int dest, int sendtag, Kokkos::View<T, U, V> recvbuf, int recvcount,
  MPI_Datatype recvtype, int source, int recvtag, MPI_Comm comm, MPI_Status *status) {
  //check size
  assert(sendbuf.span() == sendcount);
  assert(sendbuf.span_is_contiguous());
  assert(recvbuf.span() == sendcount);
  assert(recvbuf.span_is_contiguous());

  //make a local allocation if needed
  #ifdef WITH_MPI_GPU_DIRECT
    auto sendHostArray = sendbuf;
    auto recvHostArray = recvbuf;
  #else
    auto sendHostArray = Kokkos::create_mirror_view(sendbuf);
    auto recvHostArray = Kokkos::create_mirror_view(recvbuf);
    Kokkos::deep_copy(sendHostArray, sendbuf);
  #endif

  //make the communication
  const int result = MPI_Sendrecv(sendHostArray.data(), sendcount, sendtype, dest, sendtag,
               recvHostArray.data(), recvcount, recvtype, source, recvtag,
               comm, status);

  //copy back if needed
  #ifndef WITH_MPI_GPU_DIRECT
    Kokkos::deep_copy(recvbuf, recvHostArray);
  #endif

  //ok
  return result;
}

/**
 * Provide a wrapper of the standard MPI_Sendrecv by handling an idefix
 * communication array. It will automate the process of transfering the data
 * to the host if required (when WITH_MPI_GPU_DIRECT is disabled).
 */
template <class T>
int idefix_MPI_View_Sendrecv(IdefixCommArray<T> sendbuf, int sendcount, MPI_Datatype sendtype,
  int dest, int sendtag, IdefixCommArray<T> recvbuf, int recvcount,
  MPI_Datatype recvtype, int source, int recvtag, MPI_Comm comm, MPI_Status *status) {
  //check size
  assert(sendbuf.span() == sendcount);
  assert(sendbuf.span_is_contiguous());
  assert(recvbuf.span() == sendcount);
  assert(recvbuf.span_is_contiguous());

  //make a local allocation if needed
  #ifndef WITH_MPI_GPU_DIRECT
    sendbuf.syncCommData();
  #endif

  //make the communication
  const int result = MPI_Sendrecv(sendbuf.commData(), sendcount, sendtype, dest, sendtag,
               recvbuf.commData(), recvcount, recvtype, source, recvtag,
               comm, status);

  //copy back if needed
  #ifndef WITH_MPI_GPU_DIRECT
    recvbuf.syncDeviceData();
  #endif

  //ok
  return result;
}

#endif // MPI_MPIVIEW_HPP_
