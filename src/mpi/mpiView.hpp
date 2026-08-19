// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef MPI_MPIVIEW_HPP_
#define MPI_MPIVIEW_HPP_

#ifdef WITH_MPI
  #include <mpi.h>
#endif //WITH_MPI

#include <cassert>
#include <vector>
#include <array>

#include "idefix.hpp"
#include "commArray.hpp"

namespace idefix {

#ifdef WITH_MPI

/**
 * Provide a wrapper of the standard MPI_Sendrecv by handling an idefix
 * communication array. It will automate the process of transfering the data
 * to the host if required (when WITH_MPI_GPU_DIRECT is disabled).
 */
template <class T>
int MPI_Sendrecv(IdefixCommArray<T> & sendbuf, int sendcount, MPI_Datatype sendtype,
  int dest, int sendtag, IdefixCommArray<T> &recvbuf, int recvcount,
  MPI_Datatype recvtype, int source, int recvtag, MPI_Comm comm, MPI_Status *status) {
  //check size
  assert(sendbuf.span() == sendcount);
  assert(sendbuf.span_is_contiguous());
  assert(recvbuf.span() == sendcount);
  assert(recvbuf.span_is_contiguous());

  //make a local allocation if needed
  sendbuf.syncCommData();

  //make the communication
  const int result = ::MPI_Sendrecv(sendbuf.commData(), sendcount, sendtype, dest, sendtag,
               recvbuf.commData(), recvcount, recvtype, source, recvtag,
               comm, status);

  //copy back if needed
  recvbuf.syncDeviceData();

  //ok
  return result;
}

/**
 * Overlad the MPI_Send_init() function to handle an idefix communication array. It will automate the
 * CPU / GPU transfers if the array in on GPU and WITH_MPI_GPU_DIRECT is disabled.
 */
template <class T>
int MPI_Send_init(IdefixCommArray<T> & buf, int count, MPI_Datatype datatype, int dest,
                  int tag, MPI_Comm comm, Idefix_MPI_Request<T> *request) {
    //check
    assert(request != nullptr);
    assert(buf.span() == count);
    assert(buf.span_is_contiguous());

    //store the values
    request->commArray = buf;

    //call MPI
    const int result = ::MPI_Send_init(buf.commData(), count, datatype,
      dest, tag, comm, &request->request);

    //ok
    return result;
}

/**
 * Overlad the MPI_Recv_init() function to handle an idefix communication array. It will automate the
 * CPU / GPU transfers if the array in on GPU and WITH_MPI_GPU_DIRECT is disabled.
 */
template <class T>
int MPI_Recv_init(IdefixCommArray<T> & buf, int count, MPI_Datatype datatype, int
    source, int tag, MPI_Comm comm, Idefix_MPI_Request<T>* request) {
  //check
  assert(request != nullptr);
  assert(buf.span() == count);
  assert(buf.span_is_contiguous());

  //store the values
  request->commArray = buf;

  //call MPI
  const int result = ::MPI_Recv_init(buf.commData(), count, datatype,
    source, tag, comm, &request->request);

  //ok
  return result;
}

/**
 * Overlad the MPI_Request_free() function to handle an idefix communication array. It will automate the
 * CPU / GPU transfers if the array in on GPU and WITH_MPI_GPU_DIRECT is disabled.
 */
template <class T>
int MPI_Request_free(Idefix_MPI_Request<T>* request) {
  //check
  assert(request != nullptr);

  //call MPI
  return ::MPI_Request_free(&request->request);
}

/**
 * Overlad the MPI_Startall() function to handle an idefix communication array. It will automate the
 * CPU / GPU transfers if the array in on GPU and WITH_MPI_GPU_DIRECT is disabled.
 */
template <class T>
int MPI_Startall(int count, Idefix_MPI_Request<T>* request) {
  //check
  assert(request != nullptr);
  assert(count >= 0);

  //sync the copies in async mode
  for (size_t i = 0 ; i < count ; i++) {
    request[i].commArray.syncCommDataAsync();
  }

  //wait transfers to be done
  Kokkos::fence();

  //launch the comms
  int final_result = 0;
  for (size_t i = 0 ; i < count ; i++) {
    const int result = ::MPI_Start(&request[i].request);
    if (result < 0)
      final_result = result;
  }

  //ok
  return final_result;
}

/**
 * Overlad the MPI_Waitall() function to handle an idefix communication array. It will automate the
 * CPU / GPU transfers if the array in on GPU and WITH_MPI_GPU_DIRECT is disabled.
 */
template <class T>
int MPI_Waitall(int count, Idefix_MPI_Request<T>* request,
  MPI_Status *array_of_statuses) {
  //check
  assert(request != nullptr);
  assert(count >= 0);

  //wait the comms
  int final_result = 0;
  for (size_t i = 0 ; i < count ; i++) {
    const int result = ::MPI_Wait(&request[i].request, &array_of_statuses[i]);
    if (result < 0)
      final_result = result;
  }

  //sync the copies in async mode
  for (size_t i = 0 ; i < count ; i++) {
    request[i].commArray.syncDeviceDataAsync();
  }

  //wait transfers to be done
  Kokkos::fence();

  //ok
  return final_result;
}

/**
 * Overlad the MPI_Start() function to handle an idefix communication array. It will automate the
 * CPU / GPU transfers if the array in on GPU and WITH_MPI_GPU_DIRECT is disabled.
 */
template <class T>
int MPI_Start(Idefix_MPI_Request<T>* request) {
  //check
  assert(request != nullptr);

  //sync the copies in async mode
  request->commArray.syncCommData();

  //call MPI
  return ::MPI_Start(&request->request);
}

/**
 * Overlad the MPI_Wait() function to handle an idefix communication array. It will automate the
 * CPU / GPU transfers if the array in on GPU and WITH_MPI_GPU_DIRECT is disabled.
 */
template <class T>
int MPI_Wait(Idefix_MPI_Request<T>* request, MPI_Status *array_of_statuses) {
  //check
  assert(request != nullptr);

  //wait the comms
  const int result = ::MPI_Wait(&request->request, array_of_statuses);

  //sync the copies in async mode
  request->commArray.syncDeviceData();

  //ok
  return result;
}

/**
 * Overlad the MPI_Send() function to handle an idefix communication array. It will automate the
 * CPU / GPU transfers if the array in on GPU and WITH_MPI_GPU_DIRECT is disabled.
 */
template <class T>
int MPI_Send(IdefixCommArray<T> & buf, int count, MPI_Datatype datatype, int
    dest, int tag, MPI_Comm comm) {
  //check
  assert(buf.span() == count);
  assert(buf.span_is_contiguous());

  //sync
  buf.syncCommData();

  //call MPI
  const int result = ::MPI_Send(buf.commData(), count, datatype,
    dest, tag, comm);

  //ok
  return result;
}

/**
 * Overlad the MPI_Send() function to handle an idefix communication array. It will automate the
 * CPU / GPU transfers if the array in on GPU and WITH_MPI_GPU_DIRECT is disabled.
 */
template <class T, size_t U>
int MPI_Send(const std::array<T, U> & buf, int count, MPI_Datatype datatype,
    int dest, int tag, MPI_Comm comm) {
  //call MPI
  return ::MPI_Send(buf.data(), count, datatype, dest, tag, comm);
}

/**
 * Overlad the MPI_Send() function to handle an idefix communication array. It will automate the
 * CPU / GPU transfers if the array in on GPU and WITH_MPI_GPU_DIRECT is disabled.
 */
template <class T>
int MPI_Send(const Kokkos::View<T, Kokkos::LayoutRight, Kokkos::HostSpace> & buf, int count,
    MPI_Datatype datatype, int dest, int tag, MPI_Comm comm) {
  //check
  assert(buf.span() == count);
  assert(buf.span_is_contiguous());

  //call MPI
  return ::MPI_Send(buf.data(), count, datatype, dest, tag, comm);
}

/**
 * Overlad the MPI_Recv() function to handle an idefix communication array. It will automate the
 * CPU / GPU transfers if the array in on GPU and WITH_MPI_GPU_DIRECT is disabled.
 */
template <class T>
int MPI_Recv(IdefixCommArray<T> & buf, int count, MPI_Datatype datatype, int
    source, int tag, MPI_Comm comm, MPI_Status * status) {
  //check
  assert(buf.span() == count);
  assert(buf.span_is_contiguous());

  //call MPI
  const int result = ::MPI_Recv(buf.commData(), count, datatype,
    source, tag, comm, status);

  //sync
  buf.syncDeviceData();

  //ok
  return result;
}

/**
 * Overlad the MPI_Recv() function to handle an idefix communication array. It will automate the
 * CPU / GPU transfers if the array in on GPU and WITH_MPI_GPU_DIRECT is disabled.
 */
template <class T, size_t U>
int MPI_Recv(std::array<T, U> & buf, int count, MPI_Datatype datatype, int
    source, int tag, MPI_Comm comm, MPI_Status * status) {
  //call MPI
  return ::MPI_Recv(buf.data(), count, datatype, source, tag, comm, status);
}

/**
 * Overlad the MPI_Recv() function to handle an idefix communication array. It will automate the
 * CPU / GPU transfers if the array in on GPU and WITH_MPI_GPU_DIRECT is disabled.
 */
template <class T>
int MPI_Recv(Kokkos::View<T, Kokkos::LayoutRight, Kokkos::HostSpace> & buf, int count,
    MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status * status) {
  //check
  assert(buf.span() == count);
  assert(buf.span_is_contiguous());

  //call MPI
  return ::MPI_Recv(buf.data(), count, datatype, source, tag, comm, status);
}

/**
 * Overlad the MPI_Isend() function to handle an idefix communication array. It will automate the
 * CPU / GPU transfers if the array in on GPU and WITH_MPI_GPU_DIRECT is disabled.
 */
template <class T>
int MPI_Isend(IdefixCommArray<T> & buf, int count, MPI_Datatype datatype, int
    source, int tag, MPI_Comm comm, Idefix_MPI_Request<T>* request) {
  //check
  assert(request != nullptr);
  assert(buf.span() == count);
  assert(buf.span_is_contiguous());

  //store the values
  request->commArray = buf;

  //sync
  buf.syncCommData();

  //call MPI
  const int result = ::MPI_Isend(buf.commData(), count, datatype,
    source, tag, comm, &request->request);

  //ok
  return result;
}

/**
 * Overlad the MPI_Irecv() function to handle an idefix communication array. It will automate the
 * CPU / GPU transfers if the array in on GPU and WITH_MPI_GPU_DIRECT is disabled.
 */
template <class T>
int MPI_Irecv(IdefixCommArray<T> & buf, int count, MPI_Datatype datatype, int
    source, int tag, MPI_Comm comm, Idefix_MPI_Request<T>* request) {
  //check
  assert(request != nullptr);
  assert(buf.span() == count);
  assert(buf.span_is_contiguous());

  //store the values
  request->commArray = buf;

  //call MPI
  const int result = ::MPI_Irecv(buf.commData(), count, datatype,
    source, tag, comm, &request->request);

  //ok
  return result;
}

/**
 * Overlad the MPI_Bcast() function to handle an idefix communication array. It will automate the
 * CPU / GPU transfers if the array in on GPU and WITH_MPI_GPU_DIRECT is disabled.
 */
template <class T>
int MPI_Bcast(IdefixCommArray<T> & buffer, int count, MPI_Datatype datatype, int root,
  MPI_Comm comm) {
  //check
  assert(buffer.span() == count);
  assert(buffer.span_is_contiguous());

  //get current rank
  int rank = -1;
  const int status = ::MPI_Comm_rank(comm, &rank);
  assert(status == MPI_SUCCESS);

  //sync
  if (rank == root)
    buffer.syncCommData();

  //call MPI
  const int result = ::MPI_Bcast(buffer.commData(), count, datatype, root, comm);

  //sync back
  if (rank != root)
    buffer.syncDeviceData();

  //ok
  return result;
}

/**
 * Overlad the MPI_Bcast() function to handle an idefix communication array. It will automate the
 * CPU / GPU transfers if the array in on GPU and WITH_MPI_GPU_DIRECT is disabled.
 */
template <class T>
int MPI_Bcast(IdefixHostArray1D<T> & buffer, int count, MPI_Datatype datatype, int root,
  MPI_Comm comm) {
  //check
  assert(buffer.span() == count);
  assert(buffer.span_is_contiguous());

  //call MPI
  const int result = ::MPI_Bcast(buffer.data(), count, datatype, root, comm);

  //ok
  return result;
}

/**
 * Overlad the MPI_Allreduce() function to handle an idefix communication array. It will automate the
 * CPU / GPU transfers if the array in on GPU and WITH_MPI_GPU_DIRECT is disabled.
 */
template <class T>
int MPI_Allreduce(IdefixCommArray<T> & sendbuf, IdefixCommArray<T> & recvbuf, int count,
    MPI_Datatype datatype, MPI_Op op, MPI_Comm comm) {
  //check
  assert(sendbuf.span() == count);
  assert(sendbuf.span_is_contiguous());
  assert(recvbuf.span() == count);
  assert(recvbuf.span_is_contiguous());

  sendbuf.syncCommData();

  //call MPI
  const int result = ::MPI_Allreduce(sendbuf.commData(), recvbuf.commData(), count, datatype,
    op, comm);

  //sync back
  recvbuf.syncDeviceData();

  //ok
  return result;
}

/**
 * Overlad the MPI_Allreduce() function to handle an idefix communication array. It will automate the
 * CPU / GPU transfers if the array in on GPU and WITH_MPI_GPU_DIRECT is disabled.
 */
template <class T>
int MPI_Allreduce(void * sendbuf, IdefixHostArray3D<T> & recvbuf, int count,
    MPI_Datatype datatype, MPI_Op op, MPI_Comm comm) {
  //check
  assert(recvbuf.span() == count);
  assert(recvbuf.span_is_contiguous());

  //call MPI
  const int result = ::MPI_Allreduce(sendbuf, recvbuf.data(), count, datatype, op, comm);

  //ok
  return result;
}

/**
 * Overlad the MPI_Allreduce() function to handle an idefix communication array. It will automate the
 * CPU / GPU transfers if the array in on GPU and WITH_MPI_GPU_DIRECT is disabled.
 */
template <class T>
int MPI_Allreduce(void * sendbuf, IdefixCommArray<T> & recvbuf, int count,
    MPI_Datatype datatype, MPI_Op op, MPI_Comm comm) {
  //check
  assert(recvbuf.span() == count);
  assert(recvbuf.span_is_contiguous());

  //call MPI
  const int result = ::MPI_Allreduce(sendbuf, recvbuf.commData(), count, datatype, op, comm);

  //sync back
  recvbuf.syncDeviceData();

  //ok
  return result;
}

/**
 * Overlad the MPI_File_write_all() function to handle an idefix communication array. It will automate the
 * CPU / GPU transfers if the array in on GPU and WITH_MPI_GPU_DIRECT is disabled.
 */
template <class T>
int MPI_File_write_all(MPI_File fh, const IdefixHostArray4D<T> & buf, int count,
    MPI_Datatype datatype, MPI_Status* status) {
  //check
  assert(buf.span() == count);
  assert(buf.span_is_contiguous());

  //call MPI
  const int result = ::MPI_File_write_all(fh, buf.data(), count, datatype, status);

  //ok
  return result;
}

/**
 * Overlad the MPI_Gather() function to handle an idefix communication array. It will automate the
 * CPU / GPU transfers if the array in on GPU and WITH_MPI_GPU_DIRECT is disabled.
 */
template <class T>
int MPI_Gather(const void* sendbuf, int sendcount, MPI_Datatype
    sendtype, std::vector<T> & recvbuf, int recvcount, MPI_Datatype recvtype,
    int root, MPI_Comm comm) {
  assert(recvcount <= recvbuf.size());
  return ::MPI_Gather(sendbuf, sendcount, sendtype, recvbuf.data(), recvcount,
    recvtype, root, comm);
}

#endif //WITH_MPI

} // namespace idefix

#endif // MPI_MPIVIEW_HPP_
