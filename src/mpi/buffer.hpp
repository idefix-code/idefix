// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef MPI_BUFFER_HPP_
#define MPI_BUFFER_HPP_

#include "idefix.hpp"
#include "arrays.hpp"

using BoundingBox = std::array<std::array<int,2>,3>;


class Buffer {
 public:
  Buffer() = default;
  explicit Buffer(size_t size): pointer{0}, array{IdefixArray1D<real>("BufferArray",size)} {};

  // Compute the size of a bounding box
  static size_t ComputeBoxSize(BoundingBox box) {
    const int ni = box[IDIR][1]-box[IDIR][0];
    const int ninj = (box[JDIR][1]-box[JDIR][0])*ni;
    const int ninjnk = (box[KDIR][1]-box[KDIR][0])*ninj;
    return(ninjnk);
  }

  void* data() {
    return(array.data());
  }

  int Size() {
    return(array.size());
  }

  void ResetPointer() {
    this->pointer = 0;
  }

  void Pack(IdefixArray3D<real>& in,  BoundingBox box) {
    const int ni = box[IDIR][1]-box[IDIR][0];
    const int ninj = (box[JDIR][1]-box[JDIR][0])*ni;
    const int ninjnk = (box[KDIR][1]-box[KDIR][0])*ninj;
    const int ibeg = box[IDIR][0];
    const int jbeg = box[JDIR][0];
    const int kbeg = box[KDIR][0];
    const int iend = box[IDIR][1];
    const int jend = box[JDIR][1];
    const int kend = box[KDIR][1];
    const int offset = this->pointer;

    auto arr = this->array;
    idefix_for("LoadBuffer3D",kbeg,kend,jbeg,jend,ibeg,iend,
      KOKKOS_LAMBDA (int k, int j, int i) {
      arr(i-ibeg + (j-jbeg)*ni + (k-kbeg)*ninj + offset ) = in(k,j,i);
    });

    // Update pointer
    this->pointer += ninjnk;
  }

  void Pack(IdefixArray4D<real>& in,
       const int var,
       BoundingBox box) {
    const int ni = box[IDIR][1]-box[IDIR][0];
    const int ninj = (box[JDIR][1]-box[JDIR][0])*ni;
    const int ninjnk = (box[KDIR][1]-box[KDIR][0])*ninj;
    const int ibeg = box[IDIR][0];
    const int jbeg = box[JDIR][0];
    const int kbeg = box[KDIR][0];
    const int iend = box[IDIR][1];
    const int jend = box[JDIR][1];
    const int kend = box[KDIR][1];
    const int offset = this->pointer;

    auto arr = this->array;
    idefix_for("LoadBuffer4D_var",kbeg,kend,jbeg,jend,ibeg,iend,
      KOKKOS_LAMBDA (int k, int j, int i) {
      arr(i-ibeg + (j-jbeg)*ni + (k-kbeg)*ninj + offset ) = in(var, k,j,i);
    });

    // Update pointer
    this->pointer += ninjnk;
  }

  void Pack(IdefixArray4D<real>& in,
       IdefixArray1D<int>& map,
        BoundingBox box) {
    const int ni = box[IDIR][1]-box[IDIR][0];
    const int ninj = (box[JDIR][1]-box[JDIR][0])*ni;
    const int ninjnk = (box[KDIR][1]-box[KDIR][0])*ninj;
    const int ibeg = box[IDIR][0];
    const int jbeg = box[JDIR][0];
    const int kbeg = box[KDIR][0];
    const int iend = box[IDIR][1];
    const int jend = box[JDIR][1];
    const int kend = box[KDIR][1];
    const int offset = this->pointer;
    auto arr = this->array;

    idefix_for("LoadBuffer4D_map",0,map.size(),
                             kbeg,kend,
                             jbeg,jend,
                             ibeg,iend,
      KOKKOS_LAMBDA (int n, int k, int j, int i) {
      arr(i-ibeg + (j-jbeg)*ni + (k-kbeg)*ninj + n*ninjnk + offset ) = in(map(n), k,j,i);
    });

    // Update pointer
    this->pointer += ninjnk*map.size();
  }

  void Unpack(IdefixArray3D<real>& out,
       BoundingBox box) {
    const int ni = box[IDIR][1]-box[IDIR][0];
    const int ninj = (box[JDIR][1]-box[JDIR][0])*ni;
    const int ninjnk = (box[KDIR][1]-box[KDIR][0])*ninj;
    const int ibeg = box[IDIR][0];
    const int jbeg = box[JDIR][0];
    const int kbeg = box[KDIR][0];
    const int iend = box[IDIR][1];
    const int jend = box[JDIR][1];
    const int kend = box[KDIR][1];
    const int offset = this->pointer;
    auto arr = this->array;

    idefix_for("UnLoadBuffer3D",kbeg,kend,jbeg,jend,ibeg,iend,
      KOKKOS_LAMBDA (int k, int j, int i) {
        out(k,j,i) = arr(i-ibeg + (j-jbeg)*ni + (k-kbeg)*ninj + offset );
    });

    // Update pointer
    this->pointer += ninjnk;
  }

  void Unpack(IdefixArray4D<real>& out,
       const int var,
        BoundingBox box) {
      const int ni = box[IDIR][1]-box[IDIR][0];
      const int ninj = (box[JDIR][1]-box[JDIR][0])*ni;
      const int ninjnk = (box[KDIR][1]-box[KDIR][0])*ninj;
      const int ibeg = box[IDIR][0];
      const int jbeg = box[JDIR][0];
      const int kbeg = box[KDIR][0];
      const int iend = box[IDIR][1];
      const int jend = box[JDIR][1];
      const int kend = box[KDIR][1];
    const int offset = this->pointer;

    auto arr = this->array;
    idefix_for("UnLoadBuffer4D_var",kbeg,kend,jbeg,jend,ibeg,iend,
      KOKKOS_LAMBDA (int k, int j, int i) {
        out(var,k,j,i) = arr(i-ibeg + (j-jbeg)*ni + (k-kbeg)*ninj + offset );
    });

    // Update pointer
    this->pointer += ninjnk;
  }

  void Unpack(IdefixArray4D<real>& out,
       IdefixArray1D<int>& map,
       BoundingBox box) {
    const int ni = box[IDIR][1]-box[IDIR][0];
    const int ninj = (box[JDIR][1]-box[JDIR][0])*ni;
    const int ninjnk = (box[KDIR][1]-box[KDIR][0])*ninj;
    const int ibeg = box[IDIR][0];
    const int jbeg = box[JDIR][0];
    const int kbeg = box[KDIR][0];
    const int iend = box[IDIR][1];
    const int jend = box[JDIR][1];
    const int kend = box[KDIR][1];
    const int offset = this->pointer;

    auto arr = this->array;
    idefix_for("UnLoadBuffer4D_map",0,map.size(),
                              kbeg,kend,
                              jbeg,jend,
                              ibeg,iend,
      KOKKOS_LAMBDA (int n, int k, int j, int i) {
        out(map(n),k,j,i) = arr(i-ibeg + (j-jbeg)*ni + (k-kbeg)*ninj + n*ninjnk + offset );
    });

    // Update pointer
    this->pointer += ninjnk*map.size();
  }


 private:
  size_t pointer;
  IdefixArray1D<real> array;
};

#endif // MPI_BUFFER_HPP_
