// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

// This file is originally from the ASTRA code
// https://github.com/glesur/astra

#ifndef UTILS_FFT_TRANSPOSE_HPP_
#define UTILS_FFT_TRANSPOSE_HPP_

#ifdef WITH_MPI
#include <mpi.h>
#endif
#include <iostream>
#include <Kokkos_Core.hpp>


template <typename T>
class Transpose {
 public:
  explicit Transpose(std::array<int,3> n) {
    // Allocate temporary arrays for domain-splited FFTs and transposes
    int64_t n1 = n[0]*n[1]; // Block size
    int64_t nk = n[2];
    this->tempB = Kokkos::View<T**, Kokkos::LayoutRight, Device>("FFT transpose tempB", n1,nk);
    this->tempC = Kokkos::View<T**, Kokkos::LayoutRight, Device>("FFT transpose tempC", n1,nk);
  }
  void Apply(const IdefixArray3D<T>& input, const IdefixArray3D<T>& output);
  void Test();

 private:
  Kokkos::View<T**, Kokkos::LayoutRight, Device> tempB, tempC;
};



/* MPI Transposition routines.
From complex to real, we first transform along x2
Then we transpose x2 and x1 to have x1 contiguous in memory
We then transform along x1 and finally x3 (this last one being a real fft). This gives a final

Transposition is done as follows:
Star with a matrix A_ijk
we require the first two indices to be dividible by n. so that (i,j) can be divided into n^2
blocks. Hence the indices reads

A_mi'pj'k
where i'=i'/n and j'=j/n while m = i%n and p = j%n
Initially the matrix is distributed accross MPI processes along m,

The first step is to transpose locally pj'=j and i'
B_mpj'i'k = A_mi'pj'k

Then we do an MPI_Alltoall to exchange the data between processors, this transposes m and p
C_pmj'i'k = B_mpj'i'k

Finally, we do a local block transpose between j' and m
D_pj'mi'k = C_pmj'i'k

The Matrix D is the final result with mi' contiguous in memory and pj' distributed accross MPI processes.

Example:
Assuming a full domain of 6x4 elements.

Initial layout contiguous along the second index (x2) and decomposed along the first (x1):
The dots represent the block division inside each MPI process whch is just represented to guide the eye
A11 A12 . A13 A14
A21 A22 . A23 A24
A31 A32 . A33 A34
------------- MPI_Proc division -------------
A41 A42 . A43 A44
A51 A52 . A53 A54
A61 A62 . A63 A64

1st step, local transposition (in each processor):
A11 A21 A31
A12 A22 A32
 .   .   .
A13 A23 A33
A14 A24 A34
------------- MPI_Proc division -------------
A41 A51 A61
A42 A52 A62
 .   .   .
A43 A53 A63
A44 A54 A64


Then we do a MPI_AllToall to exchange the data between processors.
A11 A21 A31
A12 A22 A32
 .   .   .
A41 A51 A61
A42 A52 A62
------------- MPI_Proc division -------------
A13 A23 A33
A14 A24 A34
 .   .   .
A43 A53 A63
A44 A54 A64

Finally the block transpose
A11 A21 A31 . A41 A51 A61
A12 A22 A32 . A42 A52 A62
-------------- MPI_Proc division -------------
A13 A23 A33 . A43 A53 A63
A14 A24 A34 . A44 A54 A64


After local transpose:


From real to complex, we first transform along x3 (real fft) then x1 (which is transposed, so this
is really the second coordinate).
We then transpose x1 and x2 to have x2 contiguous in memory and finally transform along x2.


Then along x3 (real fft)
we first transform along x3 and x2.
Then we need to transpose x1 and x2 to have x1 contiguous in memory.


*/

#include "idefix.hpp"

template <typename T>
void Transpose<T>::Apply(const IdefixArray3D<T>& in, const IdefixArray3D<T>& out) {
  idfx::pushRegion("Transpose::Apply");
  #ifdef WITH_MPI
  const int64_t n = idfx::psize; // number of MPI processes
  const int64_t ni = in.extent(0); // Size of local block
  const int64_t nj = in.extent(1)/n;
  const int64_t nk = in.extent(2);
  const int64_t n1 = ni*nj*n; // Full size of the first dimension


  //Array2D<complex> tempB("FFT transpose temp", n1,nk);
  auto tempB = this->tempB;
  if(tempB.extent(0) != n1 || tempB.extent(1) != nk) {
    idfx:: cout << "!! tempB size: " << tempB.extent(0) << "," << tempB.extent(1) << std::endl;
    idfx:: cout << "!! in size: " << in.extent(0) << "," << in.extent(1) << "," << in.extent(2) << std::endl;
    throw std::runtime_error("Transpose::Apply: temporary array size does not match input size");
  }

  idefix_for("FFT::TransposeLocal1", 0,n1, 0, nk,
    KOKKOS_LAMBDA(const int64_t ij, const int64_t k) {
    // for in array, ij = j'+nj*p +i'*nj*n= j+i'*nj*n (given that j=j'+nj*p)
    // For tempB, ij=i'+ni*(j'+nj*p)
     int64_t i = ij / (nj*n);
     int64_t j = ij- i*nj*n;
     int64_t ijprime = i + ni*j;
     tempB(ijprime, k) = in(i,j,k);
    });
    Kokkos::fence();

  auto tempC = this->tempC;
  //Kokkos::View<complex**, Kokkos::LayoutRight, Device> tempC("FFT transpose temp", n1,nk);
  int ret = MPI_Alltoall(tempB.data(), ni*nj*nk*sizeof(T),
                         MPI_BYTE,
                         tempC.data(), ni*nj*nk*sizeof(T),
                         MPI_BYTE,
                         MPI_COMM_WORLD);

  idefix_for("FFT::TransposeBlock", 0,n1, 0, nk,
    KOKKOS_LAMBDA(const int64_t ij, const int64_t k) {
      // For tempC, ij = i' + ni*(j'+m*nj)
      // For out, ij = i' + m*ni + j'*ni*n (given that i = i'+m*ni)

      const int64_t j = ij / (ni*n);
      const int64_t i = ij - j*ni*n;
      const int64_t m = i / ni;
      const int64_t iprime = i - m*ni;
      const int64_t ijprime = iprime +ni*(j + m*nj);

      out(j,i,k) = tempC(ijprime,k);
    });
    Kokkos::fence();

    #endif
  idfx::popRegion();
}

template <typename T>
void Transpose<T>::Test() {
  idfx::cout << "Testing FFT Transpose" << std::endl;
  const int n = idfx::psize; // number of MPI processes
  const int rank = idfx::prank; // rank of the current process
  const int ni = 3; // Size of local block
  const int nj = 4;
  const int nk = 9;


  Transpose<T> myTranspose({ni,nj*n,nk});
  IdefixArray3D<T> in("FFT transpose test input", ni,nj*n,nk);
  IdefixArray3D<T> out("FFT transpose test output", nj,ni*n,nk);

  // Fill the input array with some test values
  idefix_for("FFT::TestTransposeFill", 0,ni, 0,nj*n, 0,nk,
    KOKKOS_LAMBDA(const int64_t i, const int64_t j, const int64_t k) {
      in(i,j,k) = k+10*j+100*(i+rank*ni);
    });
  myTranspose.Apply(in, out);
  IdefixHostArray3D<T> outHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),out);

  // Check the output values
  bool error = false;
  for(int64_t j = 0 ; j < nj ; j++) {
    for(int64_t i = 0 ; i < ni*n ; i++) {
      for(int64_t k = 0 ; k < nk ; k++) {
        if(outHost(j,i,k) != k+10*(j+rank*nj)+100*(i)) {
          idfx::cout << "rank " << rank << " Transpose error at (" << i << "," << j << "," << k
                      << "): got " << outHost(j,i,k) << " instead of " << k+10*(j+rank*nj)+100*(i) << std::endl;
          error = true;
        }
      }
    }
  }

  if(!error) {
    idfx::cout << "Transpose test passed" << std::endl;
  } else {
    throw std::runtime_error("Transpose test failed @ rank " + std::to_string(rank));
  }
}


#endif // UTILS_FFT_TRANSPOSE_HPP_
