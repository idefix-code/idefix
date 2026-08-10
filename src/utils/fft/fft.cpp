// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

// This file is originally from the ASTRA code
// https://github.com/glesur/astra

#include <memory>
#include <KokkosFFT.hpp>
#include <Kokkos_Random.hpp>

#include "idefix.hpp"
#include "fft.hpp"
#include "global.hpp"
#include "loop.hpp"
#include "transpose.hpp"

template <typename T> void ShowExtent(T array) {
  for(int i = 0 ; i < 3 ; i++) {
    idfx::cout << array.extent(i) << "  ";
  }
  idfx::cout << std::endl;
}

  // Empty constructor
FFT::FFT() {};

FFT::FFT(std::array<int,3> npr_glob, std::array<int,3> npf_glob) {
    this->npf_glob = npf_glob;
    this->npr_glob = npr_glob;
    // Local dimensions
    // We assume a decomposition along the first dimension only for now
    this->npr = npr_glob;
    this->npf = npf_glob;
    this->npr_t = npr_glob;

    #ifdef WITH_MPI
      this->npr[0] = npr_glob[0]/idfx::psize;
      this->npf[0] = npf_glob[0]/idfx::psize;
      this->npr_t[0] = npr_glob[1]/idfx::psize;
      this->npr_t[1] = npr_glob[0];
    #endif

    tempReal = IdefixArray3D<real>("FFT temp real", npr[0],npr[1],npr[2]);
    tempComplex = IdefixArray3D<complex>("FFT temp complex", npf[0],npf[1],npf[2]);

    // Create the FFT plans
    this->r2cPlan = std::make_unique<PlanR2CType>(Kokkos::DefaultExecutionSpace(),
      tempReal, tempComplex, KokkosFFT::Direction::forward, std::array<int,3>{-3,-2,-1});
    this->c2rPlan = std::make_unique<PlanC2RType>(Kokkos::DefaultExecutionSpace(),
      tempComplex, tempReal, KokkosFFT::Direction::backward, std::array<int,3>{-3,-2,-1});

    #ifdef WITH_MPI
      // Allocate temporary arrays for domain-splited FFTs and transposes
      this->tempTransposedComplex = IdefixArray3D<complex>("FFT transpose temp", npf[1]/idfx::psize, npf[0]*idfx::psize, npf[2]);
      this->tempTransposedComplex2 = IdefixArray3D<complex>("FFT transpose temp2", npf[1]/idfx::psize, npf[0]*idfx::psize, npf[2]);
      this->tempTransposedReal = IdefixArray3D<real>("FFT transpose temp2", npr[1]/idfx::psize, npr[0]*idfx::psize, npr[2]);
      this->tempT2Complex = IdefixArray3D<complex>("FFT temp2 complex", npf[0],npf[2], npf[1]);
      this->tempT2Complex2 = IdefixArray3D<complex>("FFT temp2 complex2",  npf[0],npf[2], npf[1]);
      IdefixArray3D<complex> tempComplex2 = IdefixArray3D<complex>("FFT temp complex", npf[0],npf[1],npf[2]);

     // MPI C2R Plans
      // Axis 1 transposed is axis2 for the fft library
      this->c2ciMPIPlan_axis2 = std::make_unique<PlanC2CType1D>(Kokkos::DefaultExecutionSpace(),
                                tempT2Complex, tempT2Complex2, KokkosFFT::Direction::backward, std::array<int,1>{-1});
      this->c2rMPIPlan_axis1t3 = std::make_unique<PlanC2RType2D>(Kokkos::DefaultExecutionSpace(),
                                tempTransposedComplex, tempTransposedReal, KokkosFFT::Direction::backward, std::array<int,2>{-2,-1});

      // MPI R2C plans
      this->r2cMPIPlan_axis1t3 = std::make_unique<PlanR2CType2D>(Kokkos::DefaultExecutionSpace(),
                                  tempTransposedReal, tempTransposedComplex, KokkosFFT::Direction::forward, std::array<int,2>{-2,-1});
      this->c2cfMPIPlan_axis2 = std::make_unique<PlanC2CType1D>(Kokkos::DefaultExecutionSpace(),
                                  tempT2Complex, tempT2Complex2, KokkosFFT::Direction::forward, std::array<int,1>{-1});

      this->transposeComplex = std::make_unique<Transpose<complex>>(npf);
      this->transposeReal = std::make_unique<Transpose<real>>(npr);

      #endif
    havePlan = true;
  };



void FFT::R2C_MPI(const IdefixArray3D<real> in, IdefixArray3D<complex> out, bool transpose) {
  idfx::pushRegion("FFT::R2C_MPI");

  if(transpose) {
    this->transposeReal->Apply(in,tempTransposedReal);
    idfx::pushRegion("FFT::R2C_MPI axis1t3");
    KokkosFFT::execute(*(r2cMPIPlan_axis1t3.get()), tempTransposedReal, tempTransposedComplex);
    idfx::popRegion();
  } else {
    idfx::pushRegion("FFT::R2C_MPI axis1t3");
    KokkosFFT::execute(*(r2cMPIPlan_axis1t3.get()), in, tempTransposedComplex);
    idfx::popRegion();
  }
  this->transposeComplex->Apply(tempTransposedComplex,tempComplex);
  TransposeLocal(tempComplex,tempT2Complex);
  idfx::pushRegion("FFT::R2C_MPI axis2");
  KokkosFFT::execute(*(c2cfMPIPlan_axis2.get()), tempT2Complex, tempT2Complex2);
  idfx::popRegion();
  TransposeLocal(tempT2Complex2,out);
  idfx::popRegion();
}


void FFT::C2R_MPI(const IdefixArray3D<complex> in, IdefixArray3D<real> out, bool transpose) {
  idfx::pushRegion("FFT::C2R_MPI");
  idfx::pushRegion("FFT::C2R_MPI axis2");
  TransposeLocal(in,tempT2Complex);
  KokkosFFT::execute(*(c2ciMPIPlan_axis2.get()), tempT2Complex, tempT2Complex2);
  TransposeLocal(tempT2Complex2,tempComplex);
  idfx::popRegion();
  this->transposeComplex->Apply(tempComplex,tempTransposedComplex);
  if(transpose) {
    idfx::pushRegion("FFT::C2R_MPI axis1t3");
    KokkosFFT::execute(*(c2rMPIPlan_axis1t3.get()), tempTransposedComplex, tempTransposedReal);
    idfx::popRegion();
    this->transposeReal->Apply(tempTransposedReal,out);
  } else {
    idfx::pushRegion("FFT::C2R_MPI axis1t3");
    KokkosFFT::execute(*(c2rMPIPlan_axis1t3.get()), tempTransposedComplex, out);
    idfx::popRegion();
  }
  idfx::popRegion();
}


void FFT::TestMPI() {
#ifdef WITH_MPI
  if(npr_glob[0] % idfx::psize != 0 || npr_glob[1] % idfx::psize != 0) {
    throw std::runtime_error("Global problem size must be dividible by the number of MPI processes");
  }

  // Make an array with the local problem size
  std::array<int,3> npr = npr_glob;
  npr[0] /= idfx::psize;
  std::array<int,3> npf = npr;
  npf[2] = npr[2]/2+1;
  std::array<int,3> npf_glob = npr_glob;
  npf_glob[2] = npr_glob[2]/2+1;

  Kokkos::View<real***, Kokkos::LayoutRight, Device> localReal_right("local real array", npr[0], npr[1], npr[2]);
  IdefixArray3D<real> localReal = IdefixArray3D<real>("local real array", npr[0], npr[1], npr[2]);

  Kokkos::View<real***, Kokkos::LayoutRight, Device> globalReal_right("global real array", npr_glob[0], npr_glob[1], npr_glob[2]);
  IdefixArray3D<real> globalReal = IdefixArray3D<real>("global real array", npr_glob[0], npr_glob[1], npr_glob[2]);

  // Compute a dummy real array
  if(idfx::prank == 0) {
    Kokkos::Random_XorShift64_Pool<> random_pool(12345);
    Kokkos::fill_random(globalReal_right, random_pool, 1);
  }
  // Scatter to make local arrays
  // And broadcast the global the array
  MPI_Scatter(globalReal_right.data(), npr[0]*npr[1]*npr[2],
              realMPI,
              localReal_right.data(), npr[0]*npr[1]*npr[2],
              realMPI,
              0, MPI_COMM_WORLD);

  MPI_Bcast(globalReal_right.data(), npr_glob[0]*npr_glob[1]*npr_glob[2],
            realMPI, 0, MPI_COMM_WORLD);

  // Change array Layout
  idefix_for("Reshape global",0, npr_glob[0],
                            0, npr_glob[1],
                            0, npr_glob[2],
    KOKKOS_LAMBDA(int i, int j, int k) {
      globalReal(i,j,k) = globalReal_right(i,j,k);
    });

  // Change array Layout
  idefix_for("Reshape local",0, npr[0],
                            0, npr[1],
                            0, npr[2],
    KOKKOS_LAMBDA(int i, int j, int k) {
      localReal(i,j,k) = localReal_right(i,j,k);
    });

    // Create the complex arrays
    IdefixArray3D<complex> localComplex = IdefixArray3D<complex>("local complex array", npf[0], npf[1], npf[2]);
    IdefixArray3D<complex> globalComplex = IdefixArray3D<complex>("global complex array", npf_glob[0], npf_glob[1], npf_glob[2]);

    // Compute the full serial fft
    KokkosFFT::rfftn(Kokkos::DefaultExecutionSpace(), globalReal, globalComplex);

    // Compute the parallele fft
    this->R2C_MPI(localReal, localComplex);

    // Check on a per-process that they all agree
    int offset = npf[0]*idfx::prank;

    idefix_for("Reshape local",0, npf[0],
                            0, npf[1],
                            0, npf[2],
    KOKKOS_LAMBDA(int i, int j, int k) {
      int iglob = i+offset;
      real norm = std::pow(localComplex(i,j,k).real()-globalComplex(iglob,j,k).real(),2)
                  +std::pow(localComplex(i,j,k).real()-globalComplex(iglob,j,k).real(),2);

      norm=std::sqrt(norm);

      if(norm > 1e-8) {
        Kokkos::abort("incoherent values after MPI ifft");
      }
    });

    // Compute the parallele fft
    this->C2R_MPI(localComplex, localReal);

    // Check on a per-process that they all agree
    offset = npr[0]*idfx::prank;

    idefix_for("Reshape local",0, npr[0],
                            0, npr[1],
                            0, npr[2],
    KOKKOS_LAMBDA(int i, int j, int k) {
      int iglob = i+offset;
      real norm = std::fabs(localReal(i,j,k)-globalReal(iglob,j,k));

      if(norm > 1e-8) {
        Kokkos::abort("incoherent values after MPI ifft");
      }
    });
#endif
}
