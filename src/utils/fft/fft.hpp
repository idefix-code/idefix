// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef FFT_HPP_
#define FFT_HPP_

#include <memory>
#include <type_traits>
#include <KokkosFFT.hpp>
#include "arrays.hpp"
#include "transpose.hpp"

using PlanR2CType = KokkosFFT::Plan<Kokkos::DefaultExecutionSpace, IdefixArray3D<real>, IdefixArray3D<complex>,3>;
using PlanC2RType = KokkosFFT::Plan<Kokkos::DefaultExecutionSpace, IdefixArray3D<complex>, IdefixArray3D<real>,3>;

using PlanC2CType1D = KokkosFFT::Plan<Kokkos::DefaultExecutionSpace, IdefixArray3D<complex>, IdefixArray3D<complex>,1>;
using PlanC2RType1D = KokkosFFT::Plan<Kokkos::DefaultExecutionSpace, IdefixArray3D<complex>, IdefixArray3D<real>,1>;
using PlanR2CType1D = KokkosFFT::Plan<Kokkos::DefaultExecutionSpace, IdefixArray3D<real>, IdefixArray3D<complex>,1>;

using PlanR2CType2D = KokkosFFT::Plan<Kokkos::DefaultExecutionSpace, IdefixArray3D<real>, IdefixArray3D<complex>,2>;
using PlanC2RType2D = KokkosFFT::Plan<Kokkos::DefaultExecutionSpace, IdefixArray3D<complex>, IdefixArray3D<real>,2>;

class FFT {
 public:
  FFT();
  FFT(std::array<int,3> npr_glob, std::array<int,3> npf_glob);

  template<typename InView, typename OutView>
  void R2C(const InView &in, const OutView &out, bool transpose=false);

  template<typename InView, typename OutView>
  void C2R(const InView &in, const OutView &out, bool transpose=false);

  template<typename InView, typename OutView>
  void R2C_Host(const InView &in, const OutView &out);

  template<typename InView, typename OutView>
  void C2R_Host(const InView &in, const OutView &out);

  template<typename ViewIn, typename ViewOut>
  void TransposeLocal(const ViewIn &in, const ViewOut &out);

  void R2C_MPI(const IdefixArray3D<real> in, IdefixArray3D<complex> out, bool transpose=false);
  void C2R_MPI(const IdefixArray3D<complex> in, IdefixArray3D<real> out, bool transpose=false);
  void TestMPI();

 private:
  template<typename View>
  using view_t = std::decay_t<View>;

  template<typename View>
  static constexpr bool is_view_v = Kokkos::is_view<view_t<View>>::value;

  template<typename View>
  static constexpr bool is_rank3_v = is_view_v<View> && (view_t<View>::rank == 3);

  template<typename View, typename Scalar>
  static constexpr bool has_scalar_v =
      is_view_v<View> &&
      std::is_same_v<typename view_t<View>::non_const_value_type, Scalar>;

  template<typename View>
  static constexpr bool device_accessible_v =
      is_view_v<View> &&
      Kokkos::SpaceAccessibility<
        Kokkos::DefaultExecutionSpace,
        typename view_t<View>::memory_space>::accessible;

  template<typename View>
  static constexpr bool host_accessible_v =
      is_view_v<View> &&
      Kokkos::SpaceAccessibility<
        Kokkos::DefaultHostExecutionSpace,
        typename view_t<View>::memory_space>::accessible;

 public:
  std::array<int,3> npr_glob, npf_glob;
  std::array<int,3> npr, npf, npr_t;
  bool havePlan{false};

  IdefixArray3D<real> tempReal;
  IdefixArray3D<complex> tempComplex;

#ifdef WITH_MPI
  IdefixArray3D<complex> tempTransposedComplex, tempTransposedComplex2;
  IdefixArray3D<real> tempTransposedReal;
  IdefixArray3D<complex> tempT2Complex, tempT2Complex2;

  std::unique_ptr<Transpose<complex>> transposeComplex;
  std::unique_ptr<Transpose<real>> transposeReal;

  std::unique_ptr<PlanC2CType1D> c2ciMPIPlan_axis2;
  std::unique_ptr<PlanC2RType2D> c2rMPIPlan_axis1t3;
  std::unique_ptr<PlanR2CType2D> r2cMPIPlan_axis1t3;
  std::unique_ptr<PlanC2CType1D> c2cfMPIPlan_axis2;
#endif

  std::unique_ptr<PlanR2CType> r2cPlan;
  std::unique_ptr<PlanC2RType> c2rPlan;
};

template<typename InView, typename OutView>
void FFT::R2C(const InView &in, const OutView &out, bool transpose) {
  static_assert(is_rank3_v<InView>, "FFT::R2C: input must be a rank-3 Kokkos::View");
  static_assert(is_rank3_v<OutView>, "FFT::R2C: output must be a rank-3 Kokkos::View");
  static_assert(has_scalar_v<InView, real>, "FFT::R2C: input scalar type must be real");
  static_assert(has_scalar_v<OutView, complex>, "FFT::R2C: output scalar type must be complex");
  static_assert(device_accessible_v<InView>, "FFT::R2C: input view must be accessible from DefaultExecutionSpace");
  static_assert(device_accessible_v<OutView>, "FFT::R2C: output view must be accessible from DefaultExecutionSpace");

  idfx::pushRegion("FFT::R2C");

#ifdef WITH_MPI
  Kokkos::deep_copy(tempReal, in);
  R2C_MPI(tempReal, tempComplex, transpose);
  Kokkos::deep_copy(out, tempComplex);
#else
  Kokkos::deep_copy(tempReal, in);
  KokkosFFT::execute(*(r2cPlan.get()), tempReal, tempComplex);
  Kokkos::deep_copy(out, tempComplex);
#endif

  idfx::popRegion();
}

template<typename InView, typename OutView>
void FFT::C2R(const InView &in, const OutView &out, bool transpose) {
  static_assert(is_rank3_v<InView>, "FFT::C2R: input must be a rank-3 Kokkos::View");
  static_assert(is_rank3_v<OutView>, "FFT::C2R: output must be a rank-3 Kokkos::View");
  static_assert(has_scalar_v<InView, complex>, "FFT::C2R: input scalar type must be complex");
  static_assert(has_scalar_v<OutView, real>, "FFT::C2R: output scalar type must be real");
  static_assert(device_accessible_v<InView>, "FFT::C2R: input view must be accessible from DefaultExecutionSpace");
  static_assert(device_accessible_v<OutView>, "FFT::C2R: output view must be accessible from DefaultExecutionSpace");

  idfx::pushRegion("FFT::C2R");

#ifdef WITH_MPI
  Kokkos::deep_copy(tempComplex, in);
  C2R_MPI(tempComplex, tempReal, transpose);
  Kokkos::deep_copy(out, tempReal);
#else
  Kokkos::deep_copy(tempComplex, in);
  KokkosFFT::execute(*(c2rPlan.get()), tempComplex, tempReal);
  Kokkos::deep_copy(out, tempReal);
#endif

  idfx::popRegion();
}

template<typename InView, typename OutView>
void FFT::R2C_Host(const InView &in, const OutView &out) {
  static_assert(is_rank3_v<InView>, "FFT::R2C_Host: input must be a rank-3 Kokkos::View");
  static_assert(is_rank3_v<OutView>, "FFT::R2C_Host: output must be a rank-3 Kokkos::View");
  static_assert(has_scalar_v<InView, real>, "FFT::R2C_Host: input scalar type must be real");
  static_assert(has_scalar_v<OutView, complex>, "FFT::R2C_Host: output scalar type must be complex");
  static_assert(host_accessible_v<InView>, "FFT::R2C_Host: input view must be accessible from host");
  static_assert(host_accessible_v<OutView>, "FFT::R2C_Host: output view must be accessible from host");

  IdefixArray3D<real> inDev = Kokkos::create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), in);
  IdefixArray3D<complex> outDev = Kokkos::create_mirror_view(Kokkos::DefaultExecutionSpace(), out);
  R2C(inDev, outDev, true);
  Kokkos::deep_copy(out, outDev);
}

template<typename InView, typename OutView>
void FFT::C2R_Host(const InView &in, const OutView &out) {
  static_assert(is_rank3_v<InView>, "FFT::C2R_Host: input must be a rank-3 Kokkos::View");
  static_assert(is_rank3_v<OutView>, "FFT::C2R_Host: output must be a rank-3 Kokkos::View");
  static_assert(has_scalar_v<InView, complex>, "FFT::C2R_Host: input scalar type must be complex");
  static_assert(has_scalar_v<OutView, real>, "FFT::C2R_Host: output scalar type must be real");
  static_assert(host_accessible_v<InView>, "FFT::C2R_Host: input view must be accessible from host");
  static_assert(host_accessible_v<OutView>, "FFT::C2R_Host: output view must be accessible from host");

  IdefixArray3D<complex> inDev = Kokkos::create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), in);
  IdefixArray3D<real> outDev = Kokkos::create_mirror_view(Kokkos::DefaultExecutionSpace(), out);
  C2R(inDev, outDev, true);
  Kokkos::deep_copy(out, outDev);
}

template<typename ViewIn, typename ViewOut>
void FFT::TransposeLocal(const ViewIn &in, const ViewOut &out) {
  static_assert(is_rank3_v<ViewIn>, "FFT::TransposeLocal: input must be a rank-3 Kokkos::View");
  static_assert(is_rank3_v<ViewOut>, "FFT::TransposeLocal: output must be a rank-3 Kokkos::View");
  static_assert(device_accessible_v<ViewIn>, "FFT::TransposeLocal: input view must be device accessible");
  static_assert(device_accessible_v<ViewOut>, "FFT::TransposeLocal: output view must be device accessible");

  idefix_for("TransposeLocal", 0, in.extent(0),
                               0, in.extent(1),
                               0, in.extent(2),
    KOKKOS_LAMBDA(int i, int j, int k) {
      out(i,k,j) = in(i,j,k);
    });
}

#endif // FFT_HPP_