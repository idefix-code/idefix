// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2021 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef HYDRO_HYDRO_HPP_
#define HYDRO_HYDRO_HPP_

#include <string>
#include <vector>
#include "idefix.hpp"

#include "grid.hpp"
#include "hydro_defs.hpp"
#include "electroMotiveForce.hpp"
#include "viscosity.hpp"
#include "axis.hpp"
#include "fargo.hpp"
#include "hydroboundary.hpp"

// forward class declaration
class DataBlock;


class Hydro {
 public:
  Hydro();
  void Init(Input &, Grid &, DataBlock *);
  void ConvertConsToPrim();
  void ConvertPrimToCons();
  template <int> void CalcRiemannFlux(const real);
  template <int> void CalcParabolicFlux(const real);
  template <int> void AddNonIdealMHDFlux(const real);
  template <int> void CalcRightHandSide(real, real );
  void CalcCurrent();
  void AddSourceTerms(real, real );

  real GetGamma();
  real CheckDivB();
  void ResetStage();

  // Our boundary conditions
  HydroBoundary boundary;

  // Source terms
  bool haveSourceTerms{false};

  // Parabolic terms
  bool haveExplicitParabolicTerms{false};
  bool haveRKLParabolicTerms{false};

  // Current
  bool haveCurrent{false};
  bool needExplicitCurrent{false};
  bool needRKLCurrent{false};

  // Whether gravitational potential is computed
  bool haveGravPotential{false};

  // Whether a body force is present
  bool haveBodyForce{false};

  // Nonideal MHD effects coefficients
  ParabolicModuleStatus resistivityStatus, ambipolarStatus, hallStatus;

  // Whether or not we have viscosity
  ParabolicModuleStatus viscosityStatus;

  // Viscosity object
  Viscosity viscosity;

  // Whether or not we have to treat the axis
  bool haveAxis{false};
  Axis myAxis;

  // Do we use fargo-like scheme ? (orbital advection)
  bool haveFargo{false};
  Fargo fargo;

  // Rotation vector
  bool haveRotation{false};
  real OmegaZ;

  bool haveShearingBox{false};
  // Shear rate for shearing box problems
  real sbS;
  // Box width for shearing box problems
  real sbLx;


  // Enroll user-defined boundary conditions
  void EnrollUserDefBoundary(UserDefBoundaryFunc);
  void EnrollInternalBoundary(InternalBoundaryFunc);
  void EnrollEmfBoundary(EmfBoundaryFunc);
  void EnrollFluxBoundary(UserDefBoundaryFunc);

  // Enroll user-defined gravitational potential
  void EnrollGravPotential(GravPotentialFunc);

  // Enroll user-defined body force
  void EnrollBodyForce(BodyForceFunc);

  // Add some user source terms
  void EnrollUserSourceTerm(SrcTermFunc);

  // Enroll user-defined ohmic, ambipolar and Hall diffusivities
  void EnrollOhmicDiffusivity(DiffusivityFunc);
  void EnrollAmbipolarDiffusivity(DiffusivityFunc);
  void EnrollHallDiffusivity(DiffusivityFunc);

  // Enroll user-defined isothermal sound speed
  void EnrollIsoSoundSpeed(IsoSoundSpeedFunc);

  // Riemann Solvers
#if MHD == YES
  template<const int>
    void HlldMHD();
  template<const int>
    void HllMHD();
  template<const int>
    void RoeMHD();
  template<const int>
    void TvdlfMHD();
#else
  template<const int>
    void HllcHD();
  template<const int>
    void HllHD();
  template<const int>
    void RoeHD();
  template<const int>
    void TvdlfHD();
#endif
  // Extrapolate function
  template<const int DIR>
  KOKKOS_FORCEINLINE_FUNCTION void K_ExtrapolatePrimVar
      (const int, const int, const int, const IdefixArray4D<real>&,
      const IdefixArray4D<real>&, const IdefixArray1D<real>&, real[], real[]);

  // Flux functions and converter functions
  #if MHD == YES
  KOKKOS_INLINE_FUNCTION void K_Flux(real F[], real V[], real U[], real Cs2Iso,
                                   ARG_EXPAND(const int Xn, const int Xt, const int Xb),
                                   ARG_EXPAND(const int BXn, const int BXt, const int BXb));
  KOKKOS_INLINE_FUNCTION void K_ConsToPrim(real Vc[], real Uc[], real gamma_m1);
  KOKKOS_INLINE_FUNCTION void K_PrimToCons(real Uc[], real Vc[], real gamma_m1);

  #else

  KOKKOS_INLINE_FUNCTION void K_Flux(real *KOKKOS_RESTRICT F, const real *KOKKOS_RESTRICT V,
                                   const real *KOKKOS_RESTRICT U, real Cs2Iso, const int Xn);
  KOKKOS_INLINE_FUNCTION void K_ConsToPrim(real Vc[], real Uc[],real gamma_m1);
  KOKKOS_INLINE_FUNCTION void K_PrimToCons(real *KOKKOS_RESTRICT Uc, const real *KOKKOS_RESTRICT Vc,
                                         real gamma_m1);
  #endif

  // Arrays required by the Hydro object
  IdefixArray4D<real> Vc;      // Main cell-centered primitive variables index
  IdefixArray4D<real> Vs;      // Main face-centered varariables
  IdefixArray4D<real> Uc;      // Main cell-centered conservative variables
  IdefixArray4D<real> J;       // Electrical current
                               // (only defined when non-ideal MHD effects are enabled)

  // Name of the fields (used in outputs)
  std::vector<std::string> VcName;
  std::vector<std::string> VsName;

  // Storing all of the electromotive forces
  ElectroMotiveForce emf;

  // Required by time integrator
  IdefixArray4D<real> Uc0;
  IdefixArray4D<real> Vs0;
  IdefixArray3D<real> InvDt;

  IdefixArray4D<real> FluxRiemann;
  IdefixArray3D<real> dMax;    // Maximum diffusion speed



 private:
  friend class ElectroMotiveForce;
  friend class Viscosity;
  friend class Fargo;
  friend class Axis;
  friend class RKLegendre;
  friend class HydroBoundary;

  // Isothermal EOS parameters
  real isoSoundSpeed;
  HydroModuleStatus haveIsoSoundSpeed{Disabled};
  IdefixArray3D<real> isoSoundSpeedArray;
  IsoSoundSpeedFunc isoSoundSpeedFunc{NULL};

  // Adiabatic EOS parameters
  real gamma;

  Solver mySolver;

  DataBlock *data;

  // Emf boundary conditions
  bool haveEmfBoundary{false};
  EmfBoundaryFunc emfBoundaryFunc{NULL};

  // User defined gravitational potential
  GravPotentialFunc gravPotentialFunc{NULL};

  // User defined source term
  SrcTermFunc userSourceTerm{NULL};
  bool haveUserSourceTerm{false};

  real etaO, xH, xA;  // Ohmic resistivity, Hall, ambipolar (when constant)

  // Ohmic, Hall and ambipolar diffusivity (when function-defined)
  DiffusivityFunc ohmicDiffusivityFunc{NULL};
  DiffusivityFunc ambipolarDiffusivityFunc{NULL};
  DiffusivityFunc hallDiffusivityFunc{NULL};

  IdefixArray3D<real> cMax;    // Maximum propagation speed

  // Gravitational potential
  IdefixArray3D<real> phiP;

  // Body force
  IdefixArray4D<real> bodyForceVector;
  BodyForceFunc bodyForceFunc{NULL};

  // Nonideal effect diffusion coefficient (only allocated when needed)
  IdefixArray3D<real> etaOhmic;
  IdefixArray3D<real> xHall;
  IdefixArray3D<real> xAmbipolar;
};

#endif // HYDRO_HYDRO_HPP_
