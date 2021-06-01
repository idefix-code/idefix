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

// forward class declaration
class DataBlock;


class Hydro {
 public:
  Hydro();
  void Init(Input &, Grid &, DataBlock *);
  void ConvertConsToPrim();
  void ConvertPrimToCons();
  void CalcRiemannFlux(int, const real);
  void CalcParabolicFlux(int, const real);
  void AddNonIdealMHDFlux(int, const real);
  void CalcRightHandSide(int, real, real );
  void CalcCurrent();
  void AddSourceTerms(real, real );
  void ReconstructVcField(IdefixArray4D<real> &);
  void ReconstructNormalField(int);


  void SetBoundary(real);
  void EnforceBoundaryDir(real, int);
  real GetGamma();
  real CheckDivB();
  void ResetStage();

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

  // Enroll user-defined boundary conditions
  void EnrollUserDefBoundary(UserDefBoundaryFunc);
  void EnrollInternalBoundary(InternalBoundaryFunc);
  void EnrollEmfBoundary(EmfBoundaryFunc);

  // Enroll user-defined gravitational potential
  void EnrollGravPotential(GravPotentialFunc);

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
  template<const int DIR, ARG_EXPAND(const int Xn, const int Xt, const int Xb),
         ARG_EXPAND(const int BXn, const int BXt, const int BXb)>
    void HlldMHD();
  template<const int DIR, ARG_EXPAND(const int Xn, const int Xt, const int Xb),
         ARG_EXPAND(const int BXn, const int BXt, const int BXb)>
    void HllMHD();
  template<const int DIR, ARG_EXPAND(const int Xn, const int Xt, const int Xb),
         ARG_EXPAND(const int BXn, const int BXt, const int BXb)>
    void RoeMHD();
  template<const int DIR, ARG_EXPAND(const int Xn, const int Xt, const int Xb),
                        ARG_EXPAND(const int BXn, const int BXt, const int BXb)>
    void TvdlfMHD();
#else
  template<const int DIR, const int Xn, const int Xt, const int Xb>
    void HllcHD();
  template<const int DIR, const int Xn, const int Xt, const int Xb>
    void HllHD();
  template<const int DIR, const int Xn, const int Xt, const int Xb>
    void RoeHD();
  template<const int DIR, const int Xn, const int Xt, const int Xb>
    void TvdlfHD();
#endif
  // Extrapolate function
  template<const int DIR>
  KOKKOS_FORCEINLINE_FUNCTION void K_ExtrapolatePrimVar
      (const int, const int, const int, const IdefixArray4D<real>&,
      const IdefixArray4D<real>&, real[], real[]);

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

  // Isothermal EOS parameters
  real isoSoundSpeed;
  HydroModuleStatus haveIsoSoundSpeed{Disabled};
  IdefixArray3D<real> isoSoundSpeedArray;
  IsoSoundSpeedFunc isoSoundSpeedFunc{NULL};

  // Adiabatic EOS parameters
  real gamma;

  Solver mySolver;

  DataBlock *data;

  // Rotation vector
  bool haveRotation{false};
  real OmegaZ;

  bool haveShearingBox{false};
  // Shear rate for shearing box problems
  real sbS;
  // Box width for shearing box problems
  real sbLx;


  // User defined Boundary conditions
  UserDefBoundaryFunc userDefBoundaryFunc{NULL};
  bool haveUserDefBoundary{false};

  // Internal boundary function
  bool haveInternalBoundary{false};
  InternalBoundaryFunc internalBoundaryFunc{NULL};

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

  // Nonideal effect diffusion coefficient (only allocated when needed)
  IdefixArray3D<real> etaOhmic;
  IdefixArray3D<real> xHall;
  IdefixArray3D<real> xAmbipolar;
};

#endif // HYDRO_HYDRO_HPP_
