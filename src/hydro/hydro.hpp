// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2022 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
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
#include "thermalDiffusion.hpp"
#include "axis.hpp"
#include "hydroboundary.hpp"
#include "shockFlattening.hpp"

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
  void CoarsenFlow(IdefixArray4D<real>&);
  void CoarsenMagField(IdefixArray4D<real>&);
  void CoarsenVectorPotential();
  real GetGamma();
  real CheckDivB();
  void ResetStage();
  void ShowConfig();

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

  // Nonideal MHD effects coefficients
  ParabolicModuleStatus resistivityStatus, ambipolarStatus, hallStatus;

  // Whether or not we have viscosity
  ParabolicModuleStatus viscosityStatus;

  // Whether or not we have thermal diffusion
  ParabolicModuleStatus thermalDiffusionStatus;

  // Viscosity object
  Viscosity viscosity;

  // Thermal Diffusion object
  ThermalDiffusion thermalDiffusion;

  // Whether or not we have to treat the axis
  bool haveAxis{false};
  Axis myAxis;

  // Rotation vector
  bool haveRotation{false};
  real OmegaZ;

  bool haveShearingBox{false};
  // Shear rate for shearing box problems
  real sbS;
  // Box width for shearing box problems
  real sbLx;

  // ShockFlattening
  bool haveShockFlattening{false};
  ShockFlattening shockFlattening;

  // Enroll user-defined boundary conditions
  void EnrollUserDefBoundary(UserDefBoundaryFunc);
  void EnrollInternalBoundary(InternalBoundaryFunc);
  void EnrollEmfBoundary(EmfBoundaryFunc);
  void EnrollFluxBoundary(UserDefBoundaryFunc);

  // Add some user source terms
  void EnrollUserSourceTerm(SrcTermFunc);

  // DEPRECATED gravity enrollment
  void EnrollGravPotential(GravPotentialFunc);
  void EnrollBodyForce(BodyForceFunc);

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

  // Arrays required by the Hydro object
  IdefixArray4D<real> Vc;      // Main cell-centered primitive variables index
  IdefixArray4D<real> Vs;      // Main face-centered varariables
  IdefixArray4D<real> Ve;      // Main edge-centered varariables (only when EVOLVE_VECTOR_POTENTIAL)
  IdefixArray4D<real> Uc;      // Main cell-centered conservative variables
  IdefixArray4D<real> J;       // Electrical current
                               // (only defined when non-ideal MHD effects are enabled)

  // Name of the fields (used in outputs)
  std::vector<std::string> VcName;
  std::vector<std::string> VsName;
  std::vector<std::string> VeName;

  // Storing all of the electromotive forces
  ElectroMotiveForce emf;

  // Required by time integrator
  IdefixArray4D<real> Uc0;
  IdefixArray4D<real> Vs0;
  IdefixArray4D<real> Ve0;
  IdefixArray3D<real> InvDt;

  IdefixArray4D<real> FluxRiemann;
  IdefixArray3D<real> dMax;    // Maximum diffusion speed



 private:
  friend class ElectroMotiveForce;
  friend class Viscosity;
  friend class ThermalDiffusion;
  friend class Fargo;
  friend class Axis;
  friend class RKLegendre;
  friend class HydroBoundary;
  friend class ShockFlattening;

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

  // User defined source term
  SrcTermFunc userSourceTerm{NULL};
  bool haveUserSourceTerm{false};

  real etaO, xH, xA;  // Ohmic resistivity, Hall, ambipolar (when constant)

  // Ohmic, Hall and ambipolar diffusivity (when function-defined)
  DiffusivityFunc ohmicDiffusivityFunc{NULL};
  DiffusivityFunc ambipolarDiffusivityFunc{NULL};
  DiffusivityFunc hallDiffusivityFunc{NULL};

  IdefixArray3D<real> cMax;    // Maximum propagation speed

  // Nonideal effect diffusion coefficient (only allocated when needed)
  IdefixArray3D<real> etaOhmic;
  IdefixArray3D<real> xHall;
  IdefixArray3D<real> xAmbipolar;
};

#endif // HYDRO_HYDRO_HPP_
