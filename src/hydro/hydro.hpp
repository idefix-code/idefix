// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef HYDRO_HYDRO_HPP_
#define HYDRO_HYDRO_HPP_

#include <string>
#include <vector>
#include "../idefix.hpp"

// forward class declaration
class DataBlock;


class Hydro {
 public:
  Hydro();
  void Init(Input &, Grid &, DataBlock *);
  void ConvertConsToPrim();
  void ConvertPrimToCons();
  void ExtrapolatePrimVar(int);
  void CalcRiemannFlux(int, const real);
  void CalcParabolicFlux(int, const real);
  void AddNonIdealMHDFlux(int, const real);
  void CalcRightHandSide(int, real, real );
  void CalcCurrent();
  void AddSourceTerms(real, real );
  void ReconstructVcField(IdefixArray4D<real> &);
  void ReconstructNormalField(int);
  void EvolveMagField(real, real);
  void CalcCornerEMF(real );
  void CalcNonidealEMF(real );
  void SetBoundary(real);
  real GetGamma();
  real CheckDivB();
  void ResetStage();

  // Source terms
  bool haveSourceTerms;

  // Parabolic terms
  bool haveParabolicTerms;

  // Current
  bool haveCurrent;
  bool needCurrent;

  // Whether gravitational potential is computed
  bool haveGravPotential;

  // Nonideal MHD effects coefficients
  HydroModuleStatus haveResistivity, haveAmbipolar, haveHall;

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


  // Whether or not we have viscosity
  bool haveViscosity;
  // Viscosity object
  Viscosity viscosity;

 private:
  friend class Viscosity;

  // Isothermal EOS parameters
  real isoSoundSpeed;
  HydroModuleStatus haveIsoSoundSpeed;
  IdefixArray3D<real> isoSoundSpeedArray;
  IsoSoundSpeedFunc isoSoundSpeedFunc;

  // Adiabatic EOS parameters
  real gamma;

  Solver mySolver;

  DataBlock *data;

  // Rotation vector
  bool haveRotation;
  real OmegaX1, OmegaX2, OmegaX3;

  bool haveShearingBox;
  // Shear rate for shearing box problems
  real sbS;
  // Box width for shearing box problems
  real sbLx;


  // User defined Boundary conditions
  UserDefBoundaryFunc userDefBoundaryFunc;
  bool haveUserDefBoundary;

  // Internal boundary function
  bool haveInternalBoundary;
  InternalBoundaryFunc internalBoundaryFunc;

  // Emf boundary conditions
  bool haveEmfBoundary;
  EmfBoundaryFunc emfBoundaryFunc;

  // User defined gravitational potential
  GravPotentialFunc gravPotentialFunc;

  // User defined source term
  SrcTermFunc userSourceTerm;
  bool haveUserSourceTerm;

  real etaO, xH, xA;  // Ohmic resistivity, Hall, ambipolar (when constant)

  // Ohmic, Hall and ambipolar diffusivity (when function-defined)
  DiffusivityFunc ohmicDiffusivityFunc;
  DiffusivityFunc ambipolarDiffusivityFunc;
  DiffusivityFunc hallDiffusivityFunc;

  IdefixArray3D<real> cMax;    // Maximum propagation speed
  IdefixArray3D<real> dMax;    // Maximum diffusion speed

  // Value extrapolated on faces
  IdefixArray4D<real> PrimL;
  IdefixArray4D<real> PrimR;
  IdefixArray4D<real> FluxRiemann;

  // Gravitational potential
  IdefixArray3D<real> phiP;

  // Nonideal effect diffusion coefficient (only allocated when needed)
  IdefixArray3D<real> etaOhmic;
  IdefixArray3D<real> xHall;
  IdefixArray3D<real> xAmbipolar;
};

#endif // HYDRO_HYDRO_HPP_
