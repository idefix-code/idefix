// ***********************************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************************

#ifndef HYDRO_HYDRO_HPP_
#define HYDRO_HYDRO_HPP_

#include <string>
#include <vector>
#include "../idefix.hpp"

#define     SMALL_PRESSURE_FIX      (1.0e-5)
#define     eps_UCT_CONTACT         (1.0e-6)

// forward class Datablock declaration
class DataBlock;

// Solver type
#if MHD == YES
enum Solver {TVDLF=1, HLL, HLLD, ROE};
#else
enum Solver {TVDLF=1, HLL, HLLC, ROE};
#endif

/*---- EMFs -----*/
#define ARITHMETIC   1
#define UCT0         2
#define UCT_CONTACT  3

// Default EMF_AVERAGE value
#ifndef EMF_AVERAGE
  #define EMF_AVERAGE     UCT_CONTACT
#endif

// Parabolic terms can have different status
enum ParabolicType {Disabled, Constant, UserDefFunction };

using UserDefBoundaryFunc = void (*) (DataBlock &, int dir, BoundarySide side,
                                      const real t);
using GravPotentialFunc = void (*) (DataBlock &, const real t, IdefixArray1D<real>&,
                                    IdefixArray1D<real>&, IdefixArray1D<real>&,
                                    IdefixArray3D<real> &);

using SrcTermFunc = void (*) (DataBlock &, const real t, const real dt);
using InternalBoundaryFunc = void (*) (DataBlock &, const real t);
using EmfBoundaryFunc = void (*) (DataBlock &, const real t);
using DiffusivityFunc = void (*) (DataBlock &, const real t, IdefixArray3D<real> &);


class Hydro {
 public:
  Hydro();
  Hydro(Input &, Grid &, DataBlock *);
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
  real GetC2iso();
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
  ParabolicType haveResistivity, haveAmbipolar, haveHall;

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



 private:
  friend class Viscosity;

  real C2Iso;
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

  // Whether or not we have viscosity
  bool haveViscosity;
  // Viscosity object
  Viscosity viscosity;

};

#endif // HYDRO_HYDRO_HPP_
