// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef GRAVITY_GRAVITY_HPP_
#define GRAVITY_GRAVITY_HPP_

#include "idefix.hpp"
#include "input.hpp"
#include "selfGravity.hpp"
#include "planetarySystem.hpp"

class DataBlock;

// Enrolled functions signature
using GravPotentialFunc = void (*) (DataBlock &, const real t, IdefixArray1D<real>&,
                                    IdefixArray1D<real>&, IdefixArray1D<real>&,
                                    IdefixArray3D<real> &);

using BodyForceFunc = void (*) (DataBlock &, const real t, IdefixArray4D<real>&);


class Gravity {
 public:
  Gravity(Input&, DataBlock*);
  void ComputeGravity(int );           ///< compute gravitational field at current time t

  void EnrollPotential(GravPotentialFunc);
  void EnrollBodyForce(BodyForceFunc);

  void ResetPotential();            ///< fill the potential with zeros.

  void AddCentralMassPotential();   ///< Àdd the potential due to a centrall mass

  void ShowConfig();                ///< Show the gravity configuration
  bool havePotential{false};        ///< Whether a gravitational potential is present
                                        ///< in which case, (at least) one of the following is true
  bool haveUserDefPotential{false};     ///< Whether a potential is defined by user
  bool haveCentralMassPotential{false}; ///< Whether a potential is due to the central mass
  bool havePlanetsPotential{false};     ///< Whether a potential is due to planet(s)
  bool haveSelfGravityPotential{false}; ///< Whether a potential is defined through self-gravity

  bool haveBodyForce{false};            ///< Whether a body force (=acceleration) is present

  // Gravitational potential
  IdefixArray3D<real> phiP;

  // Bodyforce
  IdefixArray4D<real> bodyForceVector;

  // Self gravity
  SelfGravity selfGravity;

  // JM : moved in public class to handle changing centralMass during computation
  real centralMass{1.0};                    ///< central mass parameter when central mass potential
                                            ///< is enabled

  // Gravitational constant G
  real gravCst{1.0};

  // Whether we should skip gravity computation every n steps
  int skipGravity{1};


 private:
  friend class PlanetarySystem;
  bool haveInitialisedPotential{false};     ///< whether a potential has already been initialised
  bool haveInitialisedBodyForce{false};     ///< whether a body force has already been initialised
  bool haveInitialisedSelfGravity{false};   ///< whether self-gravity has already been initialised

  DataBlock *data;

  // User defined gravitational potential
  GravPotentialFunc gravPotentialFunc{NULL};

  // Body force
  BodyForceFunc bodyForceFunc{NULL};

  #ifdef DEBUG_GRAVITY
  // Used to get fields usefull for debugging
  std::ofstream potTotFile;
  #endif
};

#endif // GRAVITY_GRAVITY_HPP_
