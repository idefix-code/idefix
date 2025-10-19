// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************
#ifndef FLUID_DRAG_HPP_
#define FLUID_DRAG_HPP_

#include <string>
#include "idefix.hpp"
#include "input.hpp"
#include "fluid_defs.hpp"
#include "eos.hpp"

using UserDefDragFunc = void (*) (DataBlock *, int n, real beta, IdefixArray3D<real> &gammai);

/* A class that holds the kind of drag law we intend to use */
class GammaDrag {
 public:
  enum class Type{Gamma, Tau, Size, Userdef};
  GammaDrag() = default;
  GammaDrag(Input &, std::string BlockName, int instanceNumber, DataBlock *data);
  void EnrollUserDrag(UserDefDragFunc);
  void RefreshUserDrag(DataBlock *);

  KOKKOS_INLINE_FUNCTION real GetGamma(const int k, const int j, const int i) const {
    real gamma;  // The drag coefficient
    if(type ==  Type::Gamma) {
      gamma = dragCoeff;

    } else if(type == Type::Tau) {
      // In this case, the coefficient is the stopping time (assumed constant)
      gamma = 1/(dragCoeff*VcGas(RHO,k,j,i));
    } else if(type == Type::Size) {
      real cs;
      // Assume a fixed size, hence for both Epstein or Stokes, gamma~1/rho_g/cs
      // Get the sound speed
      #if HAVE_ENERGY == 1
        cs = std::sqrt(  eos.GetGamma(VcGas(PRS,k,j,i),VcGas(RHO,k,j,i))
                        *VcGas(PRS,k,j,i)/VcGas(RHO,k,j,i));
      #else
        cs = eos.GetWaveSpeed(k,j,i);
      #endif
      gamma = cs/dragCoeff;
    } else if(type == Type::Userdef) {
      gamma = gammai(k,j,i);
    }
    return gamma;
  }

  Type type;
  real dragCoeff;
  EquationOfState eos;
  IdefixArray3D<real> gammai;
  IdefixArray4D<real> VcGas;

  int instanceNumber;
  UserDefDragFunc userDrag{NULL};
};


class Drag {
 public:
  // Different types of implementation for the drag force.
  template <typename Phys>
  Drag(Input &, Fluid<Phys> *);
  void ShowConfig();                    // print configuration
  void AddDragForce(const real);

  //////////////////////////
  // Implicit functions
  void AddImplicitBackReaction(const real, IdefixArray3D<real>);  // Add the back reaction
  void NormalizeImplicitBackReaction(const real);  // Normalize the implicit back reaction
  void AddImplicitFluidMomentum(const real);  // Add the implicit drag force on dust grains
  /////////////////////////

  void EnrollUserDrag(UserDefDragFunc);   // User defined drag function enrollment
  bool IsImplicit() const { return implicit; }  // Check if the drag is implicit

  IdefixArray4D<real> UcDust;  // Dust conservative quantities
  IdefixArray4D<real> UcGas;  // Gas conservative quantities
  IdefixArray4D<real> VcDust;  // Gas primitive quantities
  IdefixArray4D<real> VcGas;  // Gas primitive quantities
  IdefixArray3D<real> InvDt;  // The InvDt of current dust specie
  IdefixArray3D<real> implicitFactor; // The prefactor used by the implicit timestepping

  GammaDrag gammaDrag;  // The drag law

 private:
  DataBlock* data;
  bool feedback{false};
  bool implicit{false};
  int instanceNumber;
};


#include "fluid.hpp"

template<typename Phys>
Drag::Drag(Input &input, Fluid<Phys> *hydroin):
                      UcDust{hydroin->Uc},
                      UcGas{hydroin->data->hydro->Uc},
                      VcDust{hydroin->Vc},
                      VcGas{hydroin->data->hydro->Vc},
                      InvDt{hydroin->InvDt} {
  idfx::pushRegion("Drag::Drag");
  // Save the parent hydro object

  this->data = hydroin->data;

  // Check in which block we should fetch our information
  std::string blockName;
  if(Phys::dust) {
    blockName = "Dust";
  } else {
    IDEFIX_ERROR("Drag is currently implemented only for dusty fuids");
  }

  if(input.CheckEntry(blockName,"drag")>=0) {
    this->instanceNumber = hydroin->instanceNumber;
    this->gammaDrag = GammaDrag(input,blockName,instanceNumber,data);

    // Feedback is true by default, but can be switched off.
    this->feedback = input.GetOrSet<bool>(blockName,"drag_feedback",0,true);
    this->implicit = input.GetOrSet<bool>(blockName,"drag_implicit",0,false);

    if(implicit && instanceNumber == 0) {
      this->implicitFactor = IdefixArray3D<real>("ImplicitFactor",
                                  data->np_tot[KDIR],
                                  data->np_tot[JDIR],
                                  data->np_tot[IDIR]);
    }

  } else {
    IDEFIX_ERROR("A [Drag] block is required in your input file to define the drag force.");
  }


  idfx::popRegion();
}
#endif // FLUID_DRAG_HPP_
