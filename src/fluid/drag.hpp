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

using UserDefDragFunc = void (*) (DataBlock *, real beta, IdefixArray3D<real> &gammai);


class Drag {
 public:
  enum class Type{Gamma, Tau, Size, Userdef};
  // Different types of implementation for the drag force.
  template <typename Phys>
  Drag(Input &, Fluid<Phys> *);
  void ShowConfig();                    // print configuration
  void AddDragForce(const real);
  void EnrollUserDrag(UserDefDragFunc);   // User defined drag function enrollment

  IdefixArray4D<real> UcDust;  // Dust conservative quantities
  IdefixArray4D<real> UcGas;  // Gas conservative quantities
  IdefixArray4D<real> VcDust;  // Gas primitive quantities
  IdefixArray4D<real> VcGas;  // Gas primitive quantities
  IdefixArray3D<real> InvDt;  // The InvDt of current dust specie
  IdefixArray3D<real> gammai; // the drag coefficient (only used for user-defined dust grains)
  Type type;

 private:
  DataBlock* data;
  real dragCoeff;
  bool feedback{false};

  UserDefDragFunc userDrag{NULL};

  // Sound speed computation
  EquationOfState *eos;
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
  // We use the EOS of the basic fluid instance, as there is no EOS for dust!
  this->eos = hydroin->data->hydro->eos.get();

  // Check in which block we should fetch our information
  std::string BlockName;
  if(Phys::dust) {
    BlockName = "Dust";
  } else {
    IDEFIX_ERROR("Drag is currently implemented only for dusty fuids");
  }

  if(input.CheckEntry(BlockName,"drag")>=0) {
    std::string dragType = input.Get<std::string>(BlockName,"drag",0);
    if(dragType.compare("gamma") == 0) {
      this->type = Type::Gamma;
    } else if(dragType.compare("tau") == 0) {
      this->type = Type::Tau;
    } else if(dragType.compare("size") == 0) {
      this->type = Type::Size;
    } else if(dragType.compare("userdef") == 0) {
      this->type = Type::Userdef;
      this->gammai = IdefixArray3D<real>("UserDrag",
                                  data->np_tot[KDIR],
                                  data->np_tot[JDIR],
                                  data->np_tot[IDIR]);
    } else {
      std::stringstream msg;
      msg << "Unknown drag type \"" <<  dragType
          << "\" in your input file." << std::endl
          << "Allowed values are: gamma, tau, epstein, stokes, userdef." << std::endl;

      IDEFIX_ERROR(msg);
    }
    // Fetch the drag coefficient for the current specie.
    const int n = hydroin->instanceNumber;
    this->dragCoeff = input.Get<real>(BlockName,"drag",n+1);

    // Feedback is true by default, but can be switched off.
    this->feedback = input.GetOrSet<bool>(BlockName,"drag_feedback",0,true);
  } else {
    IDEFIX_ERROR("A [Drag] block is required in your input file to define the drag force.");
  }


  idfx::popRegion();
}
#endif // FLUID_DRAG_HPP_
