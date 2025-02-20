// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************
#include "drag.hpp"
#include <string>
#include "physics.hpp"

void Drag::AddDragForce(const real dt) {
  idfx::pushRegion("Drag::AddDragForce");

  auto UcGas = this->UcGas;
  auto VcGas = this->VcGas;
  auto UcDust = this->UcDust;
  auto VcDust = this->VcDust;
  auto InvDt = this->InvDt;

  bool feedback = this->feedback;

  if(implicit) {
    IDEFIX_ERROR("Add DragForce should not be called when drag is implicit");
  }

  auto gammaDrag = this->gammaDrag;
  gammaDrag.RefreshUserDrag(data);

  // Compute a drag force fd = - gamma*rhod*rhog*(vd-vg)
  // Where gamma is computed according to the choice of drag type
  idefix_for("DragForce",0,data->np_tot[KDIR],0,data->np_tot[JDIR],0,data->np_tot[IDIR],
    KOKKOS_LAMBDA (int k, int j, int i) {
      real gamma = gammaDrag.GetGamma(k,j,i);  // The drag coefficient

      real dp = dt * gamma * VcDust(RHO,k,j,i) * VcGas(RHO,k,j,i);
      for(int n = MX1 ; n < MX1+COMPONENTS ; n++) {
        real dv = VcDust(n,k,j,i) - VcGas(n,k,j,i);
        UcDust(n,k,j,i) -= dp*dv;
        if(feedback) {
          UcGas(n,k,j,i) += dp*dv;
          #if HAVE_ENERGY == 1
            // We add back the energy dissipated for the dust which is not accounted for
            // (since there is no energy equation for dust grains)

            // TODO(GL): this should be disabled in the case of a true multifluid system where
            // both fluids have a proper energy equation
            UcGas(ENG,k,j,i) += dp*dv*VcDust(n,k,j,i);
          #endif
        } // feedback
      }
      // Cfl constraint
      real idt = gamma*VcGas(RHO,k,j,i);
      if(feedback) idt += gamma*VcDust(RHO,k,j,i);
      InvDt(k,j,i) += idt;
    });
  idfx::popRegion();
}
//
// $ p_g^{(n+1)}=p_g^{(n)}+\sum_i\frac{\rho_g\gamma_i dt}{1+\rho_g \gamma_i dt}p_i^{(n)} $
// We accumulate in the array "prefactor"
// $ \sum_i\frac{\rho_i\gamma_i dt}{1+\rho_g\gamma_i dt}
//

void Drag::AddImplicitBackReaction(const real dt, IdefixArray3D<real> preFactor) {
  if(!feedback) {
    // no feedback, no need for anything
    return;
  }
  idfx::pushRegion("AddImplicitFluidMomentum");

  auto UcGas = this->UcGas;
  auto VcGas = this->VcGas;
  auto VcDust = this->VcDust;
  auto UcDust = this->UcDust;

  bool feedback = this->feedback;

  bool isFirst = this->instanceNumber == 0;


  if(!implicit) {
    IDEFIX_ERROR("AddImplicitGasMomentum should not be called when drag is explicit");
  }

  auto gammaDrag = this->gammaDrag;
  gammaDrag.RefreshUserDrag(data);

  // Compute a drag force fd = - gamma*rhod*rhog*(vd-vg)
  // Where gamma is computed according to the choice of drag type
  idefix_for("DragForce",0,data->np_tot[KDIR],0,data->np_tot[JDIR],0,data->np_tot[IDIR],
    KOKKOS_LAMBDA (int k, int j, int i) {
      real gamma = gammaDrag.GetGamma(k,j,i);  // The drag coefficient


      const real factor = UcDust(RHO,k,j,i)*gamma*dt/(1+UcGas(RHO,k,j,i)*gamma*dt);
      if(isFirst) {
        preFactor(k,j,i) = factor;
      } else {
        preFactor(k,j,i) += factor;
      }

      for(int n = MX1 ; n < MX1+COMPONENTS ; n++) {
        UcGas(n,k,j,i) +=  dt * gamma * UcGas(RHO,k,j,i) * UcDust(n,k,j,i) /
                                                                   (1 + UcGas(RHO,k,j,i)*dt*gamma);
      }
    });
  idfx::popRegion();
}

void Drag::NormalizeImplicitBackReaction(const real dt) {
  if(!feedback) {
    // no feedback, no need for anything
    return;
  }
  idfx::pushRegion("AddImplicitFluidMomentum");

  auto UcGas = this->UcGas;
  auto preFactor = this->implicitFactor;

  if(!implicit) {
    IDEFIX_ERROR("AddImplicitGasMomentum should not be called when drag is explicit");
  }

  // Compute a drag force fd = - gamma*rhod*rhog*(vd-vg)
  // Where gamma is computed according to the choice of drag type
  idefix_for("DragForce",0,data->np_tot[KDIR],0,data->np_tot[JDIR],0,data->np_tot[IDIR],
    KOKKOS_LAMBDA (int k, int j, int i) {
      const real factor = 1+preFactor(k,j,i);
      for(int n = MX1 ; n < MX1+COMPONENTS ; n++) {
        UcGas(n,k,j,i) /= factor;
      }
    });
  idfx::popRegion();
}

void Drag::AddImplicitFluidMomentum(const real dt) {
  idfx::pushRegion("AddImplicitFluidMomentum");

  auto UcGas = this->UcGas;
  auto VcGas = this->VcGas;
  auto UcDust = this->UcDust;
  auto VcDust = this->VcDust;
  auto InvDt = this->InvDt;

  bool feedback = this->feedback;

  if(!implicit) {
    IDEFIX_ERROR("AddImplicitGasMomentum should not be called when drag is explicit");
  }

  auto gammaDrag = this->gammaDrag;
  if(!feedback) gammaDrag.RefreshUserDrag(data);

  // Compute a drag force fd = - gamma*rhod*rhog*(vd-vg)
  // Where gamma is computed according to the choice of drag type
  idefix_for("DragForce",0,data->np_tot[KDIR],0,data->np_tot[JDIR],0,data->np_tot[IDIR],
    KOKKOS_LAMBDA (int k, int j, int i) {
      real gamma = gammaDrag.GetGamma(k,j,i);  // The drag coefficient

      for(int n = MX1 ; n < MX1+COMPONENTS ; n++) {
        real oldUc = UcDust(n,k,j,i);
        UcDust(n,k,j,i) = (oldUc + dt * gamma * UcDust(RHO,k,j,i) * UcGas(n,k,j,i)) /
                          (1 + UcGas(RHO,k,j,i)*dt*gamma);

        #if HAVE_ENERGY == 1
          real dp = UcDust(n,k,j,i) - oldUc;

          // We add back the energy dissipated for the dust which is not accounted for
          // (since there is no energy equation for dust grains)

          // TODO(GL): this should be disabled in the case of a true multifluid system where
          // both fluids have a proper energy equation
          if(feedback) UcGas(ENG,k,j,i) -= dp*VcDust(n,k,j,i);
        #endif
      }
    });
  idfx::popRegion();
}

void Drag::ShowConfig() {
  idfx::cout << "Drag: Using ";
  if(implicit) {
    idfx::cout << " IMPLICIT ";
  } else {
    idfx::cout << " EXPLICIT ";
  }
  switch(gammaDrag.type) {
    case GammaDrag::Type::Gamma:
      idfx::cout << "constant gamma";
      break;
    case GammaDrag::Type::Tau:
      idfx::cout << "constant stopping time";
      break;
    case GammaDrag::Type::Size:
      idfx::cout << "constant dust size";
      break;
    case GammaDrag::Type::Userdef:
      idfx::cout << "user-defined";
      break;
  }

  idfx::cout << " drag law";
  if(feedback) {
    idfx::cout << " with feedback." << std::endl;
  } else {
    idfx::cout << " without feedback." << std::endl;
  }
}

void Drag::EnrollUserDrag(UserDefDragFunc func) {
  gammaDrag.EnrollUserDrag(func);
}

////////////////////////////////////////////
// GammaDrag function definitions
////////////////////////////////////////////

void GammaDrag::RefreshUserDrag(DataBlock *data) {
  if(type == Type::Userdef) {
    if(userDrag != NULL) {
      idfx::pushRegion("GammaDrag::UserDrag");
        userDrag(data, instanceNumber, dragCoeff, gammai);
      idfx::popRegion();
    } else {
      IDEFIX_ERROR("No User-defined drag function has been enrolled");
    }
  }
}

void GammaDrag::EnrollUserDrag(UserDefDragFunc func) {
  if(type != Type::Userdef) {
    IDEFIX_ERROR("User-defined drag function requires drag entry to be set to \"userdef\"");
  }
  this->userDrag = func;
}

GammaDrag::GammaDrag(Input &input, std::string BlockName, int instanceNumber, DataBlock *data) {
  if(input.CheckEntry(BlockName,"drag")>=0) {
    std::string dragType = input.Get<std::string>(BlockName,"drag",0);
    if(dragType.compare("gamma") == 0) {
      this->type = Type::Gamma;
    } else if(dragType.compare("tau") == 0) {
      this->type = Type::Tau;
      this->VcGas = data->hydro->Vc;
    } else if(dragType.compare("size") == 0) {
      this->type = Type::Size;
      this->eos = *(data->hydro->eos.get());
      this->VcGas = data->hydro->Vc;
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
  }
  // Fetch the drag coefficient for the current specie.
  this->dragCoeff = input.Get<real>(BlockName,"drag",instanceNumber+1);
  this->instanceNumber = instanceNumber;
}
