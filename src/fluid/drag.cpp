// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************
#include "drag.hpp"
#include "physics.hpp"
void Drag::AddDragForce(const real dt) {
  idfx::pushRegion("Drag::AddDragForce");

  auto UcGas = this->UcGas;
  auto VcGas = this->VcGas;
  auto UcDust = this->UcDust;
  auto VcDust = this->VcDust;
  auto InvDt = this->InvDt;

  const Type type = this->type;
  real dragCoeff = this->dragCoeff;
  bool feedback = this->feedback;

  EquationOfState eos = *(this->eos);

  auto userGammai = this->gammai;
  if(type == Type::Userdef) {
    if(userDrag != NULL) {
      idfx::pushRegion("Drag::UserDrag");
      userDrag(data, dragCoeff, userGammai);
      idfx::popRegion();
    } else {
      IDEFIX_ERROR("No User-defined drag function has been enrolled");
    }
  }
  // Compute a drag force fd = - gamma*rhod*rhog*(vd-vg)
  // Where gamma is computed according to the choice of drag type
  idefix_for("DragForce",0,data->np_tot[KDIR],0,data->np_tot[JDIR],0,data->np_tot[IDIR],
    KOKKOS_LAMBDA (int k, int j, int i) {
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
          cs = std::sqrt(eos.GetGamma(VcGas(PRS,k,j,i),VcGas(RHO,k,j,i)
                         *VcGas(PRS,k,j,i)/VcGas(RHO,k,j,i)));
        #else
          cs = eos.GetWaveSpeed(k,j,i);
        #endif
        gamma = cs/dragCoeff;
      } else if(type == Type::Userdef) {
        gamma = userGammai(k,j,i);
      }

      real dp = dt * gamma * VcDust(RHO,k,j,i) * VcGas(RHO,k,j,i);
      for(int n = MX1 ; n < MX1+COMPONENTS ; n++) {
        real dv = VcDust(n,k,j,i) - VcGas(n,k,j,i);
        UcDust(n,k,j,i) -= dp*dv;
        if(feedback) UcGas(n,k,j,i) += dp*dv;
        #if HAVE_ENERGY == 1
          // We add back the energy dissipated for the dust which is not accounted for
          // (since there is no energy equation for dust grains)

          // TODO(GL): this should be disabled in the case of a true multifluid system where
          // both fluids have a proper energy equation
          UcGas(ENG,k,j,i) += dp*dv*VcDust(n,k,j,i);
        #endif
      }
      // Cfl constraint
      real idt = gamma*VcGas(RHO,k,j,i);
      if(feedback) idt += gamma*VcDust(RHO,k,j,i);
      InvDt(k,j,i) += idt;
    });
  idfx::popRegion();
}

void Drag::ShowConfig() {
  idfx::cout << "Drag: Using ";
  switch(type) {
    case Type::Gamma:
      idfx::cout << "constant gamma";
      break;
    case Type::Tau:
      idfx::cout << "constant stopping time";
      break;
    case Type::Size:
      idfx::cout << "constant dust size";
      break;
    case Type::Userdef:
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
  if(type != Type::Userdef) {
    IDEFIX_ERROR("User-defined drag function requires drag entry to be set to \"userdef\"");
  }
  this->userDrag = func;
}
