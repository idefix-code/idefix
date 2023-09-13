// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************


#include <string>


#include "thermalDiffusion.hpp"
#include "dataBlock.hpp"
#include "fluid.hpp"
#include "eos.hpp"



void ThermalDiffusion::ShowConfig() {
  if(status.status==Constant) {
    idfx::cout << "Thermal Diffusion: ENABLED with constant diffusivity kappa="
                    << this->kappa << " ."<< std::endl;
  } else if (status.status==UserDefFunction) {
    idfx::cout << "Thermal Diffusion: ENABLED with user-defined diffusivity function."
                   << std::endl;
    if(!diffusivityFunc) {
      IDEFIX_ERROR("No thermal diffusion function has been enrolled");
    }
  } else {
    IDEFIX_ERROR("Unknown thermal diffusion mode");
  }
  if(status.isExplicit) {
    idfx::cout << "Thermal Diffusion: uses an explicit time integration." << std::endl;
  } else if(status.isRKL) {
    idfx::cout << "Thermal Diffusion: uses a Runge-Kutta-Legendre time integration."
                << std::endl;
  } else {
    IDEFIX_ERROR("Unknown time integrator for viscosity.");
  }
}

void ThermalDiffusion::EnrollThermalDiffusivity(DiffusivityFunc myFunc) {
  if(this->status.status != UserDefFunction) {
    IDEFIX_WARNING("Thermal diffusivity enrollment requires Hydro/ThermalDiffusion "
                 "to be set to userdef in .ini file");
  }
  this->diffusivityFunc = myFunc;
}

// This function computes the viscous flux and stores it in hydro->fluxRiemann
// (this avoids an extra array)
// Associated source terms, present in non-cartesian geometry are also computed
// and stored in this->viscSrc for later use (in calcRhs).
void ThermalDiffusion::AddDiffusiveFlux(int dir, const real t, const IdefixArray4D<real> &Flux) {
  idfx::pushRegion("ThermalDiffusion::AddDiffusiveFlux");
  IdefixArray4D<real> Vc = this->Vc;
  IdefixArray3D<real> dMax = this->dMax;
  IdefixArray3D<real> kappaArr = this->kappaArr;
  IdefixArray1D<real> dx = this->data->dx[dir];

  EquationOfState eos = *(this->eos);
  #if GEOMETRY == POLAR
    IdefixArray1D<real> x1 = this->data->x[IDIR];
  #endif
  #if GEOMETRY == SPHERICAL
    IdefixArray1D<real> rt   = this->data->rt;
    IdefixArray1D<real> dmu  = this->data->dmu;
    IdefixArray1D<real> dx2 = this->data->dx[JDIR];
  #endif

  HydroModuleStatus haveThermalDiffusion = this->status.status;

  // Compute thermal diffusion if needed
  if(haveThermalDiffusion == UserDefFunction && dir == IDIR) {
    if(diffusivityFunc) {
      idfx::pushRegion("UserDef::ThermalDiffusivityFunction");
      diffusivityFunc(*this->data, t, kappaArr);
      idfx::popRegion();
    } else {
      IDEFIX_ERROR("No user-defined thermal diffusion function has been enrolled");
    }
  }


  int ibeg, iend, jbeg, jend, kbeg, kend;
  ibeg = this->data->beg[IDIR];
  iend = this->data->end[IDIR];
  jbeg = this->data->beg[JDIR];
  jend = this->data->end[JDIR];
  kbeg = this->data->beg[KDIR];
  kend = this->data->end[KDIR];

  real kappaConstant = this->kappa;

  if(dir==IDIR) iend++;
  if(dir==JDIR) jend++;
  if(dir==KDIR) kend++;

  const int ioffset = (dir==IDIR) ? 1 : 0;
  const int joffset = (dir==JDIR) ? 1 : 0;
  const int koffset = (dir==KDIR) ? 1 : 0;

  idefix_for("ThermalDiffusionFlux",kbeg, kend, jbeg, jend, ibeg, iend,
      KOKKOS_LAMBDA (int k, int j, int i) {
        // Compute gradT
        real gradT;

        gradT = Vc(PRS,k,j,i) / Vc(RHO,k,j,i)
               - Vc(PRS,k-koffset,j-joffset,i-ioffset) / Vc(RHO,k-koffset,j-joffset,i-ioffset);

        // index along dir
        const int ig = ioffset*i + joffset*j + koffset*k;

        // dx at the interface is the averaged between the two adjacent centered dx
        real dl = HALF_F*(dx(ig-1) + dx(ig));
        #if GEOMETRY == POLAR
        if(dir==JDIR)
          dl = dl*x1(i);

        #elif GEOMETRY == SPHERICAL
          if(dir==JDIR)
            dl = dl*rt(i);
          else
            if(dir==KDIR)
              dl = dl*rt(i)*dmu(j)/dx2(j);
        #endif // GEOMETRY

        gradT = gradT/dl;

        // Compute diffusion coefficient at the interface
        real kappa;
        if(haveThermalDiffusion == UserDefFunction) {
          kappa = HALF_F*(kappaArr(k,j,i) +  kappaArr(k-koffset,j-joffset,i-ioffset));
        } else {
          kappa = kappaConstant;
        }

        // Add thermal diffusion to the flux
        Flux(ENG,k,j,i) += -kappa*gradT;

        // Compute total diffusion coefficient
        real gamma = eos.GetGamma(Vc(PRS,k,j,i),Vc(RHO,k,j,i));

        real locdmax = kappa * (gamma-ONE_F) /
                        (HALF_F * ( Vc(RHO,k,j,i) + Vc(RHO,k-koffset,j-joffset,i-ioffset)));
        dMax(k,j,i) = FMAX(dMax(k,j,i) , locdmax);
      });
  idfx::popRegion();
}
