#include "idefix.hpp"
#include "coagulation.hpp"
#include "npy.hpp"

Coagulation *CoagulationGlob;

Coagulation::Coagulation(DataBlock &data, int kernel, real eps, real coagCFL) {
    this->kernel = kernel;
    this->eps = eps;
    this->coagCFL = coagCFL;

    //==============================================//
    //==============================================//
    // In this part, we load the coagulation tables //
    //==============================================//
    //==============================================//
    const int nSpecies = data.dust.size();
    //===== Load the mass grid =====// 
    std::string filename_massgrid = "data/test-massgrid-kernel-"+std::to_string(kernel)+".npy";
    std::vector<uint64_t> shape_massgrid;
    bool fortran_order_massgrid;
    std::vector<double> massgrid;
    npy::LoadArrayFromNumpy(filename_massgrid, shape_massgrid, fortran_order_massgrid, massgrid);
      //===== verification =====//
    if(massgrid.size() != static_cast<size_t>(nSpecies+1)) {
      IDEFIX_ERROR("Mass grid has an unexpected size");
    }
      //========================//

    this->binWidth = IdefixArray1D<real>("binWidth", nSpecies);
    IdefixHostArray1D<real> binWidthHost = Kokkos::create_mirror_view(this->binWidth);
    for(int n = 0; n <nSpecies; n++) {
      binWidthHost(n) = massgrid[n+1] - massgrid[n];
    }
    Kokkos::deep_copy(this->binWidth, binWidthHost);

    //===== Loading the Flux Table =====//
    std::string filename_coagTabFlux = "data/test-coagTabFlux-kernel-"+std::to_string(kernel)+".npy";
    std::vector<uint64_t> shape_coagTabFlux;
    bool fortran_order_coagTabFlux;
    std::vector<double> coagTabFlux;
    npy::LoadArrayFromNumpy(filename_coagTabFlux, shape_coagTabFlux, fortran_order_coagTabFlux, coagTabFlux);
      //===== verification =====//
    const size_t expectedSize = static_cast<size_t>(nSpecies)*static_cast<size_t>(nSpecies)*static_cast<size_t>(nSpecies);
    if(coagTabFlux.size() != expectedSize) {
      IDEFIX_ERROR("Coagulation flux table has an unexpected size");
    }
      //========================//

    this->coagTabFlux = IdefixArray3D<real>("coagTabFlux", nSpecies, nSpecies, nSpecies);
    IdefixHostArray3D<real> coagTabFluxHost = Kokkos::create_mirror_view(this->coagTabFlux);
    for(int i=0; i<nSpecies; i++) {
      for(int j=0; j<nSpecies; j++) {
        for(int k=0; k<nSpecies; k++) {
          int index = i*nSpecies*nSpecies + j*nSpecies + k;
          coagTabFluxHost(i,j,k) = coagTabFlux[index];
        }
      }
    }
    Kokkos::deep_copy(this->coagTabFlux, coagTabFluxHost);
}

void CoagSourceTerm(DataBlock &dataRef, const real t, const real dtin) {
  //==============================================//
  /*
  Coagulation source term : compute the coagulation product to finally upload the dust densities. 
  This function is implemented via operator splitting : coagulation is integrated in a step.
  Hence, enroll the function with `data.EnrollUserStepFirst(&CoagSourceTerm);`.
  This version does not preserve bins momentum yet.
  This version does not work on GPU yet.
  */
  //==============================================//
  DataBlock *data = &dataRef;
  data->PrimToCons();
  IdefixArray1D<real> h = CoagulationGlob->binWidth;
  IdefixArray3D<real> coagTabFlux = CoagulationGlob->coagTabFlux;
  const int nSpecies = data->dust.size();
  const int test_size = 20;  // Because we can't give nSpecies in the KOKKOS_LAMBDA

  const real eps = CoagulationGlob->eps;
  const real coagCFL = CoagulationGlob->coagCFL;
  const bool physicalKernel = (CoagulationGlob->kernel==3);

  idefix_for("CoagSourceTerm",0,data->np_tot[KDIR],0,data->np_tot[JDIR],0,data->np_tot[IDIR],
            KOKKOS_LAMBDA (int k, int j, int i) {
              real rho[test_size];
              real g0[test_size];
              real vx1[test_size];
              real vx2[test_size];
              real vx3[test_size];
              real deltaV[test_size][test_size];
              real Flux[test_size+1];

              // Initial State
              for(int n=0; n<nSpecies; n++) {
                rho[n] = data->dust[n]->Uc(RHO,k,j,i);
                vx1[n] = data->dust[n]->Uc(MX1,k,j,i)/rho[n];
                vx2[n] = data->dust[n]->Uc(MX2,k,j,i)/rho[n];
                vx3[n] = data->dust[n]->Uc(MX3,k,j,i)/rho[n];

                /*
                // To have a constant drift velocity that does not evolve in sub-cycling
                // WARNING: dirty here! check that it corresponds to Compute_test_velocity in utils.py
                real vx = 1e-2 + (1e-1 - 1e-2)*n/(nSpecies - 1);
                real vy = 1.3*vx;
                real vz = 1.7*vx;
                vx1[n] = vx;
                vx2[n] = vy;
                vx3[n] = vz;
                */
              }
              
              // Relative velocity
              if(physicalKernel) {
                for(int lp=0; lp<nSpecies; lp++) {
                  for(int l=0; l<nSpecies; l++) {
                    real dvx = vx1[lp] - vx1[l];
                    real dvy = vx2[lp] - vx2[l];
                    real dvz = vx3[lp] - vx3[l];
                    deltaV[lp][l] = sqrt(dvx*dvx + dvy*dvy + dvz*dvz);
                  }
                }
              }

              //===== COAGULATION SUBCYCLING =====//
              real time = 0.0;
              while(time<dtin) {
                Flux[0] = 0.0;  // No mass flux crossing the left boundary of the distribution
                for(int n=0; n<nSpecies; n++) {
                  g0[n] = rho[n]/h(n);
                }
                //===== FLUX CALCULATION =====//
                for(int n=0; n<nSpecies; n++) {
                  real flux = 0.0;
                  for(int lp=0; lp <= n; lp++) {
                    real g0_lp = g0[lp];
                    for(int l=0; l<nSpecies; l++) {
                      real g0_l = g0[l];
                      real delta_v = physicalKernel ? deltaV[lp][l] : 1.0;
                      flux += g0_lp*g0_l*coagTabFlux(n,lp,l)*delta_v;
                    }
                  }
                  Flux[n+1] = flux;
                }
                //===== CFL CALCULATION =====//
                real dt_sub = std::numeric_limits<double>::max();
                for(int n=0; n<nSpecies; n++) {
                  real dflux = Flux[n+1] - Flux[n];
                  if(std::abs(dflux)>0.0 && rho[n]>eps*h(n)) {
                    dt_sub = std::min(dt_sub, std::abs(rho[n]/dflux));
                  }
                }
                real dt_coag = dt_sub*coagCFL;
                dt_coag = std::min(dt_coag, dtin-time);  // Last substep
                if(dt_coag<=0.0) {
                  IDEFIX_ERROR("Coagulation time step is non-positive");
                }
                //===== EULER UPDATE =====//
                for(int n=0; n<nSpecies; n++) {
                  rho[n] = std::max(rho[n] - dt_coag*(Flux[n+1] - Flux[n]), eps*h(n));
                }
                time += dt_coag;
              }
              //===== UPDATE THE FINAL DENSITIES =====//
              for(int n=0; n<nSpecies; n++) {
                data->dust[n]->Uc(RHO,k,j,i) = rho[n];

                // To keep a constant velocity, just for the test with kernel = 3
                data->dust[n]->Uc(MX1,k,j,i) = rho[n]*vx1[n];
                data->dust[n]->Uc(MX2,k,j,i) = rho[n]*vx2[n];
                data->dust[n]->Uc(MX3,k,j,i) = rho[n]*vx3[n];
              }
            }
  );
  data->ConsToPrim();
}