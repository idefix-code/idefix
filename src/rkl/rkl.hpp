// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef RKL_RKL_HPP_
#define RKL_RKL_HPP_

#include <string>
#include <vector>

#include "idefix.hpp"
#include "input.hpp"
#include "dataBlock.hpp"
#include "viscosity.hpp"
#include "bragViscosity.hpp"
#ifdef WITH_MPI
#include "mpi.hpp"
#endif

// Functors being used by RKL
template <typename Phys>
struct RKLegendre_ResetStageFunctor;

template<typename Phys>
class RKLegendre {
 public:
  RKLegendre(Input &, Fluid<Phys>*);
  void Cycle();
  void ResetStage();
  void ResetFlux();
  void EvolveStage(real);
  template <int> void CalcParabolicRHS(real);
  void ComputeDt();
  void ShowConfig();
  void Copy(IdefixArray4D<real>&, IdefixArray4D<real>&);

  IdefixArray4D<real> dU;      // variation of main cell-centered conservative variables
  IdefixArray4D<real> dU0;      // dU of the first stage
  IdefixArray4D<real> Uc0;      // Uc at initial stage
  IdefixArray4D<real> Uc1;      // Uc of the previous stage, Uc1 = Uc(stage-1)

  IdefixArray4D<real> dB;      // Variation of cell-centered magnetic variables
  IdefixArray4D<real> dB0;     // dB of the first stage
  IdefixArray4D<real> Vs0;     // Vs of initial stage
  IdefixArray4D<real> Vs1;     // Vs of previous stage

  #ifdef EVOLVE_VECTOR_POTENTIAL
  IdefixArray4D<real> dA;      // Variation of edge-centered vector potential
  IdefixArray4D<real> dA0;     // dA of the first stage
  IdefixArray4D<real> Ve0;     // Ve of initial stage
  IdefixArray4D<real> Ve1;     // Ve of previous stage
  #endif

  IdefixArray1D<int> varList;  // List of variables which should be evolved
  int nvarRKL{0};               // # of active variables

  real dt, cfl_rkl, rmax_par;
  int stage{0};

 private:
  friend struct RKLegendre_ResetStageFunctor<Phys>;
  void SetBoundaries(real);        // Enforce boundary conditions on the variables solved by RKL

  DataBlock *data;
  Fluid<Phys> *hydro;

#ifdef WITH_MPI
  Mpi mpi;                      // RKL-specific MPI layer
#endif

  bool haveVs{false};           // Whether we have (and need to compute) face-centered variables
  bool haveVc{false};           // Whether we need to compute cell-centered variables
  void AddVariable(int, std::vector<int> & );

  bool checkNan{false};         // whether we should look for Nans when RKL is running

 private:
  template<int> void LoopDir(real);   // Dimensional loop
};

#include "fluid.hpp"
#include "calcParabolicFlux.hpp"

#ifndef RKL_ORDER
  #define RKL_ORDER       2
#endif


template<typename Phys>
void RKLegendre<Phys>::AddVariable(int var, std::vector<int> &varListHost ) {
  bool haveit{false};

  // Check whether we have this variable in the list
  for(int i = 0 ; i < varListHost.size() ; i++) {
    if(varListHost[i] == var) haveit=true;
  }

  // We don't have it, then add it to the list
  if(!haveit) {
    varListHost.push_back(var);
  }
}

// Copy just the variables required by the RK scheme
template<typename Phys>
void RKLegendre<Phys>::Copy(IdefixArray4D<real> &out, IdefixArray4D<real> &in) {
  IdefixArray1D<int> vars = this->varList;

  idefix_for("RKL_Copy",
             0, nvarRKL,
             0, data->np_tot[KDIR],
             0, data->np_tot[JDIR],
             0, data->np_tot[IDIR],
             KOKKOS_LAMBDA(int n, int k, int j, int i) {
               const int var = vars(n);
               out(var,k,j,i) = in(var,k,j,i);
             });
}

template<typename Phys>
RKLegendre<Phys>::RKLegendre(Input &input, Fluid<Phys>* hydroin) {
  idfx::pushRegion("RKLegendre::Init");

  // Save the datablock to which we are attached from now on
  this->data = hydroin->data;
  this->hydro = hydroin;

  cfl_rkl = input.GetOrSet<real> ("RKL","cfl",0, 0.5);
  rmax_par = input.GetOrSet<real> ("RKL","rmax_par",0, 100.0);

  // By default check nans in debug mode
  #ifdef DEBUG
  this->checkNan = true;
  #endif

  this->checkNan = input.GetOrSet<bool>("RKL","check_nan",0, this->checkNan);

  // Make a list of variables

  std::vector<int> varListHost;
  // Create a list of variables
  // Viscosity
  if(hydro->viscosityStatus.isRKL) {
    haveVc = true;
    EXPAND( AddVariable(MX1, varListHost);   ,
            AddVariable(MX2, varListHost);   ,
            AddVariable(MX3, varListHost);   )

    #if HAVE_ENERGY
      AddVariable(ENG, varListHost);
    #endif
  }
  // BragViscosity
  if(hydro->bragViscosityStatus.isRKL) {
    haveVc = true;
    EXPAND( AddVariable(MX1, varListHost);   ,
            AddVariable(MX2, varListHost);   ,
            AddVariable(MX3, varListHost);   )
    #if HAVE_ENERGY
      AddVariable(ENG, varListHost);
    #endif
  }

  // Thermal diffusion
  #if HAVE_ENERGY
    if(hydro->thermalDiffusionStatus.isRKL) {
      haveVc = true;
      AddVariable(ENG, varListHost);
    }
    // Braginskii Thermal diffusion
    if(hydro->bragThermalDiffusionStatus.isRKL) {
      haveVc = true;
      AddVariable(ENG, varListHost);
    }
  #endif
  // Ambipolar diffusion
  if(hydro->ambipolarStatus.isRKL || hydro->resistivityStatus.isRKL) {
    #if COMPONENTS == 3 && DIMENSIONS < 3
      haveVc = true;
      AddVariable(BX3, varListHost);
    #endif
    #if COMPONENTS >= 2 && DIMENSIONS < 2
      haveVc = true;
      AddVariable(BX2, varListHost);
    #endif
    #if HAVE_ENERGY
      haveVc = true;
      AddVariable(ENG, varListHost);
    #endif
    haveVs = true;
  }

  // Copy the list on the device
  varList = idfx::ConvertVectorToIdefixArray(varListHost);
  nvarRKL = varListHost.size();

  #ifdef WITH_MPI
    mpi.Init(data->mygrid, varListHost, data->nghost.data(), data->np_int.data(), haveVs);
  #endif


  // Variable allocation

  dU = IdefixArray4D<real>("RKL_dU", NVAR,
                           data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
  dU0 = IdefixArray4D<real>("RKL_dU0", NVAR,
                           data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
  Uc0 = IdefixArray4D<real>("RKL_Uc0", NVAR,
                           data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
  Uc1 = IdefixArray4D<real>("RKL_Uc1", NVAR,
                           data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);

  if(haveVs) {
    #ifdef EVOLVE_VECTOR_POTENTIAL
      dA = IdefixArray4D<real>("RKL_dA", AX3e+1,
                      data->np_tot[KDIR]+KOFFSET,
                      data->np_tot[JDIR]+JOFFSET,
                      data->np_tot[IDIR]+IOFFSET);
      dA0 = IdefixArray4D<real>("RKL_dA0", AX3e+1,
                        data->np_tot[KDIR]+KOFFSET,
                        data->np_tot[JDIR]+JOFFSET,
                        data->np_tot[IDIR]+IOFFSET);
      Ve0 = IdefixArray4D<real>("RKL_Ve0", AX3e+1,
                        data->np_tot[KDIR]+KOFFSET,
                        data->np_tot[JDIR]+JOFFSET,
                        data->np_tot[IDIR]+IOFFSET);
      Ve1 = IdefixArray4D<real>("RKL_Ve1", AX3e+1,
                        data->np_tot[KDIR]+KOFFSET,
                        data->np_tot[JDIR]+JOFFSET,
                        data->np_tot[IDIR]+IOFFSET);
    #else
      dB = IdefixArray4D<real>("RKL_dB", DIMENSIONS,
                        data->np_tot[KDIR]+KOFFSET,
                        data->np_tot[JDIR]+JOFFSET,
                        data->np_tot[IDIR]+IOFFSET);
      dB0 = IdefixArray4D<real>("RKL_dB0", DIMENSIONS,
                        data->np_tot[KDIR]+KOFFSET,
                        data->np_tot[JDIR]+JOFFSET,
                        data->np_tot[IDIR]+IOFFSET);
      Vs0 = IdefixArray4D<real>("RKL_Vs0", DIMENSIONS,
                        data->np_tot[KDIR]+KOFFSET,
                        data->np_tot[JDIR]+JOFFSET,
                        data->np_tot[IDIR]+IOFFSET);
      Vs1 = IdefixArray4D<real>("RKL_Vs1", DIMENSIONS,
                        data->np_tot[KDIR]+KOFFSET,
                        data->np_tot[JDIR]+JOFFSET,
                        data->np_tot[IDIR]+IOFFSET);
    #endif
  }

  idfx::popRegion();
}

template<typename Phys>
void RKLegendre<Phys>::ShowConfig() {
  #if RKL_ORDER == 1
    idfx::cout << "RKLegendre: 1st order scheme ENABLED." << std::endl;
  #elif RKL_ORDER == 2
    idfx::cout << "RKLegendre: 2nd order scheme ENABLED." << std::endl;
  #else
    IDEFIX_ERROR("Unknown RKL scheme order");
  #endif
  idfx::cout << "RKLegendre: RKL cfl set to " << cfl_rkl <<  "." << std::endl;
  idfx::cout << "RKLegendre: maximum ratio hyperbolic/parabolic timestep "
             << rmax_par <<  "." << std::endl;
  if(haveVc) {
     idfx::cout << "RKLegendre: will evolve cell-centered fields Vc." << std::endl;
  }
  if(haveVs) {
     idfx::cout << "RKLegendre: will evolve face-centered fields Vs." << std::endl;
  }
  if(checkNan) {
    idfx::cout << "RKLegendre: will check consistency of solution in the integrator (slow!)."
               << std::endl;
  }
}

template<typename Phys>
void RKLegendre<Phys>::Cycle() {
  idfx::pushRegion("RKLegendre::Cycle");

  IdefixArray4D<real> dU = this->dU;
  IdefixArray4D<real> dU0 = this->dU0;
  IdefixArray4D<real> Uc = hydro->Uc;
  IdefixArray4D<real> Uc0 = this->Uc0;
  IdefixArray4D<real> Uc1 = this->Uc1;

  IdefixArray4D<real> dB = this->dB;
  IdefixArray4D<real> dB0 = this->dB0;
  IdefixArray4D<real> Vs = hydro->Vs;
  IdefixArray4D<real> Vs0 = this->Vs0;
  IdefixArray4D<real> Vs1 = this->Vs1;

  #ifdef EVOLVE_VECTOR_POTENTIAL
  IdefixArray4D<real> dA = this->dA;
  IdefixArray4D<real> dA0 = this->dA0;
  IdefixArray4D<real> Ve = hydro->Ve;
  IdefixArray4D<real> Ve0 = this->Ve0;
  IdefixArray4D<real> Ve1 = this->Ve1;
  #endif

  IdefixArray1D<int> varList = this->varList;
  real time = data->t;

  real dt_hyp = data->dt;

  // Tell the datablock that we're performing the RKL cycle
  data->rklCycle = true;

  // first RKL stage
  stage = 1;

  // Apply Boundary conditions on the full set of variables
  hydro->boundary->SetBoundaries(time);

  // Convert current state into conservative variable
  hydro->ConvertPrimToCons();

  // Coarsen the conservative variables if needed
  if(data->haveGridCoarsening) {
    data->Coarsen();
  }

  // Store the result in Uc0
  Copy(Uc0,Uc);
  if(haveVs) {
    #ifdef EVOLVE_VECTOR_POTENTIAL
      Kokkos::deep_copy(Ve0,Ve);
    #else
      Kokkos::deep_copy(Vs0,Vs);
    #endif
  }

  // evolve RKL stage
  EvolveStage(time);

  ComputeDt();

  Copy(dU0,dU);

  if(haveVs) {
    #ifdef EVOLVE_VECTOR_POTENTIAL
      Kokkos::deep_copy(dA0,dA);
    #else
      Kokkos::deep_copy(dB0,dB);
    #endif
  }

  // Compute number of RKL steps
  real nrkl;
  real scrh =  dt_hyp/dt;
#if RKL_ORDER == 1
  // Solution of quadratic Eq.
  // 2*dt_hyp/dt_exp = s^2 + s
  nrkl = 4.0*scrh / (1.0 + std::sqrt(1.0 + 8.0*scrh));
#elif RKL_ORDER == 2
  // Solution of quadratic Eq.
  // 4*dt_hyp/dt_exp = s^2 + s - 2
  nrkl = 4.0*(1.0 + 2.0*scrh)
          / (1.0 + sqrt(9.0 + 16.0*scrh));
#else
  //#error Invalid RKL_ORDER
#endif
  int rklstages = 1 + floor(nrkl);

  // Compute coefficients
  real w1, mu_tilde_j;
#if RKL_ORDER == 1
  w1 = 2.0/(rklstages*rklstages + rklstages);
  mu_tilde_j = w1;
#elif RKL_ORDER == 2
  real b_j, b_jm1, b_jm2, a_jm1;
  w1 = 4.0/(rklstages*rklstages + rklstages - 2.0);
  mu_tilde_j = w1/3.0;

  b_j = b_jm1 = b_jm2 = 1.0/3.0;
  a_jm1 = 1.0 - b_jm1;
#endif

#if RKL_ORDER == 1
  time = data->t + 0.5*dt_hyp*(stage*stage+stage)*w1;
#elif RKL_ORDER == 2
  time = data->t + 0.25*dt_hyp*(stage*stage+stage-2)*w1;
#endif
  if(haveVc) {
    idefix_for("RKL_Cycle_InitUc1",
              0, nvarRKL,
              data->beg[KDIR],data->end[KDIR],
              data->beg[JDIR],data->end[JDIR],
              data->beg[IDIR],data->end[IDIR],
      KOKKOS_LAMBDA (int n, int k, int j, int i) {
        int nv = varList(n);
        Uc1(nv,k,j,i) = Uc(nv,k,j,i);
        Uc(nv,k,j,i) = Uc1(nv,k,j,i) + mu_tilde_j*dt_hyp*dU0(nv,k,j,i);
      }
    );
  }
  if(haveVs) {
    #ifdef EVOLVE_VECTOR_POTENTIAL
      idefix_for("RKL_Cycle_InitVe1",
              0, AX3e+1,
              data->beg[KDIR],data->end[KDIR]+KOFFSET,
              data->beg[JDIR],data->end[JDIR]+JOFFSET,
              data->beg[IDIR],data->end[IDIR]+IOFFSET,
      KOKKOS_LAMBDA (int n, int k, int j, int i) {
        Ve1(n,k,j,i) = Ve(n,k,j,i);
        Ve(n,k,j,i) = Ve1(n,k,j,i) + mu_tilde_j*dt_hyp*dA0(n,k,j,i);
      });
      hydro->emf->ComputeMagFieldFromA(Ve,Vs);
    #else
      idefix_for("RKL_Cycle_InitVs1",
              0, DIMENSIONS,
              data->beg[KDIR],data->end[KDIR]+KOFFSET,
              data->beg[JDIR],data->end[JDIR]+JOFFSET,
              data->beg[IDIR],data->end[IDIR]+IOFFSET,
      KOKKOS_LAMBDA (int n, int k, int j, int i) {
        Vs1(n,k,j,i) = Vs(n,k,j,i);
        Vs(n,k,j,i) = Vs1(n,k,j,i) + mu_tilde_j*dt_hyp*dB0(n,k,j,i);
      });
    #endif
    hydro->boundary->ReconstructVcField(Uc);
  }

  // Coarsen conservative variables once they have been evolved
  if(data->haveGridCoarsening) {
    data->Coarsen();
  }

  // Convert current state into primitive variable
  hydro->ConvertConsToPrim();
  if(checkNan) {
    if(data->CheckNan()>0) {
      throw std::runtime_error(std::string("Nan found during RKL stage 1"));
    }
  }

  real mu_j, nu_j, gamma_j;
  // subStages loop
  for(stage=2; stage <= rklstages ; stage++) {
    //idfx::cout << "RKL: looping stages" << std::endl;
    // compute RKL coefficients
#if RKL_ORDER == 1
    mu_j       = (2.0*stage -1.0)/stage;
    mu_tilde_j = w1*mu_j;
    nu_j       = -(stage -1.0)/stage;
#elif RKL_ORDER == 2
    mu_j       = (2.0*stage -1.0)/stage * b_j/b_jm1;
    mu_tilde_j = w1*mu_j;
    gamma_j    = -a_jm1*mu_tilde_j;
    nu_j       = -(stage -1.0)*b_j/(stage*b_jm2);

    b_jm2 = b_jm1;
    b_jm1 = b_j;
    a_jm1 = 1.0 - b_jm1;
    b_j   = 0.5*(stage*stage+3.0*stage)/(stage*stage+3.0*stage+2.0);
#endif

    // Apply Boundary conditions
    this->SetBoundaries(time);

    // evolve RKL stage
    EvolveStage(time);
    if(haveVc) {
      // update Uc
      idefix_for("RKL_Cycle_UpdateUc",
              0, nvarRKL,
              data->beg[KDIR],data->end[KDIR],
              data->beg[JDIR],data->end[JDIR],
              data->beg[IDIR],data->end[IDIR],
        KOKKOS_LAMBDA (int n, int k, int j, int i) {
          const int nv = varList(n);
          real Y = mu_j*Uc(nv,k,j,i) + nu_j*Uc1(nv,k,j,i);
          Uc1(nv,k,j,i) = Uc(nv,k,j,i);
  #if RKL_ORDER == 1
          Uc(nv,k,j,i) = Y + dt_hyp*mu_tilde_j*dU(nv,k,j,i);
  #elif RKL_ORDER == 2
          Uc(nv,k,j,i) = Y + (1.0 - mu_j - nu_j)*Uc0(nv,k,j,i)
                                  + dt_hyp*mu_tilde_j*dU(nv,k,j,i)
                                  + gamma_j*dt_hyp*dU0(nv,k,j,i);
  #endif
          });
    }
    if(haveVs) {
      #ifdef EVOLVE_VECTOR_POTENTIAL
        // update Ve
        idefix_for("RKL_Cycle_UpdateVe",
                0, AX3e+1,
                data->beg[KDIR],data->end[KDIR]+KOFFSET,
                data->beg[JDIR],data->end[JDIR]+JOFFSET,
                data->beg[IDIR],data->end[IDIR]+IOFFSET,
          KOKKOS_LAMBDA (int n, int k, int j, int i) {
            real Y = mu_j*Ve(n,k,j,i) + nu_j*Ve1(n,k,j,i);
            Ve1(n,k,j,i) = Ve(n,k,j,i);
            #if RKL_ORDER == 1
              Ve(n,k,j,i) = Y + dt_hyp*mu_tilde_j*dA(n,k,j,i);
            #elif RKL_ORDER == 2
              Ve(n,k,j,i) = Y + (1.0 - mu_j - nu_j)*Ve0(n,k,j,i)
                                    + dt_hyp*mu_tilde_j*dA(n,k,j,i)
                                    + gamma_j*dt_hyp*dA0(n,k,j,i);
            #endif
          });
        hydro->emf->ComputeMagFieldFromA(Ve,Vs);
      #else
        // update Vs
        idefix_for("RKL_Cycle_UpdateVs",
                0, DIMENSIONS,
                data->beg[KDIR],data->end[KDIR]+KOFFSET,
                data->beg[JDIR],data->end[JDIR]+JOFFSET,
                data->beg[IDIR],data->end[IDIR]+IOFFSET,
          KOKKOS_LAMBDA (int n, int k, int j, int i) {
            real Y = mu_j*Vs(n,k,j,i) + nu_j*Vs1(n,k,j,i);
            Vs1(n,k,j,i) = Vs(n,k,j,i);
            #if RKL_ORDER == 1
              Vs(n,k,j,i) = Y + dt_hyp*mu_tilde_j*dB(n,k,j,i);
            #elif RKL_ORDER == 2
              Vs(n,k,j,i) = Y + (1.0 - mu_j - nu_j)*Vs0(n,k,j,i)
                                    + dt_hyp*mu_tilde_j*dB(n,k,j,i)
                                    + gamma_j*dt_hyp*dB0(n,k,j,i);
            #endif
          });
      #endif  // EVOLVE_VECTOR_POTENTIAL

      hydro->boundary->ReconstructVcField(Uc);
    }

    // Coarsen the flow if needed
    if(data->haveGridCoarsening) {
      data->Coarsen();
    }
    // Convert current state into primitive variable
    hydro->ConvertConsToPrim();

    if(checkNan) {
      if(data->CheckNan()>0) {
        throw std::runtime_error(std::string("Nan found during RKL stage ")+std::to_string(stage));
      }
    }

    // increment time
#if RKL_ORDER == 1
    time = data->t + 0.5*dt_hyp*(stage*stage+stage)*w1;
#elif RKL_ORDER == 2
    time = data->t + 0.25*dt_hyp*(stage*stage+stage-2)*w1;
#endif
  }

  // Tell the datablock that we're done
  data->rklCycle = false;
  idfx::popRegion();
}



template<typename Phys>
void RKLegendre<Phys>::ResetFlux() {
  idfx::pushRegion("RKLegendre::ResetFlux");
  IdefixArray4D<real> Flux = hydro->FluxRiemann;
  IdefixArray1D<int> vars = this->varList;
  idefix_for("RKL_ResetFlux",
             0,nvarRKL,
             0,data->np_tot[KDIR],
             0,data->np_tot[JDIR],
             0,data->np_tot[IDIR],
    KOKKOS_LAMBDA (int n, int k, int j, int i) {
      const int nv = vars(n);
      Flux(nv,k,j,i) = ZERO_F;
    }
  );
  idfx::popRegion();
}

template<typename Phys>
struct RKLegendre_ResetStageFunctor {
  explicit RKLegendre_ResetStageFunctor(RKLegendre<Phys> *rkl) {
    dU = rkl->dU;
    Flux = rkl->hydro->FluxRiemann;
    vars = rkl->varList;
    stage = rkl->stage;
    nvar = rkl->nvarRKL;
    haveVc = rkl->haveVc || (rkl->stage ==1 );
    haveVs = rkl->haveVs;
    invDt = rkl->hydro->InvDt;
    if constexpr(Phys::mhd) {
      #ifdef EVOLVE_VECTOR_POTENTIAL
        dA = rkl->dA;
      #else
        dB = rkl->dB;
      #endif
      ex = rkl->hydro->emf->ex;
      ey = rkl->hydro->emf->ey;
      ez = rkl->hydro->emf->ez;
    }
  }

  IdefixArray4D<real> dU;
  IdefixArray4D<real> Flux;
  IdefixArray1D<int> vars;
  IdefixArray4D<real> dA, dB;
  IdefixArray3D<real> ex,ey,ez;
  IdefixArray3D<real> invDt;
  int stage, nvar;
  bool haveVs, haveVc;

  KOKKOS_INLINE_FUNCTION void operator() (const int k, const int j,  const int i) const {
    if(haveVc) {
      for(int n = 0 ; n < nvar ; n++) {
        const int nv = vars(n);
        dU(nv,k,j,i) = ZERO_F;
      }
    }
    if(stage == 1)   invDt(k,j,i) = ZERO_F;
    if constexpr(Phys::mhd) {
      if(haveVs) {
        #ifdef EVOLVE_VECTOR_POTENTIAL
          for(int n=0; n < AX3e+1; n++) {
            dA(n,k,j,i) = ZERO_F;
          }
        #else
          for(int n=0; n < DIMENSIONS; n++) {
            dB(n,k,j,i) = ZERO_F;
          }
        #endif
        D_EXPAND( ez(k,j,i) = 0.0;    ,
                                      ,
                  ex(k,j,i) = 0.0;
                  ey(k,j,i) = 0.0;    )
      }
    }
  }
};

template<typename Phys>
void RKLegendre<Phys>::ResetStage() {
  idfx::pushRegion("RKLegendre::ResetStage");

  auto func = RKLegendre_ResetStageFunctor(this);


  idefix_for("RKL_ResetStage",
             0,data->np_tot[KDIR],
             0,data->np_tot[JDIR],
             0,data->np_tot[IDIR],
             func);

  idfx::popRegion();
}


template<typename Phys>
void RKLegendre<Phys>::ComputeDt() {
  idfx::pushRegion("RKLegendre::ComputeDt");

  IdefixArray3D<real> invDt = hydro->InvDt;

  real newinvdt = ZERO_F;
  Kokkos::parallel_reduce("RKL_Timestep_reduction",
    Kokkos::MDRangePolicy<Kokkos::Rank<3, Kokkos::Iterate::Right, Kokkos::Iterate::Right>>
    ({0,0,0},{data->end[KDIR],data->end[JDIR],data->end[IDIR]}),
    KOKKOS_LAMBDA (int k, int j, int i, real &invdt) {
      invdt = std::fmax(invDt(k,j,i), invdt);
    },
    Kokkos::Max<real>(newinvdt)
  );

#ifdef WITH_MPI
  if(idfx::psize>1) {
    MPI_SAFE_CALL(MPI_Allreduce(MPI_IN_PLACE, &newinvdt, 1, realMPI, MPI_MAX, MPI_COMM_WORLD));
  }
#endif

  dt = 1.0/newinvdt;
  dt = (cfl_rkl*dt)/2.0; // parabolic time step

  idfx::popRegion();
}

template<typename Phys>
template<int dir>
void RKLegendre<Phys>::LoopDir(real t) {
    ResetFlux();

    // CalcParabolicFlux
    hydro->template CalcParabolicFlux<dir>(t);

    // Calc Right Hand Side
    CalcParabolicRHS<dir>(t);

    // Recursive: do next dimension
    if constexpr(dir+1<DIMENSIONS) {
      LoopDir<dir+1>(t);
    }
}


template<typename Phys>
void RKLegendre<Phys>::EvolveStage(real t) {
  idfx::pushRegion("RKLegendre::EvolveStage");

  ResetStage();

  if(haveVs && hydro->needRKLCurrent) hydro->CalcCurrent();

  // Loop on dimensions for the parabolic fluxes and RHS, starting from IDIR
  if(haveVc || stage == 1) LoopDir<IDIR>(t);

  if(haveVs) {
    hydro->emf->CalcNonidealEMF(t);
    hydro->emf->EnforceEMFBoundary();
    real dt=1.0;
    #ifdef EVOLVE_VECTOR_POTENTIAL
      hydro->emf->EvolveVectorPotential(dt, this->dA);
    #else
      hydro->emf->EvolveMagField(t, dt, this->dB);
    #endif
  }
  idfx::popRegion();
}

template<typename Phys>
template <int dir>
void RKLegendre<Phys>::CalcParabolicRHS(real t) {
  idfx::pushRegion("RKLegendre::CalcParabolicRHS");

  IdefixArray4D<real> Flux = hydro->FluxRiemann;
  IdefixArray3D<real> A    = data->A[dir];
  IdefixArray3D<real> dV   = data->dV;
  IdefixArray1D<real> x1m  = data->xl[IDIR];
  IdefixArray1D<real> x1   = data->x[IDIR];
  IdefixArray1D<real> sm   = data->sinx2m;
  IdefixArray1D<real> rt   = data->rt;
  IdefixArray1D<real> dmu  = data->dmu;
  IdefixArray1D<real> s    = data->sinx2;
  IdefixArray1D<real> dx   = data->dx[dir];
  IdefixArray1D<real> dx2  = data->dx[JDIR];
  IdefixArray3D<real> invDt = hydro->InvDt;
  IdefixArray3D<real> dMax = hydro->dMax;
  IdefixArray4D<real> viscSrc;
  IdefixArray4D<real> dU = this->dU;
  IdefixArray1D<int> varList = this->varList;

  bool haveViscosity = hydro->viscosityStatus.isRKL;
  if(haveViscosity) viscSrc = hydro->viscosity->viscSrc;
  IdefixArray4D<real> bragViscSrc;
  bool haveBragViscosity = hydro->bragViscosityStatus.isRKL;
  if(haveBragViscosity) bragViscSrc = hydro->bragViscosity->bragViscSrc;

  int ioffset,joffset,koffset;
  ioffset=joffset=koffset=0;
  // Determine the offset along which we do the extrapolation
  if(dir==IDIR) ioffset=1;
  if(dir==JDIR) joffset=1;
  if(dir==KDIR) koffset=1;


  idefix_for("CalcTotalFlux",
             0, this->nvarRKL,
             data->beg[KDIR],data->end[KDIR]+koffset,
             data->beg[JDIR],data->end[JDIR]+joffset,
             data->beg[IDIR],data->end[IDIR]+ioffset,
    KOKKOS_LAMBDA (int n, int k, int j, int i) {
      real Ax = A(k,j,i);

#if GEOMETRY != CARTESIAN
      if(Ax<SMALL_NUMBER)
        Ax=SMALL_NUMBER;    // Essentially to avoid singularity around poles
#endif

      const int nv = varList(n);

      Flux(nv,k,j,i) = Flux(nv,k,j,i) * Ax;

      // Curvature terms
#if    (GEOMETRY == POLAR       && COMPONENTS >= 2) \
    || (GEOMETRY == CYLINDRICAL && COMPONENTS == 3)
      if(dir==IDIR && nv==iMPHI) {
        // Conserve angular momentum, hence flux is R*Vphi
        Flux(iMPHI,k,j,i) = Flux(iMPHI,k,j,i) * FABS(x1m(i));
      }
#endif // GEOMETRY==POLAR OR CYLINDRICAL

#if GEOMETRY == SPHERICAL && COMPONENTS == 3
      if(dir==IDIR && nv==iMPHI) {
        Flux(iMPHI,k,j,i) = Flux(iMPHI,k,j,i) * FABS(x1m(i));
      } else if(dir==JDIR && nv==iMPHI) {
        Flux(iMPHI,k,j,i) = Flux(iMPHI,k,j,i) * FABS(sm(j));
      }
#endif // GEOMETRY == SPHERICAL && COMPONENTS == 3
    }
  );


  idefix_for("CalcRightHandSide",
             0, this->nvarRKL,
             data->beg[KDIR],data->end[KDIR],
             data->beg[JDIR],data->end[JDIR],
             data->beg[IDIR],data->end[IDIR],
    KOKKOS_LAMBDA (int n, int k, int j, int i) {
      constexpr const int ioffset = (dir==IDIR) ? 1 : 0;
      constexpr const int joffset = (dir==JDIR) ? 1 : 0;
      constexpr const int koffset = (dir==KDIR) ? 1 : 0;
      real rhs;

      const int nv = varList(n);

      rhs = -  ( Flux(nv, k+koffset, j+joffset, i+ioffset)
                     - Flux(nv, k, j, i))/dV(k,j,i);

      // Viscosity source terms
      if( haveViscosity && (nv-VX1 < COMPONENTS) && (nv-VX1>=0)) {
        rhs += viscSrc(nv-VX1,k,j,i);
      }
      // Braginskii Viscosity source terms
      if (haveBragViscosity && (nv-VX1 < COMPONENTS) && (nv-VX1>=0)) {
        rhs += bragViscSrc(nv-VX1,k,j,i);
      }

#if GEOMETRY != CARTESIAN
  #ifdef iMPHI
      if((dir==IDIR) && (nv == iMPHI)) {
        rhs /= x1(i);
      }
    #if (GEOMETRY == SPHERICAL) && (COMPONENTS == 3)
      if((dir==JDIR) && (nv == iMPHI)) {
        rhs /= FABS(s(j));
      }
    #endif // GEOMETRY
      // Nothing for KDIR
  #endif  // iMPHI
#endif // GEOMETRY != CARTESIAN

      // store the field components
      dU(nv,k,j,i) += rhs;
    });

  // Compute hyperbolic timestep only if we're in the first stage of the RKL loop
  if(stage==1) {
    // Grid coarsening
    bool haveGridCoarsening = false;
    IdefixArray2D<int> coarseningLevel;

    if(data->haveGridCoarsening) {
      haveGridCoarsening = data->coarseningDirection[dir];
      coarseningLevel = data->coarseningLevel[dir];
    }

    idefix_for("CalcDt",
             data->beg[KDIR],data->end[KDIR],
             data->beg[JDIR],data->end[JDIR],
             data->beg[IDIR],data->end[IDIR],
             KOKKOS_LAMBDA (int k, int j, int i) {
                constexpr const int ioffset = (dir==IDIR) ? 1 : 0;
                constexpr const int joffset = (dir==JDIR) ? 1 : 0;
                constexpr const int koffset = (dir==KDIR) ? 1 : 0;
               // Compute dt from max signal speed
                const int ig = ioffset*i + joffset*j + koffset*k;
                real dl = dx(ig);

                if(haveGridCoarsening) {
                  int factor;
                  //factor = 2^(coarsening-1)
                  if(dir==IDIR) {
                    factor = 1 << (coarseningLevel(k,j) - 1);
                  }
                  if(dir==JDIR) {
                    factor = 1 << (coarseningLevel(k,i) - 1);
                  }
                  if(dir==KDIR) {
                    factor = 1 << (coarseningLevel(j,i) - 1);
                  }
                  dl = dl * factor;
                }
                #if GEOMETRY == POLAR
                  if (dir==JDIR)
                    dl = dl*x1(i);

                #elif GEOMETRY == SPHERICAL
                  if (dir==JDIR)
                    dl = dl*rt(i);
                  else
                    if (dir==KDIR)
                      dl = dl*rt(i)*dmu(j)/dx2(j);
                 #endif

                invDt(k,j,i) += 0.5 * std::fmax(dMax(k+koffset,j+joffset,i+ioffset),
                                                        dMax(k,j,i)) / (dl*dl);
      });
  }

  idfx::popRegion();
}

template<typename Phys>
void RKLegendre<Phys>::SetBoundaries(real t) {
  idfx::pushRegion("RKLegendre::SetBoundaries");
  if(data->haveGridCoarsening) {
    hydro->CoarsenFlow(hydro->Vc);
    if constexpr(Phys::mhd) {
      hydro->CoarsenMagField(hydro->Vs);
    }
  }

  // set internal boundary conditions
  // Disabled since this might affect fields that are NOT updated
  // by the MPI instance of RKLegendre
  //if(hydro->boundary->haveInternalBoundary)
  //   hydro->boundary->internalBoundaryFunc(*data, t);
  for(int dir=0 ; dir < DIMENSIONS ; dir++ ) {
      // MPI Exchange data when needed
      // We use the RKL instance MPI object to ensure that we only exchange the data
      // solved by RKL
    #ifdef WITH_MPI
    if(data->mygrid->nproc[dir]>1) {
      switch(dir) {
        case 0:
          this->mpi.ExchangeX1(hydro->Vc, hydro->Vs);
          break;
        case 1:
          this->mpi.ExchangeX2(hydro->Vc, hydro->Vs);
          break;
        case 2:
          this->mpi.ExchangeX3(hydro->Vc, hydro->Vs);
          break;
      }
    }
    #endif
    hydro->boundary->EnforceBoundaryDir(t, dir);
    if constexpr(Phys::mhd) {
      // Reconstruct the normal field component when using CT
      if(haveVs) {
        hydro->boundary->ReconstructNormalField(dir);
      }
    }
  } // Loop on dimension ends

  if constexpr(Phys::mhd) {
    // Remake the cell-centered field.
    if(haveVs) {
      hydro->boundary->ReconstructVcField(hydro->Vc);
    }
  }
  idfx::popRegion();
}


#endif // RKL_RKL_HPP_
