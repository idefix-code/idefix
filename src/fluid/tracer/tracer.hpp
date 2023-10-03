// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef FLUID_TRACER_TRACER_HPP_
#define FLUID_TRACER_TRACER_HPP_

#include <string>
#include "idefix.hpp"
#include "slopeLimiter.hpp"

// Forward class hydro declaration
template <typename Phys> class Fluid;
class DataBlock;


class Tracer {
 public:
  template <typename Phys> Tracer(Fluid<Phys> *, int n);
  void ConvertConsToPrim();
  void ConvertPrimToCons();
  template <int, typename> void CalcFlux(IdefixArray4D<real> &);
  template <int, typename> void CalcRightHandSide(IdefixArray4D<real> &, real, real);

 private:
  IdefixArray4D<real> Vc;  // Vector of primitive variables for the passive tracer
  IdefixArray4D<real> Uc;  // Vector of conservative variables for the passive tracer

  std::string prefix;

  DataBlock *data;

  int nTracer;
  int nVar;
};


#include "fluid.hpp"
#include "dataBlock.hpp"

template <typename Phys>
Tracer::Tracer(Fluid<Phys> *fluid, int n) {
  idfx::pushRegion("Tracer::Tracer");
  nTracer = n;
  Vc = fluid->Vc;
  Uc = fluid->Uc;
  data = fluid->data;
  nVar = Phys::nvar;
  idfx::popRegion();
}


// Compute the upwinded flux
template <int dir, typename Phys>
void Tracer::CalcFlux(IdefixArray4D<real> &Flux) {
  idfx::pushRegion("Tracer::CalcFlux");

  IdefixArray4D<real> Vc = this->Vc;
  IdefixArray4D<real> Uc = this->Uc;
  IdefixArray3D<real> A    = data->A[dir];

  constexpr int ioffset = (dir==IDIR ? 1 : 0);
  constexpr int joffset = (dir==JDIR ? 1 : 0);
  constexpr int koffset = (dir==KDIR ? 1 : 0);

  idefix_for("ComputeTracerFlux",
             Phys::nvar, Phys::nvar+nTracer,   // Loop on the index where tracers are lying
             data->beg[KDIR],data->end[KDIR]+koffset,
             data->beg[JDIR],data->end[JDIR]+joffset,
             data->beg[IDIR],data->end[IDIR]+ioffset,
    KOKKOS_LAMBDA (int nv, int k, int j, int i) {
      real vface;
      if(Flux(RHO,k,j,i) > 0) {
        // Interpolate from the left
        real dvm = Vc(nv,k-koffset,j-joffset,i-ioffset)
                  -Vc(nv,k-2*koffset,j-2*joffset,i-2*ioffset);
        real dvp = Vc(nv,k,j,i)-Vc(nv,k-koffset,j-joffset,i-ioffset);

        real dv = SlopeLimiter<>::PLMLim(dvp,dvm);

        vface = Vc(nv,k-koffset,j-joffset,i-ioffset) + HALF_F*dv;
      } else {
        // interpolation from the right
        real dvm = Vc(nv,k,j,i)
                  -Vc(nv,k-koffset,j-joffset,i-ioffset);
        real dvp = Vc(nv,k+koffset,j+joffset,i+ioffset) - Vc(nv,k,j,i);

        real dv = SlopeLimiter<>::PLMLim(dvp,dvm);

        vface =  Vc(nv,k,j,i) - HALF_F*dv;
      }

      // Compute flux*Area (! different from the fluid flux)
      Flux(nv,k,j,i) = Flux(RHO,k,j,i) * vface * A(k,j,i);
  });
  idfx::popRegion();
}

template <int dir, typename Phys>
void Tracer::CalcRightHandSide(IdefixArray4D<real> &Flux, real t, real dt) {
  idfx::pushRegion("Tracer::ComputeRHS");

  IdefixArray4D<real> Uc = this->Uc;
  IdefixArray3D<real> dV  = data->dV;

  constexpr int ioffset = (dir==IDIR ? 1 : 0);
  constexpr int joffset = (dir==JDIR ? 1 : 0);
  constexpr int koffset = (dir==KDIR ? 1 : 0);

  idefix_for("ComputeTracerRHS",
             Phys::nvar, Phys::nvar+nTracer,   // Loop on the index where tracers are lying
             data->beg[KDIR],data->end[KDIR],
             data->beg[JDIR],data->end[JDIR],
             data->beg[IDIR],data->end[IDIR],
    KOKKOS_LAMBDA (int nv, int k, int j, int i) {
      Uc(nv,k,j,i) += -dt / dV(k,j,i) * (Flux(nv,k+koffset,j+joffset,i+ioffset) - Flux(nv,k,j,i));
  });

  idfx::popRegion();
}

#endif // FLUID_TRACER_TRACER_HPP_
