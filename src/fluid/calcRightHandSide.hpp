// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef FLUID_CALCRIGHTHANDSIDE_HPP_
#define FLUID_CALCRIGHTHANDSIDE_HPP_

#include "fluid.hpp"
#include "dataBlock.hpp"
#include "gravity.hpp"

template<typename Phys, int dir>
struct Fluid_CorrectFluxFunctor {
  // Correct the flux to take into account non-cartesian geometries and Fargo
  //*****************************************************************
  // Functor constructor
  //*****************************************************************
  explicit Fluid_CorrectFluxFunctor (Fluid<Phys> *hydro, real dt) {
    Uc   = hydro->Uc;
    Vc   = hydro->Vc;
    Flux = hydro->FluxRiemann;
    A    = hydro->data->A[dir];
    dV   = hydro->data->dV;
    x1m  = hydro->data->xl[IDIR];
    x1   = hydro->data->x[IDIR];
    sinx2m   = hydro->data->sinx2m;
    sinx2 = hydro->data->sinx2;

    this->dt = dt;
    // Fargo
    haveFargo  = hydro->data->haveFargo;
    if(haveFargo) {
      fargoVelocity = hydro->data->fargo->meanVelocity;
      fargoType = hydro->data->fargo->type;
    }

    // Rotation
    haveRotation = hydro->haveRotation;
    // disable rotation in cartesian geometry, as Coriolis is then treated as a source term
    #if GEOMETRY == CARTESIAN
      haveRotation = false;
    #endif
    Omega = hydro->OmegaZ;

    // Shearing box shear rate
    sbS = hydro->sbS;
  }

  //*****************************************************************
  // Functor Variables
  //*****************************************************************
  IdefixArray4D<real> Uc;
  IdefixArray4D<real> Vc;
  IdefixArray4D<real> Flux;
  IdefixArray3D<real> A;
  IdefixArray3D<real> dV;
  IdefixArray1D<real> x1m;
  IdefixArray1D<real> x1;
  IdefixArray1D<real> sinx2m;
  IdefixArray1D<real> sinx2;

  // Fargo
  IdefixArray2D<real> fargoVelocity;
  Fargo::FargoType fargoType;
  bool haveFargo;


  //Rotation
  bool haveRotation;
  real Omega;

  // shearingBox
  real sbS;

  // timestep
  real dt;

  //*****************************************************************
  // Functor Operator
  //*****************************************************************
  KOKKOS_INLINE_FUNCTION void operator() (const int k, const int j,  const int i) const {
      // Add Fargo velocity to the fluxes
      if(haveFargo || haveRotation) {
        // Set mean advection direction
        #if (GEOMETRY == CARTESIAN || GEOMETRY == POLAR) && DIMENSIONS >=2
          const int meanDir = JDIR;
        #elif GEOMETRY == SPHERICAL && DIMENSIONS == 3
          const int meanDir = KDIR;
        #else
          const int meanDir = 0;
        #endif

        // set mean advection velocity
        real meanV = ZERO_F;
        if constexpr(dir == IDIR) {
          #if (GEOMETRY == CARTESIAN || GEOMETRY == POLAR) && DIMENSIONS >=2
            if(haveFargo) {
              if(fargoType==Fargo::userdef) {
                meanV = HALF_F*(fargoVelocity(k,i-1)+fargoVelocity(k,i));
              } else if(fargoType==Fargo::shearingbox) {
                meanV = sbS * x1m(i);
              }
            }
            #if GEOMETRY != CARTESIAN
            if(haveRotation) {
              meanV += x1m(i)*Omega;
            }
            #endif
          #elif GEOMETRY == SPHERICAL && DIMENSIONS == 3
            if(haveFargo) {
              meanV = HALF_F*(fargoVelocity(j,i-1)+fargoVelocity(j,i));
            }
            if(haveRotation) {
              meanV += x1m(i)*sinx2(j)*Omega;
            }
          #endif// GEOMETRY
        }
        #if GEOMETRY == SPHERICAL && DIMENSIONS == 3
          if constexpr (dir == JDIR) {
            if(haveFargo) {
              meanV = HALF_F*(fargoVelocity(j-1,i)+fargoVelocity(j,i));
            }
            if(haveRotation) {
              meanV += x1(i)*sinx2m(j)*Omega;
            }
          }
        #elif (GEOMETRY == CARTESIAN || GEOMETRY == POLAR) && DIMENSIONS >=2
          if constexpr (dir == KDIR) {
            if(haveFargo) {
              if(fargoType==Fargo::userdef) {
                meanV = HALF_F*(fargoVelocity(k-1,i)+fargoVelocity(k,i));
              } else if (fargoType==Fargo::shearingbox) {
                meanV = sbS*x1(i);
              }
            }
          }
        #endif // GEOMETRY

        // Should do nothing along meanV direction, automatically satisfied
        // since in that case meanV=0
        if constexpr(Phys::pressure) {
          // Mignone (2012): second and third term of rhs of (25)
          Flux(ENG,k,j,i) += meanV * (HALF_F*meanV*Flux(RHO,k,j,i) + Flux(MX1+meanDir,k,j,i));
        }
        // Mignone+2012: second term of rhs of (24)
        Flux(MX1+meanDir,k,j,i) += meanV * Flux(RHO,k,j,i);
      } // Fargo & Rotation corrections

      real Ax = A(k,j,i);

      for(int nv = 0 ; nv < Phys::nvar ; nv++) {
        Flux(nv,k,j,i) = Flux(nv,k,j,i) * Ax;
      }

      // Curvature terms
#if    (GEOMETRY == POLAR       && COMPONENTS >= 2) \
    || (GEOMETRY == CYLINDRICAL && COMPONENTS == 3)
      if constexpr (dir==IDIR) {
        // Conserve angular momentum, hence flux is R*Bphi
        Flux(iMPHI,k,j,i) = Flux(iMPHI,k,j,i) * FABS(x1m(i));
        if constexpr(Phys::mhd) {
          if(Ax<SMALL_NUMBER) Ax=SMALL_NUMBER;    //avoid singularity around poles
          // No area for this one
          Flux(iBPHI,k,j,i) = Flux(iBPHI,k,j,i) / Ax;
        }
      }
#endif // GEOMETRY==POLAR OR CYLINDRICAL

#if GEOMETRY == SPHERICAL
      if constexpr(dir==IDIR) {
  #if COMPONENTS == 3
        Flux(iMPHI,k,j,i) = Flux(iMPHI,k,j,i) * FABS(x1m(i));
  #endif // COMPONENTS == 3
        if constexpr(Phys::mhd) {
          if(Ax<SMALL_NUMBER) Ax=SMALL_NUMBER;    // avoid singularity around poles
          EXPAND(                                            ,
              Flux(iBTH,k,j,i)  = Flux(iBTH,k,j,i) * x1m(i) / Ax;  ,
              Flux(iBPHI,k,j,i) = Flux(iBPHI,k,j,i) * x1m(i) / Ax; )
        }
      } else if constexpr (dir==JDIR) {
  #if COMPONENTS == 3
        Flux(iMPHI,k,j,i) = Flux(iMPHI,k,j,i) * FABS(sinx2m(j));
        if constexpr(Phys::mhd) {
          if(Ax<SMALL_NUMBER) Ax=SMALL_NUMBER;    // avoid singularity around poles
          Flux(iBPHI,k,j,i) = Flux(iBPHI,k,j,i)  / Ax;
        }
  #endif // COMPONENTS = 3
      }
#endif // GEOMETRY == SPHERICAL
    }
};



template<typename Phys, int dir>
struct Fluid_CalcRHSFunctor {
  //*****************************************************************
  // Functor constructor
  //*****************************************************************
  explicit Fluid_CalcRHSFunctor (Fluid<Phys> *hydro, real dt) {
    Uc   = hydro->Uc;
    Vc   = hydro->Vc;
    Flux = hydro->FluxRiemann;
    A    = hydro->data->A[dir];
    dV   = hydro->data->dV;
    x1m  = hydro->data->xl[IDIR];
    x1   = hydro->data->x[IDIR];

    rt   = hydro->data->rt;
    dmu  = hydro->data->dmu;
    sinx2m   = hydro->data->sinx2m;
    sinx2 = hydro->data->sinx2;
    dx   = hydro->data->dx[dir];
    dx2  = hydro->data->dx[JDIR];
    invDt = hydro->InvDt;
    cMax = hydro->cMax;
    dMax = hydro->dMax;
    this->dt = dt;

    // Grid coarsening
    if(hydro->data->haveGridCoarsening) {
      haveGridCoarsening = hydro->data->coarseningDirection[dir];
      coarseningLevel = hydro->data->coarseningLevel[dir];
    }

    if(hydro->data->haveGravity) {
      // Gravitational potential
      phiP = hydro->data->gravity->phiP;
      needPotential = hydro->data->gravity->havePotential;

      // BodyForce
      bodyForce = hydro->data->gravity->bodyForceVector;
      needBodyForce = hydro->data->gravity->haveBodyForce;
    }

    // parabolic terms
    haveParabolicTerms = hydro->haveExplicitParabolicTerms;

    // Viscosity  (source term only when non-cartesian geometry)
    haveViscosity = hydro->viscosityStatus.isExplicit;
    if(haveViscosity) viscSrc = hydro->viscosity->viscSrc;

    // BragViscosity  (source term only when non-cartesian geometry)
    haveBragViscosity = hydro->bragViscosityStatus.isExplicit;
    if(haveBragViscosity) bragViscSrc = hydro->bragViscosity->bragViscSrc;

    // Fargo
    haveFargo  = hydro->data->haveFargo;
    if(haveFargo) {
      fargoVelocity = hydro->data->fargo->meanVelocity;
      fargoType = hydro->data->fargo->type;
    }

    // Rotation
    haveRotation = hydro->haveRotation;
    // disable rotation in cartesian geometry, as Coriolis is then treated as a source term
    #if GEOMETRY == CARTESIAN
      haveRotation = false;
    #endif
    Omega = hydro->OmegaZ;

    // Shearing box shear rate
    sbS = hydro->sbS;
  }
  //*****************************************************************
  // Functor Variables
  //*****************************************************************
  IdefixArray4D<real> Uc;
  IdefixArray4D<real> Vc;
  IdefixArray4D<real> Flux;
  IdefixArray3D<real> A;
  IdefixArray3D<real> dV;
  IdefixArray1D<real> x1m;
  IdefixArray1D<real> x1;

  IdefixArray1D<real> rt;
  IdefixArray1D<real> dmu;
  IdefixArray1D<real> sinx2m;
  IdefixArray1D<real> sinx2;
  IdefixArray1D<real> dx;
  IdefixArray1D<real> dx2;
  IdefixArray3D<real> invDt;
  IdefixArray3D<real> cMax;
  IdefixArray3D<real> dMax;
  IdefixArray4D<real> viscSrc;

  // Grid coarsening
  bool haveGridCoarsening{false};
  IdefixArray2D<int> coarseningLevel;

  // Gravitational potential
  IdefixArray3D<real> phiP;
  bool needPotential{false};

  // BodyForce
  IdefixArray4D<real> bodyForce;
  bool needBodyForce{false};

  // parabolic terms
  bool haveParabolicTerms{false};

  // Viscosity  (source term only when non-cartesian geometry)
  bool haveViscosity{false};

  // BragViscosity  (source term only when non-cartesian geometry)
  IdefixArray4D<real> bragViscSrc;
  bool haveBragViscosity{false};

  // Fargo
  IdefixArray2D<real> fargoVelocity;
  Fargo::FargoType fargoType;
  bool haveFargo;


  //Rotation
  bool haveRotation;
  real Omega;

  // shearingBox
  real sbS;

  // timestep
  real dt;

  //*****************************************************************
  // Functor Operator
  //*****************************************************************
  KOKKOS_INLINE_FUNCTION void operator() (const int k, const int j,  const int i) const {
    const int ioffset = (dir==IDIR) ? 1 : 0;
    const int joffset = (dir==JDIR) ? 1 : 0;
    const int koffset = (dir==KDIR) ? 1 : 0;

    real dtdV=dt / dV(k,j,i);
    real rhs[Phys::nvar];

    #pragma unroll
    for(int nv = 0 ; nv < Phys::nvar ; nv++) {
      rhs[nv] = -  dtdV*(Flux(nv, k+koffset, j+joffset, i+ioffset) - Flux(nv, k, j, i));
    }

    #if GEOMETRY != CARTESIAN
      if constexpr(dir==IDIR) {
        #ifdef iMPHI
          rhs[iMPHI] = rhs[iMPHI] / x1(i);
        #endif
        if constexpr(Phys::mhd) {
          #if (GEOMETRY == POLAR || GEOMETRY == CYLINDRICAL) &&  (defined iBPHI)
            rhs[iBPHI] = - dt / dx(i) * (Flux(iBPHI, k, j, i+1) - Flux(iBPHI, k, j, i) );

          #elif (GEOMETRY == SPHERICAL)
            real q = dt / (x1(i)*dx(i));
            EXPAND(                                                                       ,
                  rhs[iBTH]  = -q * ((Flux(iBTH, k, j, i+1)  - Flux(iBTH, k, j, i) ));  ,
                  rhs[iBPHI] = -q * ((Flux(iBPHI, k, j, i+1) - Flux(iBPHI, k, j, i) )); )
          #endif
        } // MHD
      } else if constexpr(dir==JDIR) {
        #if (GEOMETRY == SPHERICAL) && (COMPONENTS == 3)
          rhs[iMPHI] /= FABS(sinx2(j));
          if constexpr(Phys::mhd) {
            rhs[iBPHI] = -dt / (rt(i)*dx(j)) * (Flux(iBPHI, k, j+1, i) - Flux(iBPHI, k, j, i));
          } // MHD
        #endif // GEOMETRY
      }
      // Nothing for KDIR
    #endif // GEOMETRY != CARTESIAN


    // Fargo terms to enfore conservation (actually substract back what was added in
    // Totalflux loop)
    if(haveFargo || haveRotation) {
      // fetch fargo velocity when required
      real meanV = ZERO_F;
      #if (GEOMETRY == POLAR || GEOMETRY == CARTESIAN) && DIMENSIONS >=2
        if((dir==IDIR || dir == KDIR) && haveFargo) {
          if(fargoType==Fargo::userdef) {
            meanV = fargoVelocity(k,i);
          } else if(fargoType==Fargo::shearingbox) {
            meanV = sbS * x1(i);
          }
        }
        #if GEOMETRY != CARTESIAN
          if((dir==IDIR) && haveRotation) {
            meanV += Omega*x1(i);
          }
        #endif
        const int meanDir = JDIR;
      #elif GEOMETRY == SPHERICAL && DIMENSIONS ==3
        if((dir==IDIR || dir == JDIR) && haveFargo) meanV = fargoVelocity(j,i);
        if((dir==IDIR || dir == JDIR) && haveRotation) {
          meanV += Omega*x1(i)*sinx2(j);
        }
        const int meanDir = KDIR;
      #else
        const int meanDir = 0;
      #endif
      // NB: MX1+meanDir = iMPHI
      // Mignone+2012, rhs of eq. 24, last term
      rhs[MX1+meanDir] -= meanV*rhs[RHO];
      if constexpr(Phys::pressure) {
        // Mignone+2012, rhs of eq.25, 4th and 5th term. NB: rhs[MX1+meanDir] is alreay the
        // divergence of the *total* Momentum Flux F_my+w F_rho
        rhs[ENG] -= meanV * ( HALF_F*meanV*rhs[RHO] + rhs[MX1+meanDir]);
      }
    }

    #if GEOMETRY != CARTESIAN
      // Viscosity source terms
      if(haveViscosity) {
        #pragma unroll
        for(int nv = 0 ; nv < COMPONENTS ; nv++) {
          rhs[nv + VX1] += dt*viscSrc(nv,k,j,i);
        }
      } else if(haveBragViscosity) {
        #pragma unroll
        for(int nv = 0 ; nv < COMPONENTS ; nv++) {
          rhs[nv + VX1] += dt*bragViscSrc(nv,k,j,i);
        }
      }
    #endif // GEOMETRY != CARTESIAN

    // elmentary length for gradient computations
    const int ig = ioffset*i + joffset*j + koffset*k;
    real dl = dx(ig);
    #if GEOMETRY == POLAR
      if constexpr (dir==JDIR)
        dl = dl*x1(i);

    #elif GEOMETRY == SPHERICAL
      if constexpr(dir==JDIR)
        dl = dl*rt(i);
      else if constexpr(dir==KDIR)
          dl = dl*rt(i)*dmu(j)/dx2(j);
    #endif

    // Potential terms
    if(needPotential) {
      real dphi;
      if constexpr (dir==IDIR) {
        // Gravitational force in direction i
        dphi = - 1.0/12.0 * (
                      - phiP(k,j,i+2) + 8.0 * phiP(k,j,i+1)
                      - 8.0*phiP(k,j,i-1) + phiP(k,j,i-2));
      }
      if constexpr (dir==JDIR) {
        // Gravitational force in direction j
        dphi = - 1.0/12.0 * (
                      - phiP(k,j+2,i) + 8.0 * phiP(k,j+1,i)
                      - 8.0*phiP(k,j-1,i) + phiP(k,j-2,i));
      }
      if constexpr (dir==KDIR) {
        // Gravitational force in direction k
        dphi = - 1.0/12.0 * (
                      - phiP(k+2,j,i) + 8.0 * phiP(k+1,j,i)
                      - 8.0*phiP(k-1,j,i) + phiP(k-2,j,i));
      }
      rhs[MX1+dir] += dt * Vc(RHO,k,j,i) * dphi /dl;

      if constexpr(Phys::pressure) {
        // Add gravitational force work as a source term
        // This is equivalent to rho * v . nabla(phi)
        // (note that Flux has already been multiplied by A)
        rhs[ENG] += HALF_F * dtdV  *
                  (Flux(RHO,k,j,i) + Flux(RHO, k+koffset, j+joffset, i+ioffset)) * dphi;
      }
    }

    // Body force
    if(needBodyForce) {
      rhs[MX1+dir] += dt * Vc(RHO,k,j,i) * bodyForce(dir,k,j,i);
      if constexpr(Phys::pressure) {
        //  rho * v . f, where rhov is taken as a  volume average of Flux(RHO)
        rhs[ENG] += HALF_F * dtdV * dl *
                      (Flux(RHO,k,j,i) + Flux(RHO, k+koffset, j+joffset, i+ioffset)) *
                        bodyForce(dir,k,j,i);
      } // Pressure

      // Particular cases if we do not sweep all of the components
      #if DIMENSIONS == 1 && COMPONENTS > 1
        EXPAND(                                                           ,
                  rhs[MX2] += dt * Vc(RHO,k,j,i) * bodyForce(JDIR,k,j,i);   ,
                  rhs[MX3] += dt * Vc(RHO,k,j,i) * bodyForce(KDIR,k,j,i);    )
        if constexpr(Phys::pressure) {
          rhs[ENG] += dt * (EXPAND( ZERO_F                                              ,
                                    + Vc(RHO,k,j,i) * Vc(VX2,k,j,i) * bodyForce(JDIR,k,j,i)   ,
                                    + Vc(RHO,k,j,i) * Vc(VX3,k,j,i) * bodyForce(KDIR,k,j,i) ));
        }
      #endif
      #if DIMENSIONS == 2 && COMPONENTS == 3
        // Only add this term once!
        if constexpr (dir==JDIR) {
          rhs[MX3] += dt * Vc(RHO,k,j,i) * bodyForce(KDIR,k,j,i);
          if constexpr(Phys::pressure) {
            rhs[ENG] += dt * Vc(RHO,k,j,i) * Vc(VX3,k,j,i) * bodyForce(KDIR,k,j,i);
          }
        }
      #endif
    }


    // Timestep computation
    // Change elementary grid spacing according to local coarsening level.
    if(haveGridCoarsening) {
      int factor;
      //factor = 2^(coarsening-1)
      if constexpr (dir==IDIR) {
        factor = 1 << (coarseningLevel(k,j) - 1);
      }
      if constexpr (dir==JDIR) {
        factor = 1 << (coarseningLevel(k,i) - 1);
      }
      if constexpr (dir==KDIR) {
        factor = 1 << (coarseningLevel(j,i) - 1);
      }
      dl = dl * factor;
    }

    // Compute dt from max signal speed
    invDt(k,j,i) = invDt(k,j,i) + HALF_F*(cMax(k+koffset,j+joffset,i+ioffset)
                  + cMax(k,j,i)) / (dl);

    if(haveParabolicTerms) {
      invDt(k,j,i) = invDt(k,j,i) + TWO_F* FMAX(dMax(k+koffset,j+joffset,i+ioffset),
                                                  dMax(k,j,i)) / (dl*dl);
    }


    // Evolve the field components
    #pragma unroll
    for(int nv = 0 ; nv < Phys::nvar ; nv++) {
      // Do not evolve the field components if they are computed by CT (i.e. if they are in Vs)
      D_EXPAND( if(nv == BX1) { continue; }  ,
                if(nv == BX2) { continue; }  ,
                if(nv == BX3) { continue; }  )


      Uc(nv,k,j,i) = Uc(nv,k,j,i) + rhs[nv];
    }
  }
};



// Compute the right handside in direction dir from conservative equation, with timestep dt
template<typename Phys>
template<int dir>
void Fluid<Phys>::CalcRightHandSide(real t, real dt) {
  idfx::pushRegion("Fluid::CalcRightHandSide");

  // Update fargo velocity when needed
  if(data->haveFargo && data->fargo->type == Fargo::userdef) {
    data->fargo->GetFargoVelocity(t);
  }

  auto fluxCorrection = Fluid_CorrectFluxFunctor<Phys,dir>(this,dt);

  /////////////////////////////////////////////////////////////////////////////
  // Flux correction (for fargo/non-cartesian geometry)
  /////////////////////////////////////////////////////////////////////////////
  const int ioffset = (dir==IDIR) ? 1 : 0;
  const int joffset = (dir==JDIR) ? 1 : 0;
  const int koffset = (dir==KDIR) ? 1 : 0;
  idefix_for("Correct Flux",
             data->beg[KDIR],data->end[KDIR]+koffset,
             data->beg[JDIR],data->end[JDIR]+joffset,
             data->beg[IDIR],data->end[IDIR]+ioffset,
              fluxCorrection);


  // If user has requested specific flux functions for the boundaries, here they come
  if(boundary->haveFluxBoundary) boundary->EnforceFluxBoundaries(dir,t);

  auto calcRHS = Fluid_CalcRHSFunctor<Phys,dir>(this,dt);
  /////////////////////////////////////////////////////////////////////////////
  // Final conserved quantity budget from fluxes divergence
  /////////////////////////////////////////////////////////////////////////////
  idefix_for("CalcRightHandSide",
             data->beg[KDIR],data->end[KDIR],
             data->beg[JDIR],data->end[JDIR],
             data->beg[IDIR],data->end[IDIR],
              calcRHS);


  idfx::popRegion();
}
#endif // FLUID_CALCRIGHTHANDSIDE_HPP_
