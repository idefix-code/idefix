// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2022 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef HYDRO_CALCRIGHTHANDSIDE_HPP_
#define HYDRO_CALCRIGHTHANDSIDE_HPP_

#include "hydro.hpp"
#include "dataBlock.hpp"


// Compute the right handside in direction dir from conservative equation, with timestep dt
template<int dir>
void Hydro::CalcRightHandSide(real t, real dt) {
  idfx::pushRegion("Hydro::CalcRightHandSide");

  IdefixArray4D<real> Uc   = this->Uc;
  IdefixArray4D<real> Vc   = this->Vc;
  IdefixArray4D<real> Flux = this->FluxRiemann;
  IdefixArray3D<real> A    = data->A[dir];
  IdefixArray3D<real> dV   = data->dV;
  IdefixArray1D<real> x1m  = data->xl[IDIR];
  IdefixArray1D<real> x1   = data->x[IDIR];

  IdefixArray1D<real> rt   = data->rt;
  IdefixArray1D<real> dmu  = data->dmu;
  IdefixArray1D<real> sinx2m   = data->sinx2m;
  IdefixArray1D<real> sinx2 = data->sinx2;
  IdefixArray1D<real> dx   = data->dx[dir];
  IdefixArray1D<real> dx2  = data->dx[JDIR];
  IdefixArray3D<real> invDt = this->InvDt;
  IdefixArray3D<real> cMax = this->cMax;
  IdefixArray3D<real> dMax = this->dMax;
  IdefixArray4D<real> viscSrc = this->viscosity.viscSrc;
  IdefixArray2D<real> fargoVelocity = data->fargo.meanVelocity;

  // Grid coarsening
  bool haveGridCoarsening = false;
  IdefixArray2D<int> coarseningLevel;

  if(data->haveGridCoarsening) {
    haveGridCoarsening = data->coarseningDirection[dir];
    coarseningLevel = data->coarseningLevel[dir];
  }

  // Gravitational potential
  IdefixArray3D<real> phiP = data->gravity.phiP;
  bool needPotential = data->gravity.havePotential;

  // BodyForce
  IdefixArray4D<real> bodyForce = data->gravity.bodyForceVector;
  bool needBodyForce = data->gravity.haveBodyForce;

  // parabolic terms
  bool haveParabolicTerms = this->haveExplicitParabolicTerms;

  // Viscosity  (source term only when non-cartesian geometry)
  [[maybe_unused]] bool haveViscosity = this->viscosityStatus.isExplicit;

  // Fargo
  bool haveFargo  = data->haveFargo;
  Fargo::FargoType fargoType = data->fargo.type;

  //Rotation
  bool haveRotation = this->haveRotation;
  [[maybe_unused]] real Omega = this->OmegaZ;
  // disable rotation in cartesian geometry, as Coriolis is then treated as a source term
  #if GEOMETRY == CARTESIAN
    haveRotation = false;
  #endif

  // shearingBox
  [[maybe_unused]] bool haveShearingBox = this->haveShearingBox;
  [[maybe_unused]] real sbS = this->sbS;

  if(haveFargo && fargoType == Fargo::userdef) {
    data->fargo.GetFargoVelocity(t);
  }

  constexpr const int ioffset = (dir==IDIR) ? 1 : 0;
  constexpr const int joffset = (dir==JDIR) ? 1 : 0;
  constexpr const int koffset = (dir==KDIR) ? 1 : 0;

  idefix_for("CalcTotalFlux",
             data->beg[KDIR],data->end[KDIR]+koffset,
             data->beg[JDIR],data->end[JDIR]+joffset,
             data->beg[IDIR],data->end[IDIR]+ioffset,
    KOKKOS_LAMBDA (int k, int j, int i) {
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
        if(dir == IDIR) {
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
          if((dir == JDIR)) {
            if(haveFargo) {
              meanV = HALF_F*(fargoVelocity(j-1,i)+fargoVelocity(j,i));
            }
            if(haveRotation) {
              meanV += x1(i)*sinx2m(j)*Omega;
            }
          }
        #elif (GEOMETRY == CARTESIAN || GEOMETRY == POLAR) && DIMENSIONS >=2
          if((dir == KDIR) && haveFargo) {
            if(fargoType==Fargo::userdef) {
              meanV = HALF_F*(fargoVelocity(k-1,i)+fargoVelocity(k,i));
            } else if (fargoType==Fargo::shearingbox) {
              meanV = sbS*x1(i);
            }
          }
        #endif // GEOMETRY

        // Should do nothing along meanV direction, automatically satisfied
        // since in that case meanV=0

        #if HAVE_ENERGY
          Flux(ENG,k,j,i) += meanV * (HALF_F*meanV*Flux(RHO,k,j,i) + Flux(MX1+meanDir,k,j,i));
        #endif
        Flux(MX1+meanDir,k,j,i) += meanV * Flux(RHO,k,j,i);
      } // have Fargo

      real Ax = A(k,j,i);

#if GEOMETRY != CARTESIAN
      if(Ax<SMALL_NUMBER)
        Ax=SMALL_NUMBER;    // Essentially to avoid singularity around poles
#endif

      for(int nv = 0 ; nv < NVAR ; nv++) {
        Flux(nv,k,j,i) = Flux(nv,k,j,i) * Ax;
      }


      // Curvature terms
#if    (GEOMETRY == POLAR       && COMPONENTS >= 2) \
    || (GEOMETRY == CYLINDRICAL && COMPONENTS == 3)
      if(dir==IDIR) {
        // Conserve angular momentum, hence flux is R*Bphi
        Flux(iMPHI,k,j,i) = Flux(iMPHI,k,j,i) * FABS(x1m(i));
  #if MHD == YES
        // No area for this one
        Flux(iBPHI,k,j,i) = Flux(iBPHI,k,j,i) / Ax;
  #endif // MHD
      }
#endif // GEOMETRY==POLAR OR CYLINDRICAL

#if GEOMETRY == SPHERICAL
      if(dir==IDIR) {
  #if COMPONENTS == 3
        Flux(iMPHI,k,j,i) = Flux(iMPHI,k,j,i) * FABS(x1m(i));
  #endif // COMPONENTS == 3
  #if MHD == YES
        EXPAND(                                            ,
            Flux(iBTH,k,j,i)  = Flux(iBTH,k,j,i) * x1m(i) / Ax;  ,
            Flux(iBPHI,k,j,i) = Flux(iBPHI,k,j,i) * x1m(i) / Ax; )
  #endif // MHD
      } else if(dir==JDIR) {
  #if COMPONENTS == 3
        Flux(iMPHI,k,j,i) = Flux(iMPHI,k,j,i) * FABS(sinx2m(j));
    #if MHD == YES
        Flux(iBPHI,k,j,i) = Flux(iBPHI,k,j,i)  / Ax;
    #endif // MHD
  #endif // COMPONENTS = 3
      }
#endif // GEOMETRY == SPHERICAL
    }
  );

  // If user has requested specific flux functions for the boundaries, here they come
  if(boundary.haveFluxBoundary) boundary.EnforceFluxBoundaries(dir);

  idefix_for("CalcRightHandSide",
             data->beg[KDIR],data->end[KDIR],
             data->beg[JDIR],data->end[JDIR],
             data->beg[IDIR],data->end[IDIR],
    KOKKOS_LAMBDA (int k, int j, int i) {
      real dtdV=dt / dV(k,j,i);
      real rhs[NVAR];

#pragma unroll
      for(int nv = 0 ; nv < NVAR ; nv++) {
        rhs[nv] = -  dtdV*(Flux(nv, k+koffset, j+joffset, i+ioffset) - Flux(nv, k, j, i));
      }

#if GEOMETRY != CARTESIAN
      if(dir==IDIR) {
  #ifdef iMPHI
        rhs[iMPHI] = rhs[iMPHI] / x1(i);
  #endif

  #if (GEOMETRY == POLAR || GEOMETRY == CYLINDRICAL) &&  (defined iBPHI) && (MHD == YES)
        rhs[iBPHI] = - dt / dx(i) * (Flux(iBPHI, k, j, i+1) - Flux(iBPHI, k, j, i) );

  #elif (GEOMETRY == SPHERICAL) && (MHD == YES)
        real q = dt / (x1(i)*dx(i));
        EXPAND(                                                                       ,
                rhs[iBTH]  = -q * ((Flux(iBTH, k, j, i+1)  - Flux(iBTH, k, j, i) ));  ,
                rhs[iBPHI] = -q * ((Flux(iBPHI, k, j, i+1) - Flux(iBPHI, k, j, i) )); )
  #endif
      } else if(dir==JDIR) {
  #if (GEOMETRY == SPHERICAL) && (COMPONENTS == 3)
        rhs[iMPHI] /= FABS(sinx2(j));
    #if MHD == YES
        rhs[iBPHI] = -dt / (rt(i)*dx(j)) * (Flux(iBPHI, k, j+1, i) - Flux(iBPHI, k, j, i));
    #endif // MHD
  #endif // GEOMETRY
      }
      // Nothing for KDIR

      // Viscosity source terms
      if(haveViscosity) {
#pragma unroll
        for(int nv = 0 ; nv < COMPONENTS ; nv++) {
          rhs[nv + VX1] += dt*viscSrc(nv,k,j,i);
        }
      }

#endif // GEOMETRY != CARTESIAN

      // Compute dt from max signal speed
      const int ig = ioffset*i + joffset*j + koffset*k;
      real dl = dx(ig);

      // Change elementary grid spacing according to local coarsening level.

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
      if(dir==JDIR)
        dl = dl*x1(i);

#elif GEOMETRY == SPHERICAL
      if(dir==JDIR)
        dl = dl*rt(i);
      else
        if(dir==KDIR)
          dl = dl*rt(i)*dmu(j)/dx2(j);
#endif

      invDt(k,j,i) = invDt(k,j,i) + HALF_F*(cMax(k+koffset,j+joffset,i+ioffset)
                    + cMax(k,j,i)) / (dl);

      if(haveParabolicTerms) {
        invDt(k,j,i) = invDt(k,j,i) + TWO_F* FMAX(dMax(k+koffset,j+joffset,i+ioffset),
                                                   dMax(k,j,i)) / (dl*dl);
      }

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
        rhs[MX1+meanDir] -= meanV*rhs[RHO];
        #if HAVE_ENERGY
          rhs[ENG] -= meanV * ( HALF_F*meanV*rhs[RHO] + rhs[MX1+meanDir]);
        #endif
      }

      // Potential terms
      if(needPotential) {
        real dphi;
        if (dir==IDIR) {
          // Gravitational force in direction i
          dphi = - 1.0/12.0 * (
                        - phiP(k,j,i+2) + 8.0 * phiP(k,j,i+1)
                        - 8.0*phiP(k,j,i-1) + phiP(k,j,i-2));
        }
        if (dir==JDIR) {
          // Gravitational force in direction j
          dphi = - 1.0/12.0 * (
                        - phiP(k,j+2,i) + 8.0 * phiP(k,j+1,i)
                        - 8.0*phiP(k,j-1,i) + phiP(k,j-2,i));
        }
        if (dir==KDIR) {
          // Gravitational force in direction k
          dphi = - 1.0/12.0 * (
                        - phiP(k+2,j,i) + 8.0 * phiP(k+1,j,i)
                        - 8.0*phiP(k-1,j,i) + phiP(k-2,j,i));
        }
        rhs[MX1+dir] += dt * Vc(RHO,k,j,i) * dphi /dl;

        #if HAVE_ENERGY
          // Add gravitational force work as a source term
          // This is equivalent to rho * v . nabla(phi)
          // (note that Flux has already been multiplied by A)
          rhs[ENG] += HALF_F * dtdV  *
                    (Flux(RHO,k,j,i) + Flux(RHO, k+koffset, j+joffset, i+ioffset)) * dphi;
        #endif
      }

      // Body force
      if(needBodyForce) {
        rhs[MX1+dir] += dt * Vc(RHO,k,j,i) * bodyForce(dir,k,j,i);
        #if HAVE_ENERGY
          //  rho * v . f, where rhov is taken as a  volume average of Flux(RHO)
          rhs[ENG] += HALF_F * dtdV * dl *
                        (Flux(RHO,k,j,i) + Flux(RHO, k+koffset, j+joffset, i+ioffset)) *
                         bodyForce(dir,k,j,i);
        #endif // HAVE_ENERGY

        // Particular cases if we do not sweep all of the components
        #if DIMENSIONS == 1 && COMPONENTS > 1
          EXPAND(                                                           ,
                    rhs[MX2] += dt * Vc(RHO,k,j,i) * bodyForce(JDIR,k,j,i);   ,
                    rhs[MX3] += dt * Vc(RHO,k,j,i) * bodyForce(KDIR,k,j,i);    )
          #if HAVE_ENERGY
            rhs[ENG] += dt * (EXPAND( ZERO_F                                                   ,
                                      + Vc(RHO,k,j,i) * Vc(VX2,k,j,i) * bodyForce(JDIR,k,j,i)  ,
                                      + Vc(RHO,k,j,i) * Vc(VX3,k,j,i) * bodyForce(KDIR,k,j,i) ));
          #endif
        #endif
        #if DIMENSIONS == 2 && COMPONENTS == 3
          // Only add this term once!
          if(dir==JDIR) {
            rhs[MX3] += dt * Vc(RHO,k,j,i) * bodyForce(KDIR,k,j,i);
            #if HAVE_ENERGY
              rhs[ENG] += dt * Vc(RHO,k,j,i) * Vc(VX3,k,j,i) * bodyForce(KDIR,k,j,i);
            #endif
          }
        #endif
      }


      // Evolve the field components
#pragma unroll
      for(int nv = 0 ; nv < NVAR ; nv++) {
        // Do not evolve the field components if they are computed by CT (i.e. if they are in Vs)
        D_EXPAND( if(nv == BX1) { continue; }  ,
                  if(nv == BX2) { continue; }  ,
                  if(nv == BX3) { continue; }  )


        Uc(nv,k,j,i) = Uc(nv,k,j,i) + rhs[nv];
      }
    }
  );

  idfx::popRegion();
}
#endif // HYDRO_CALCRIGHTHANDSIDE_HPP_
