// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2021 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include "hydro.hpp"
#include "dataBlock.hpp"


// Compute the right handside in direction dir from conservative equation, with timestep dt
void Hydro::CalcRightHandSide(int dir, real t, real dt) {
  idfx::pushRegion("Hydro::CalcRightHandSide");

  IdefixArray4D<real> Uc   = this->Uc;
  IdefixArray4D<real> Vc   = this->Vc;
  IdefixArray4D<real> Flux = this->FluxRiemann;
  IdefixArray3D<real> A    = data->A[dir];
  IdefixArray3D<real> dV   = data->dV;
  IdefixArray1D<real> x1m  = data->xl[IDIR];
  IdefixArray1D<real> x1   = data->x[IDIR];
  IdefixArray1D<real> sm   = data->sm;
  IdefixArray1D<real> rt   = data->rt;
  IdefixArray1D<real> dmu  = data->dmu;
  IdefixArray1D<real> s    = data->s;
  IdefixArray1D<real> dx   = data->dx[dir];
  IdefixArray1D<real> dx2  = data->dx[JDIR];
  IdefixArray3D<real> invDt = this->InvDt;
  IdefixArray3D<real> cMax = this->cMax;
  IdefixArray3D<real> dMax = this->dMax;
  IdefixArray4D<real> viscSrc = this->viscosity.viscSrc;

  // Gravitational potential
  IdefixArray3D<real> phiP = this->phiP;
  bool needPotential = this->haveGravPotential;

  // parabolic terms
  bool haveParabolicTerms = this->haveParabolicTerms;

  // Viscosity
  bool haveViscosity = this->haveViscosity;

  if(needPotential) {
    IdefixArray1D<real> x1,x2,x3;

    if(dir==IDIR)
      x1 = data->xl[IDIR];
    else
      x1 = data->x[IDIR];

    if(dir==JDIR)
      x2 = data->xl[JDIR];
    else
      x2 = data->x[JDIR];

    if(dir==KDIR)
      x3 = data->xl[KDIR];
    else
      x3 = data->x[KDIR];

    if(this->gravPotentialFunc == nullptr)
      IDEFIX_ERROR("Gravitational potential is enabled, "
                   "but no user-defined potential has been enrolled.");

    gravPotentialFunc(*data, t, x1, x2, x3, phiP);
  }


  int ioffset,joffset,koffset;
  ioffset=joffset=koffset=0;
  // Determine the offset along which we do the extrapolation
  if(dir==IDIR) ioffset=1;
  if(dir==JDIR) joffset=1;
  if(dir==KDIR) koffset=1;


  idefix_for("CalcTotalFlux",
             data->beg[KDIR],data->end[KDIR]+koffset,
             data->beg[JDIR],data->end[JDIR]+joffset,
             data->beg[IDIR],data->end[IDIR]+ioffset,
    KOKKOS_LAMBDA (int k, int j, int i) {
      // TODO(lesurg): Should add gravitational potential here and Fargo source terms when needed
#if HAVE_ENERGY
      if(needPotential)
        Flux(ENG, k, j, i) += Flux(RHO, k, j, i) * phiP(k,j,i);  // Potential at the cell face
#endif

      real Ax = A(k,j,i);

#if GEOMETRY != CARTESIAN
      if(Ax<SMALL_NUMBER)
        Ax=SMALL_NUMBER;    // Essentially to avoid singularity around poles
#endif

      for(int nv = 0 ; nv < NVAR ; nv++) {
        Flux(nv,k,j,i) = Flux(nv,k,j,i) * Ax;
        /*if(Flux(nv,k,j,i)!=Flux(nv,k,j,i)) {
          printf("Shit Flux at dir=%d var=%d (%d,%d)\n",dir,nv,i,j);
          exit(1);
        }*/
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
        Flux(iMPHI,k,j,i) = Flux(iMPHI,k,j,i) * FABS(sm(j));
    #if MHD == YES
        Flux(iBPHI,k,j,i) = Flux(iBPHI,k,j,i)  / Ax;
    #endif // MHD
  #endif // COMPONENTS = 3
      }
#endif // GEOMETRY == SPHERICAL
    }
  );


  idefix_for("CalcRightHandSide",
             data->beg[KDIR],data->end[KDIR],
             data->beg[JDIR],data->end[JDIR],
             data->beg[IDIR],data->end[IDIR],
    KOKKOS_LAMBDA (int k, int j, int i) {
      //
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
        rhs[iMPHI] /= FABS(s(j));
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
                    + cMax(k,j,i)) / dl;

      if(haveParabolicTerms) {
        invDt(k,j,i) = invDt(k,j,i) + TWO_F * FMAX(dMax(k+koffset,j+joffset,i+ioffset),
                                                   dMax(k,j,i)) / (dl*dl);
      }


      // Potential terms
      if(needPotential) {
        switch(dir) {
          case IDIR:
            // Gravitational force in direction i
            rhs[MX1] -= dt/dl * Vc(RHO,k,j,i) * (phiP(k,j,i+1) - phiP(k,j,i));
            break;
#if DIMENSIONS >=2
          case JDIR:
            // Gravitational force in direction j
            rhs[MX2] -= dt/dl * Vc(RHO,k,j,i) * (phiP(k,j+1,i) - phiP(k,j,i));
            break;
#endif
#if DIMENSIONS == 3
          case KDIR:
            // Gravitational force in direction k
            rhs[MX3] -= dt/dl * Vc(RHO,k,j,i) * (phiP(k+1,j,i) - phiP(k,j,i));
            break;
#endif
        }

#if HAVE_ENERGY
        // We conserve total energy without potential
        rhs[ENG] -=  HALF_F * (phiP(k+koffset,j+joffset,i+ioffset) + phiP(k,j,i)) * rhs[RHO];
#endif
      }

      // Evolve the field components
#pragma unroll
      for(int nv = 0 ; nv < NVAR ; nv++) {
        // Do not evolve the field components if they are computed by CT (i.e. if they are in Vs)
        D_EXPAND( if(nv == BX1) { continue; }  ,
                  if(nv == BX2) { continue; }  ,
                  if(nv == BX3) { continue; }  )

/*
          if(rhs[nv]!=rhs[nv]) {
              printf("Wrong Riemann RHS for var %d in dir=%d\n",nv,dir);
              exit(1);
          }
*/
        Uc(nv,k,j,i) = Uc(nv,k,j,i) + rhs[nv];
      }
    }
  );

  idfx::popRegion();
}
