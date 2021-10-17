// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2021 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include <string>

#include "idefix.hpp"
#include "hydro.hpp"
#include "dataBlock.hpp"
#include "fargo.hpp"

// If no high order fargo, then choose the order according to reconstruction
#ifndef HIGH_ORDER_FARGO
  #if ORDER == 3
  #define HIGH_ORDER_FARGO
  #endif
#endif

#ifdef HIGH_ORDER_FARGO
KOKKOS_INLINE_FUNCTION real FargoFlux(const IdefixArray4D<real> &Vin, int n, int k, int j, int i,
                                      int so, int ds, int sbeg, real eps) {
  // compute shifted indices, taking into account the fact that we're periodic
  int sop1 = so+1;
  if(sop1-sbeg >= ds) sop1 = sop1-ds;
  int sop2 = sop1+1;
  if(sop2-sbeg >= ds) sop2 = sop2-ds;

  int som1 = so-1;
  if(som1-sbeg< 0 ) som1 = som1+ds;
  int som2 = som1-1;
  if(som2-sbeg< 0 ) som2 = som2+ds;

  int sign = (eps>=0) ? 1 : -1;
  real dqm2,dqm1,dqp1,dqp2, q0,qm1, qp1;
  #if GEOMETRY == CARTESIAN || GEOMETRY == POLAR
    q0 = Vin(n,k,so,i);
    qm1 = Vin(n,k,som1,i);
    qp1 = Vin(n,k,sop1,i);
    dqm2 = qm1 - Vin(n,k,som2,i);
    dqm1 = q0 - qm1;
    dqp1 = qp1 - q0;
    dqp2 = Vin(n,k,sop2,i) - qp1;
  #elif GEOMETRY == SPHERICAL
    q0 = Vin(n,so,j,i);
    qm1 = Vin(n,som1,j,i);
    qp1 = Vin(n,sop1,j,i);
    dqm2 = qm1 - Vin(n,som2,j,i);
    dqm1 = q0 - qm1;
    dqp1 = qp1 - q0;
    dqp2 = Vin(n,sop2,j,i) - qp1;
  #endif
    // slope limited values around the reference point
    real dqlm = (dqm2*dqm1 > ZERO_F ? TWO_F*dqm2*dqm1/(dqm2 + dqm1) : ZERO_F);
    real dql0 = (dqm1*dqp1 > ZERO_F ? TWO_F*dqm1*dqp1/(dqm1 + dqp1) : ZERO_F);
    real dqlp = (dqp2*dqp1 > ZERO_F ? TWO_F*dqp2*dqp1/(dqp2 + dqp1) : ZERO_F);

    real dqp = 0.5 * dqp1 - (dqlp - dql0) / 6.0;
    real dqm = -0.5 * dqm1 - (dql0 - dqlm) / 6.0;

    if(dqp*dqm>0.0) {
       dqp = dqm = 0.0;
    } else {
      if(FABS(dqp) >= 2.0*FABS(dqm)) dqp = -2.0*dqm;
      if(FABS(dqm) >= 2.0*FABS(dqp)) dqm = -2.0*dqp;
    }

    real qp = q0 + dqp;
    real qm = q0 + dqm;

    real dqc = dqp - dqm;
    real d2q = dqp + dqm;

    real F;
    if(eps > 0.0) {
      F = eps*(qp - 0.5*eps*(dqc + d2q*(3.0 - 2.0*eps)));
    } else {
      F = eps*(qm - 0.5*eps*(dqc - d2q*(3.0 + 2.0*eps)));
    }

  return(F);
}

#else// HIGH_ORDER_FARGO
KOKKOS_INLINE_FUNCTION real FargoFlux(const IdefixArray4D<real> &Vin, int n, int k, int j, int i,
                                      int so, int ds, int sbeg, real eps) {
  // compute shifted indices, taking into account the fact that we're periodic
  int sop1 = so+1;
  if(sop1-sbeg >= ds) sop1 = sop1-ds;
  int som1 = so-1;
  if(som1-sbeg< 0 ) som1 = som1+ds;

  int sign = (eps>=0) ? 1 : -1;
  real F, dqm, dqp, dq, V0;
  #if GEOMETRY == CARTESIAN || GEOMETRY == POLAR
    V0 = Vin(n,k,so,i);
    dqm = V0 - Vin(n,k,som1,i);
    dqp = Vin(n,k,sop1,i) - V0;
  #elif GEOMETRY == SPHERICAL
    V0 = Vin(n,so,j,i);
    dqm = V0 - Vin(n,som1,j,i);
    dqp = Vin(n,sop1,j,i) - V0;
  #endif
    dq = (dqp*dqm > ZERO_F ? TWO_F*dqp*dqm/(dqp + dqm) : ZERO_F);
    F = eps*(V0 + sign*0.5*dq*(1.0-sign*eps));
  return(F);
}

#endif // HIGH_ORDER_FARGO



void Fargo::Init(Input &input, Grid &grid, Hydro *hydro) {
  idfx::pushRegion("Fargo::Init");
  // Todo(lesurg): should work on the shearing box version of this.
  this->hydro = hydro;
  this->data = hydro->data;
  if(input.CheckEntry("Hydro","fargo")>=0) {
    std::string opType = input.GetString("Hydro","fargo",0);
    if(opType.compare("userdef")==0) {
      this->type=userdef;
    } else if(opType.compare("shearingbox")==0) {
      this->type=shearingbox;
      #if GEOMETRY != CARTESIAN
        IDEFIX_ERROR("Fargo+shearingbox only compatible with cartesian geometry");
      #endif
    } else {
      IDEFIX_ERROR("Unknown fargo option in the input file. "
      "Only userdef and shearingbox are allowed");
    }
  } else {
    IDEFIX_ERROR("Fargo should be enabled in the input file with either userdef "
                  "or shearingbox options");
  }

  #if GEOMETRY == CARTESIAN || GEOMETRY == POLAR
    // Check there is no domain decomposition in the intended fargo direction
    if(data->mygrid->nproc[JDIR]>1) {
      IDEFIX_ERROR("Fargo is not yet compatible with MPI decomposition along the X2 direction");
    }
    if(this->type==userdef)
      this->meanVelocity = IdefixArray2D<real>("FargoVelocity",data->np_tot[KDIR],
                                                             data->np_tot[IDIR]);
  #elif GEOMETRY == SPHERICAL
    // Check there is no domain decomposition in the intended fargo direction
    if(data->mygrid->nproc[KDIR]>1) {
      IDEFIX_ERROR("Fargo is not yet compatible with MPI decomposition along the X3 direction");
    }
    this->meanVelocity = IdefixArray2D<real>("FargoVelocity",data->np_tot[JDIR],
                                                             data->np_tot[IDIR]);
  #else
    IDEFIX_ERROR("Fargo is not compatible with the GEOMETRY you intend to use");
  #endif

  // Initialise our scratch space
  this->scratch = IdefixArray4D<real>("FargoScratchSpace",NVAR,data->np_tot[KDIR]
                                                              ,data->np_tot[JDIR]
                                                              ,data->np_tot[IDIR]);
  if(type==userdef) {
    idfx::cout << "Fargo: Enabled with user-defined velocity function" << std::endl;
  } else if(type==shearingbox) {
    idfx::cout << "Fargo: Enabled with shearing-box velocity function" << std::endl;
  } else {
    IDEFIX_ERROR("Something went wrong during Fargo initialisation");
  }
  idfx::popRegion();
}

void Fargo::EnrollVelocity(FargoVelocityFunc myFunc) {
  if(this->type!=userdef) {
    IDEFIX_ERROR("Fargo velocity function enrollment requires Hydro/Fargo "
                 "to be set to userdef in .ini file");
  }
  this->fargoVelocityFunc = myFunc;
  idfx::cout << "Fargo: User-defined velocity function has been enrolled." << std::endl;
}

// This function fetches Fargo velocity when required
void Fargo::GetFargoVelocity(real t) {
  idfx::pushRegion("Fargo::GetFargoVelocity");
  if(type != userdef)
    IDEFIX_ERROR("Fargo::GetFargoVelocity should not be called when fargo != userdef");
  if(velocityHasBeenComputed == false) {
    if(fargoVelocityFunc== NULL) {
      IDEFIX_ERROR("No Fargo velocity function has been defined");
    }
    fargoVelocityFunc(*data, meanVelocity);
    velocityHasBeenComputed = true;
  }

  idfx::popRegion();
}

void Fargo::AddVelocity(const real t) {
  idfx::pushRegion("Fargo::AddVelocity");
  if(type==userdef) {
    GetFargoVelocity(t);
  }
  IdefixArray1D<real> x1 = data->x[IDIR];
  IdefixArray4D<real> Vc = hydro->Vc;
  IdefixArray2D<real> meanV = this->meanVelocity;
  FargoType fargoType = type;
  real sbS = hydro->sbS;

  idefix_for("FargoAddVelocity",
              0,data->np_tot[KDIR],
              0,data->np_tot[JDIR],
              0,data->np_tot[IDIR],
              KOKKOS_LAMBDA(int k, int j, int i) {
                #if GEOMETRY == CARTESIAN || GEOMETRY == POLAR
                  if(fargoType == userdef) {
                    Vc(VX2,k,j,i) += meanV(k,i);
                  } else if(fargoType == shearingbox) {
                    Vc(VX2,k,j,i) += sbS*x1(i);
                  }
                #elif GEOMETRY == SPHERICAL
                  Vc(VX3,k,j,i) += meanV(j,i);
                #endif
              });

  idfx::popRegion();
}

void Fargo::SubstractVelocity(const real t) {
  idfx::pushRegion("Fargo::SubstractVelocity");
  if(type==userdef) {
    GetFargoVelocity(t);
  }
  IdefixArray1D<real> x1 = data->x[IDIR];
  IdefixArray4D<real> Vc = hydro->Vc;
  IdefixArray2D<real> meanV = this->meanVelocity;
  FargoType fargoType = type;
  real sbS = hydro->sbS;

  idefix_for("FargoSubstractVelocity",
              0,data->np_tot[KDIR],
              0,data->np_tot[JDIR],
              0,data->np_tot[IDIR],
              KOKKOS_LAMBDA(int k, int j, int i) {
                #if GEOMETRY == CARTESIAN || GEOMETRY == POLAR
                  if(fargoType == userdef) {
                    Vc(VX2,k,j,i) -= meanV(k,i);
                  } else if(fargoType == shearingbox) {
                    Vc(VX2,k,j,i) -= sbS*x1(i);
                  }
                #elif GEOMETRY == SPHERICAL
                  Vc(VX3,k,j,i) -= meanV(j,i);
                #endif
              });

  idfx::popRegion();
}

void Fargo::ShiftSolution(const real t, const real dt) {
  idfx::pushRegion("Fargo::ShiftSolution");

  // Refresh the fargo velocity function
  if(type==userdef) {
    GetFargoVelocity(t);
  }

  IdefixArray4D<real> Uc = hydro->Uc;
  IdefixArray4D<real> scrh = this->scratch;
  IdefixArray2D<real> meanV = this->meanVelocity;
  IdefixArray1D<real> x1 = data->x[IDIR];
  IdefixArray1D<real> x2 = data->x[JDIR];
  IdefixArray1D<real> dx2 = data->dx[JDIR];
  IdefixArray1D<real> dx3 = data->dx[KDIR];
  IdefixArray1D<real> sinx2 = data->sinx2;
  IdefixArray1D<real> sinx2m = data->sinx2m;
  FargoType fargoType = type;
  real sbS = hydro->sbS;

  real Lphi;
  int sbeg, send;
  #if GEOMETRY == CARTESIAN || GEOMETRY == POLAR
    Lphi = data->mygrid->xend[JDIR] - data->mygrid->xbeg[JDIR];
    sbeg = data->beg[JDIR];
    send = data->end[JDIR];
  #elif GEOMETRY == SPHERICAL
    Lphi = data->mygrid->xend[KDIR] - data->mygrid->xbeg[KDIR];
    sbeg = data->beg[KDIR];
    send = data->end[KDIR];
  #else
    Lphi = 1.0;   // Do nothing, but initialize this.
  #endif

  idefix_for("Fargo:ShiftVc",
              0,NVAR,
              data->beg[KDIR],data->end[KDIR],
              data->beg[JDIR],data->end[JDIR],
              data->beg[IDIR],data->end[IDIR],
              KOKKOS_LAMBDA(int n, int k, int j, int i) {
                real w,dphi;
                int s;
                #if GEOMETRY == CARTESIAN
                 if(fargoType==userdef) {
                  w = meanV(k,i);
                 } else if(fargoType==shearingbox) {
                  w = sbS*x1(i);
                 }
                 dphi = dx2(j);
                 s = j;
                #elif GEOMETRY == POLAR
                 w = meanV(k,i)/x1(i);
                 dphi = dx2(j);
                 s = j;
                #elif GEOMETRY == SPHERICAL
                 w = meanV(j,i)/(x1(i)*sinx2(j));
                 dphi = dx3(k);
                 s = k;
                #endif

                // Compute the offset in phi, modulo the full domain size
                real dL = std::fmod(w*dt, Lphi);

                // Translate this into # of cells
                int m = static_cast<int> (std::floor(dL/dphi+HALF_F));

                // get the remainding shift
                real eps = dL/dphi - m;

                // origin index before the shift
                // Note the trick to get a positive module i%%n = (i%n + n)%n;
                int ds = send-sbeg;

                // so is the "origin" index
                int so = sbeg + ((s-m-sbeg)%ds+ds)%ds;

                // Define Left and right fluxes
                // Fluxes are defined from slop-limited interpolation
                // Using Van-leer slope limiter (consistently with the main advection scheme)
                real Fl,Fr;

                if(eps>=ZERO_F) {
                  int som1 = so-1;
                  if(som1-sbeg< 0 ) som1 = som1+ds;
                  Fl = FargoFlux(Uc, n, k, j, i, som1, ds, sbeg, eps);
                  Fr = FargoFlux(Uc, n, k, j, i, so, ds, sbeg, eps);
                } else {
                  int sop1 = so+1;
                  if(sop1-sbeg >= ds) sop1 = sop1-ds;
                  Fl = FargoFlux(Uc, n, k, j, i, so, ds, sbeg, eps);
                  Fr = FargoFlux(Uc, n, k, j, i, sop1, ds, sbeg, eps);
                }

                #if GEOMETRY == CARTESIAN || GEOMETRY == POLAR
                  scrh(n,k,s,i) = Uc(n,k,so,i) - (Fr - Fl);
                #elif GEOMETRY == SPHERICAL
                  scrh(n,s,j,i) = Uc(n,so,j,i) - (Fr - Fl);
                #endif
              });

  // move back scratch space into Uc
  Kokkos::deep_copy(Uc,scrh);

#if MHD == YES
  IdefixArray4D<real> Vs = hydro->Vs;
  IdefixArray3D<real> ex = hydro->emf.Ex1;
  IdefixArray3D<real> ey = hydro->emf.Ex2;
  IdefixArray3D<real> ez = hydro->emf.Ex3;
  IdefixArray1D<real> x1m = data->xl[IDIR];
  IdefixArray1D<real> x2m = data->xl[JDIR];
  IdefixArray1D<real> dmu = data->dmu;
  IdefixArray1D<real> dx1 = data->dx[IDIR];



  #if GEOMETRY == CARTESIAN || GEOMETRY == POLAR
    IdefixArray3D<real> ek = ez;
  #elif GEOMETRY == SPHERICAL
    IdefixArray3D<real> ek = ey;
  #endif

  idefix_for("Fargo:ComputeEk",
    data->beg[KDIR],data->end[KDIR]+KOFFSET,
    data->beg[JDIR],data->end[JDIR]+JOFFSET,
    data->beg[IDIR],data->end[IDIR]+IOFFSET,
    KOKKOS_LAMBDA(int k, int j, int i) {
      real w,dphi;
      int s;
      #if GEOMETRY == CARTESIAN
        if(fargoType==userdef) {
         w = 0.5*(meanV(k,i-1)+meanV(k,i));
        } else if(fargoType==shearingbox) {
         w = sbS*x1m(i);
        }
        dphi = dx2(j);
        s = j;
      #elif GEOMETRY == POLAR
        w = 0.5*(meanV(k,i-1)+meanV(k,i))/x1m(i);
        dphi = dx2(j);
        s = j;
      #elif GEOMETRY == SPHERICAL
        w = 0.5*(meanV(j,i-1)/x1(i-1)+meanV(j,i)/x1(i))/sinx2(j);
        dphi = dx3(k);
        s = k;
      #endif

      // Compute the offset in phi, modulo the full domain size
      real dL = std::fmod(w*dt, Lphi);

      // Translate this into # of cells
      int m = static_cast<int> (std::floor(dL/dphi+HALF_F));

      // get the remainding shift
      real eps = dL/dphi - m;

      // origin index before the shift
      // Note the trick to get a positive module i%%n = (i%n + n)%n;
      int n = send-sbeg;

      // so is the "origin" index
      int so = sbeg + ((s-m-sbeg)%n+n)%n;

      #if GEOMETRY == CARTESIAN || GEOMETRY == POLAR
        if(eps>=ZERO_F) {
          int som1 = sbeg + ((so-1-sbeg)%n+n)%n;
          ek(k,s,i) = FargoFlux(Vs, BX1s, k, j, i, som1, n, sbeg, eps);

        } else {
          ek(k,s,i) = FargoFlux(Vs, BX1s, k, j, i, so, n, sbeg, eps);
        }
        if(m>0) {
          for(int ss = s-m ; ss < s ; ss++) {
            int sc = sbeg + ((ss-sbeg)%n+n)%n;
            ek(k,s,i) += Vs(BX1s,k,sc,i);
          }
        } else {
          for(int ss = s ; ss < s-m ; ss++) {
            int sc = sbeg + ((ss-sbeg)%n+n)%n;
            ek(k,s,i) -= Vs(BX1s,k,sc,i);
          }
        }
      #elif GEOMETRY == SPHERICAL
        if(eps>=ZERO_F) {
          int som1 = sbeg + ((so-1-sbeg)%n+n)%n;
          ek(s,j,i) = FargoFlux(Vs, BX1s, k, j, i, som1, n, sbeg, eps);
        } else {
          ek(s,j,i) = FargoFlux(Vs, BX1s, k, j, i, so, n, sbeg, eps);
        }
        if(m>0) {
          for(int ss = s-m ; ss < s ; ss++) {
            int sc = sbeg + ((ss-sbeg)%n+n)%n;
            ek(s,j,i) += Vs(BX1s,sc,j,i);
          }
        } else {
          for(int ss = s ; ss < s-m ; ss++) {
            int sc = sbeg + ((ss-sbeg)%n+n)%n;
            ek(s,j,i) -= Vs(BX1s,sc,j,i);
          }
        }
      #endif  // GEOMETRY

      ek(k,j,i) *= dphi;
    });

#if DIMENSIONS == 3
  // In cartesian and polar coordinates, ei is actually -Ex
  IdefixArray3D<real> ei = ex;

  idefix_for("Fargo:ComputeEi",
    data->beg[KDIR],data->end[KDIR]+KOFFSET,
    data->beg[JDIR],data->end[JDIR]+JOFFSET,
    data->beg[IDIR],data->end[IDIR]+IOFFSET,
    KOKKOS_LAMBDA(int k, int j, int i) {
      real w,dphi;
      int s;
      #if GEOMETRY == CARTESIAN
        if(fargoType==userdef) {
         w = 0.5*(meanV(k,i)+meanV(k-1,i));
        } else if(fargoType==shearingbox) {
         w = sbS*x1(i);
        }

        dphi = dx2(j);
        s = j;
      #elif GEOMETRY == POLAR
        w = 0.5*(meanV(k-1,i)+meanV(k,i))/x1(i);
        dphi = dx2(j);
        s = j;
      #elif GEOMETRY == SPHERICAL
        w = 0.5*(meanV(j-1,i)/sinx2(j-1)+meanV(j,i)/sinx2(j))/(x1(i));
        dphi = dx3(k);
        s = k;
      #endif

      // Compute the offset in phi, modulo the full domain size
      real dL = std::fmod(w*dt, Lphi);

      // Translate this into # of cells
      int m = static_cast<int> (std::floor(dL/dphi+HALF_F));

      // get the remainding shift
      real eps = dL/dphi - m;

      // origin index before the shift
      // Note the trick to get a positive module i%%n = (i%n + n)%n;
      int n = send-sbeg;

      // so is the "origin" index
      int so = sbeg + ((s-m-sbeg)%n+n)%n;

      // compute shifted indices, taking into account the fact that we're periodic
      int sop1 = sbeg + (so+1-sbeg)%(send-sbeg);
      int som1 = sbeg + (so-1-sbeg)%(send-sbeg);
      int som2 = sbeg + (so-2-sbeg)%(send-sbeg);

      // Compute EMF due to the shift via second order reconstruction
      real dqm, dqp, dq;

      #if GEOMETRY == CARTESIAN || GEOMETRY == POLAR
        if(eps>=ZERO_F) {
          int som1 = sbeg + ((so-1-sbeg)%n+n)%n;
          ei(k,s,i) = FargoFlux(Vs, BX3s, k, j, i, som1, n, sbeg, eps);
        } else {
          ei(k,s,i) = FargoFlux(Vs, BX3s, k, j, i, so, n, sbeg, eps);
        }
        if(m>0) {
          for(int ss = s-m ; ss < s ; ss++) {
            int sc = sbeg + ((ss-sbeg)%n+n)%n;
            ei(k,s,i) += Vs(BX3s,k,sc,i);
          }
        } else {
          for(int ss = s ; ss < s-m ; ss++) {
            int sc = sbeg + ((ss-sbeg)%n+n)%n;
            ei(k,s,i) -= Vs(BX3s,k,sc,i);
          }
        }
      #elif GEOMETRY == SPHERICAL
      if(eps>=ZERO_F) {
        int som1 = sbeg + ((so-1-sbeg)%n+n)%n;
        ei(s,j,i) = FargoFlux(Vs, BX2s, k, j, i, som1, n, sbeg, eps);
      } else {
        ei(s,j,i) = FargoFlux(Vs, BX2s, k, j, i, so, n, sbeg, eps);
      }
      if(m>0) {
        for(int ss = s-m ; ss < s ; ss++) {
          int sc = sbeg + ((ss-sbeg)%n+n)%n;
          ei(s,j,i) += Vs(BX2s,sc,j,i);
        }
      } else {
        for(int ss = s ; ss < s-m ; ss++) {
          int sc = sbeg + ((ss-sbeg)%n+n)%n;
          ei(s,j,i) -= Vs(BX2s,sc,j,i);
        }
      }

      #endif  // GEOMETRY

      ei(k,j,i) *= dphi;
    });
#endif

  // Update field components according to the computed EMFS

  idefix_for("Fargo::EvolvMagField",
             data->beg[KDIR],data->end[KDIR]+KOFFSET,
             data->beg[JDIR],data->end[JDIR]+JOFFSET,
             data->beg[IDIR],data->end[IDIR]+IOFFSET,
    KOKKOS_LAMBDA (int k, int j, int i) {
      real rhsx1, rhsx2, rhsx3;

#if GEOMETRY == CARTESIAN || GEOMETRY == POLAR
      Vs(BX1s,k,j,i) += -  (ez(k,j+1,i) - ez(k,j,i) ) / dx2(j);

#elif GEOMETRY == SPHERICAL
      Vs(BX1s,k,j,i) += -  sinx2(j)*dx2(j)/dmu(j)*(ey(k+1,j,i) - ey(k,j,i) ) / dx3(k);
#endif

#if GEOMETRY == CARTESIAN
      Vs(BX2s,k,j,i) += D_EXPAND(    0.0                          ,
                            + (ez(k,j,i+1) - ez(k,j,i)) / dx1(i)  ,
                            + (ex(k+1,j,i) - ex(k,j,i)) / dx3(k)  );
#elif GEOMETRY == POLAR
      Vs(BX2s,k,j,i) += D_EXPAND(    0.0                                          ,
                            + (x1m(i+1) * ez(k,j,i+1) - x1m(i)*ez(k,j,i)) / dx1(i)  ,
                            + x1(i) *  (ex(k+1,j,i) - ex(k,j,i)) / dx3(k)  );
#elif GEOMETRY == SPHERICAL
    #if DIMENSIONS == 3
      Vs(BX2s,k,j,i) += - (ex(k+1,j,i) - ex(k,j,i)) / dx3(k);
    #endif
#endif // GEOMETRY

#if DIMENSIONS == 3
  #if GEOMETRY == CARTESIAN || GEOMETRY == POLAR
    Vs(BX3s,k,j,i) -=  (ex(k,j+1,i) - ex(k,j,i) ) / dx2(j);
  #elif GEOMETRY == SPHERICAL
    real A1p = x1m(i+1)*x1m(i+1);
    real A1m = x1m(i)*x1m(i);
    real A2m = FABS(sinx2m(j));
    real A2p = FABS(sinx2m(j+1));
    Vs(BX3s,k,j,i) += sinx2(j) * (A1p * ey(k,j,i+1) - A1m * ey(k,j,i))/(x1(i)*dx1(i))
                      + (A2p * ex(k,j+1,i) - A2m * ex(k,j,i))/dx2(j);
  #endif
#endif// DIMENSIONS
  });

  // Rebuild the cell-centered field components
  this->hydro->boundary.ReconstructVcField(Uc);

#endif // MHD

  idfx::popRegion();
}
