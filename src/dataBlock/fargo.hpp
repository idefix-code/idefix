// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef DATABLOCK_FARGO_HPP_
#define DATABLOCK_FARGO_HPP_

#include <vector>
#include "idefix.hpp"
#ifdef WITH_MPI
  #include "mpi.hpp"
#endif

#include "physics.hpp"
#include "slopeLimiter.hpp"

// Forward class hydro declaration
template <typename Phys> class Fluid;
using Hydro = Fluid<DefaultPhysics>;
class DataBlock;

using FargoVelocityFunc = void (*) (DataBlock &, IdefixArray2D<real> &);



class Fargo {
 public:
  enum FargoType {none, userdef, shearingbox};
  Fargo(Input &, int, DataBlock*);  // Initialisation
  void ShiftSolution(const real t, const real dt);  // Effectively shift the solution
  void SubstractVelocity(const real);
  void AddVelocity(const real);
  void EnrollVelocity(FargoVelocityFunc);
  void CheckMaxDisplacement();
  void ShowConfig();

// For internal use
  template <typename Phys>
  void AddVelocityFluid(const real, Fluid<Phys>* );

  template <typename Phys>
  void SubstractVelocityFluid(const real, Fluid<Phys>* );

  template <typename Phys>
  void ShiftFluid(const real t, const real dt, Fluid<Phys>* );

  template <typename Phys>
  void StoreToScratch(Fluid<Phys>*);

  void GetFargoVelocity(real);

  IdefixArray2D<real> meanVelocity;
  FargoType type{none};                 // By default, Fargo is disabled

 private:
  friend Hydro;
  DataBlock *data;

  IdefixArray4D<real> scrhUc;
  IdefixArray4D<real> scrhVs;

#ifdef WITH_MPI
  Mpi mpi;                      // Fargo-specific MPI layer
#endif

  std::array<int,3> beg;
  std::array<int,3> end;
  std::array<int,3> nghost;
  int maxShift;                         //< maximum number of cells along which we plan to shift.
  real dtMax{0};                        //< Maximum allowable dt for a given Fargo velocity
                                        //< when domain decomposition is enabled
  bool velocityHasBeenComputed{false};
  bool haveDomainDecomposition{false};

  FargoVelocityFunc fargoVelocityFunc{NULL};  // The user-defined fargo velocity function
};


// If no high order fargo, then choose the order according to reconstruction
#ifndef HIGH_ORDER_FARGO
  #if ORDER >= 3
  #define HIGH_ORDER_FARGO
  #endif
#endif

#ifdef HIGH_ORDER_FARGO
KOKKOS_INLINE_FUNCTION real FargoFlux(const IdefixArray4D<real> &Vin, int n, int k, int j, int i,
                                      int so, int ds, int sbeg, real eps,
                                      bool haveDomainDecomposition) {
  // compute shifted indices, taking into account the fact that we're periodic
  int sop1 = so+1;
  if(!haveDomainDecomposition && (sop1-sbeg >= ds)) sop1 = sop1-ds;
  int sop2 = sop1+1;
  if(!haveDomainDecomposition && (sop2-sbeg >= ds)) sop2 = sop2-ds;

  int som1 = so-1;
  if(!haveDomainDecomposition && (som1-sbeg< 0 )) som1 = som1+ds;
  int som2 = som1-1;
  if(!haveDomainDecomposition && (som2-sbeg< 0 )) som2 = som2+ds;

  real q0,qm1, qp1, qm2, qp2;
  #if GEOMETRY == CARTESIAN || GEOMETRY == POLAR
    q0 = Vin(n,k,so,i);
    qm1 = Vin(n,k,som1,i);
    qp1 = Vin(n,k,sop1,i);
    qm2 = Vin(n,k,som2,i);
    qp2 = Vin(n,k,sop2,i);
  #elif GEOMETRY == SPHERICAL
    q0 = Vin(n,so,j,i);
    qm1 = Vin(n,som1,j,i);
    qp1 = Vin(n,sop1,j,i);
    qm2 = Vin(n,som2,j,i);
    qp2 = Vin(n,sop2,j,i);
  #endif

  real qp, qm;
  SlopeLimiter<>::getPPMStates(qm2,qm1,q0,qp1,qp2, qm, qp);

  real dqp = qp-q0;
  real dqm = qm-q0;

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
                                      int so, int ds, int sbeg, real eps,
                                      bool haveDomainDecomposition) {
  // compute shifted indices, taking into account the fact that we're periodic
  int sop1 = so+1;
  if(!haveDomainDecomposition && (sop1-sbeg >= ds)) sop1 = sop1-ds;
  int som1 = so-1;
  if(!haveDomainDecomposition && (som1-sbeg< 0 )) som1 = som1+ds;
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

KOKKOS_INLINE_FUNCTION int modPositive(int x, int divisor) {
  int m = x % divisor;
  return m + ((m >> 31) & divisor); // equivalent to m + (m < 0 ? divisor : 0);
}


template <typename Phys>
void Fargo::AddVelocityFluid(const real t, Fluid<Phys>* hydro ) {
  idfx::pushRegion("Fargo::AddVelocityFluid");
  if(type==userdef) {
    GetFargoVelocity(t);
  }
  IdefixArray1D<real> x1 = data->x[IDIR];
  IdefixArray4D<real> Vc = hydro->Vc;
  IdefixArray2D<real> meanV = this->meanVelocity;
  [[maybe_unused]] FargoType fargoType = type;
  [[maybe_unused]] real sbS = hydro->sbS;
  #if GEOMETRY == CARTESIAN || GEOMETRY == POLAR
    constexpr int Vadv = VX2;
  #elif GEOMETRY == SPHERICAL
    constexpr int Vadv = VX3;
  #endif
  if constexpr (Phys::nvar > Vadv) {
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
  } // if constexpr
  idfx::popRegion();
}

template <typename Phys>
void Fargo::SubstractVelocityFluid(const real t, Fluid<Phys>* hydro) {
  idfx::pushRegion("Fargo::SubstractVelocityFluid");
  if(type==userdef) {
    GetFargoVelocity(t);
  }
  IdefixArray1D<real> x1 = data->x[IDIR];
  IdefixArray4D<real> Vc = hydro->Vc;
  [[maybe_unused]] IdefixArray2D<real> meanV = this->meanVelocity;
  [[maybe_unused]] FargoType fargoType = type;
  [[maybe_unused]] real sbS = hydro->sbS;
  #if GEOMETRY == CARTESIAN || GEOMETRY == POLAR
    constexpr int Vadv = VX2;
  #elif GEOMETRY == SPHERICAL
    constexpr int Vadv = VX3;
  #endif
  if constexpr (Phys::nvar > Vadv) {
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
    }
  idfx::popRegion();
}

template<typename Phys>
void Fargo::StoreToScratch(Fluid<Phys>* hydro) {
  IdefixArray4D<real> Uc = hydro->Uc;
  IdefixArray4D<real> scrhUc = this->scrhUc;
  bool haveDomainDecomposition = this->haveDomainDecomposition;
  int maxShift = this->maxShift;

  idefix_for("Fargo:StoreUc",
            0,Phys::nvar+hydro->nTracer,
            data->beg[KDIR],data->end[KDIR],
            data->beg[JDIR],data->end[JDIR],
            data->beg[IDIR],data->end[IDIR],
            KOKKOS_LAMBDA(int n, int k, int j, int i) {
              if(!haveDomainDecomposition) {
                scrhUc(n,k,j,i) = Uc(n,k,j,i);
              } else {
                #if GEOMETRY==POLAR || GEOMETRY==CARTESIAN
                  scrhUc(n,k,j+maxShift,i) = Uc(n,k,j,i);
                #elif GEOMETRY == SPHERICAL
                  scrhUc(n,k+maxShift,j,i) = Uc(n,k,j,i);
                #endif
              }
            });

  if constexpr(Phys::mhd) {
    #ifdef EVOLVE_VECTOR_POTENTIAL
      // Update Vs to its latest
      hydro->emf->ComputeMagFieldFromA(hydro->Ve,hydro->Vs);
    #endif
    // in MHD mode, we need to copy Vs only when there is domain decomposition, otherwise,
    // we just make a reference (this is already done by init)
    if(haveDomainDecomposition) {
      IdefixArray4D<real> Vs = hydro->Vs;
      IdefixArray4D<real> scrhVs = this->scrhVs;
      idefix_for("Fargo:StoreVs",
              0,DIMENSIONS,
              data->beg[KDIR],data->end[KDIR]+KOFFSET,
              data->beg[JDIR],data->end[JDIR]+JOFFSET,
              data->beg[IDIR],data->end[IDIR]+IOFFSET,
              KOKKOS_LAMBDA(int n, int k, int j, int i) {
                  #if GEOMETRY==POLAR || GEOMETRY==CARTESIAN
                    scrhVs(n,k,j+maxShift,i) = Vs(n,k,j,i);
                  #elif GEOMETRY == SPHERICAL
                    scrhVs(n,k+maxShift,j,i) = Vs(n,k,j,i);
                  #endif
              });
    }
  }
  #if WITH_MPI
    if(haveDomainDecomposition) {
      #if GEOMETRY == CARTESIAN || GEOMETRY == POLAR
        this->mpi.ExchangeX2(scrhUc, scrhVs);
      #elif GEOMETRY == SPHERICAL
        this->mpi.ExchangeX3(scrhUc, scrhVs);
      #endif
    }
  #endif
}

template<typename Phys>
void Fargo::ShiftFluid(const real t, const real dt, Fluid<Phys>* hydro) {
  idfx::pushRegion("Fargo::ShiftFluid");

  #if GEOMETRY == CYLINDRICAL
    IDEFIX_ERROR("Fargo is not compatible with cylindrical geometry "
                 "(which is intended to be 2D axisymmetric)");
  #else

  // Refresh the fargo velocity function
  if(type==userdef) {
    GetFargoVelocity(t);
  }
  if(haveDomainDecomposition && dt>dtMax) {
    std::stringstream message;
    message << "Your dt is too large with your domain decomposition and Fargo." << std::endl
            << "Got dt=" << dt << " and Fargo:dtmax=" << dtMax << "." << std::endl
            << "Try to increase [Fargo]:maxShift to a value larger than "
            << static_cast<int>(ceil(dt/(dtMax/maxShift))) << std::endl;
    IDEFIX_ERROR(message);
  }

  IdefixArray4D<real> Uc = hydro->Uc;
  IdefixArray4D<real> scrh = this->scrhUc;
  IdefixArray2D<real> meanV = this->meanVelocity;
  IdefixArray1D<real> x1 = data->x[IDIR];
  IdefixArray1D<real> x2 = data->x[JDIR];
  IdefixArray1D<real> dx2 = data->dx[JDIR];
  IdefixArray1D<real> dx3 = data->dx[KDIR];
  IdefixArray1D<real> sinx2 = data->sinx2;
  IdefixArray1D<real> sinx2m = data->sinx2m;
  [[maybe_unused]] FargoType fargoType = type;
  [[maybe_unused]] real sbS = hydro->sbS;
  bool haveDomainDecomposition = this->haveDomainDecomposition;
  int maxShift = this->maxShift;

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

  // move Uc to scratch, and fill the ghost zones if required.
  StoreToScratch(hydro);

  idefix_for("Fargo:ShiftVc",
              0,Phys::nvar+hydro->nTracer,
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
                int so;
                if(haveDomainDecomposition) {
                  so = s-m + maxShift;    // maxshift corresponds to the offset between
                                          // the indices in scrh and in Uc
                } else {
                  so = sbeg + modPositive(s-m-sbeg, ds);
                }

                // Define Left and right fluxes
                // Fluxes are defined from slope-limited interpolation
                // Using Van-leer slope limiter (consistently with the main advection scheme)
                real Fl,Fr;

                if(eps>=ZERO_F) {
                  int som1 = so-1;
                  if(!haveDomainDecomposition && som1-sbeg< 0 ) som1 = som1+ds;
                  Fl = FargoFlux(scrh, n, k, j, i, som1, ds, sbeg, eps, haveDomainDecomposition);
                  Fr = FargoFlux(scrh, n, k, j, i, so, ds, sbeg, eps, haveDomainDecomposition);
                } else {
                  int sop1 = so+1;
                  if(!haveDomainDecomposition && sop1-sbeg >= ds) sop1 = sop1-ds;
                  Fl = FargoFlux(scrh, n, k, j, i, so, ds, sbeg, eps, haveDomainDecomposition);
                  Fr = FargoFlux(scrh, n, k, j, i, sop1, ds, sbeg, eps, haveDomainDecomposition);
                }

                #if GEOMETRY == CARTESIAN || GEOMETRY == POLAR
                  Uc(n,k,s,i) = scrh(n,k,so,i) - (Fr - Fl);
                #elif GEOMETRY == SPHERICAL
                  Uc(n,s,j,i) = scrh(n,so,j,i) - (Fr - Fl);
                #endif
              });

  if constexpr(Phys::mhd) {
    IdefixArray4D<real> scrhVs = this->scrhVs;
    IdefixArray3D<real> ex = hydro->emf->Ex1;
    IdefixArray3D<real> ey = hydro->emf->Ex2;
    IdefixArray3D<real> ez = hydro->emf->Ex3;
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
        int so;
        if(haveDomainDecomposition) {
          so = s-m + maxShift;    // maxshift corresponds to the offset between
                                  // the indices in scrh and in Uc
        } else {
          so = sbeg + modPositive(s-m-sbeg,n);
        }

        #if GEOMETRY == CARTESIAN || GEOMETRY == POLAR
          if(eps>=ZERO_F) {
            int som1;
            if(haveDomainDecomposition) {
              som1 = so - 1;
            } else {
              som1 = sbeg + modPositive(so-1-sbeg,n);
            }
            ek(k,s,i) = FargoFlux(scrhVs, BX1s, k, j, i, som1,
                                  n, sbeg, eps, haveDomainDecomposition);

          } else {
            ek(k,s,i) = FargoFlux(scrhVs, BX1s, k, j, i, so,
                                  n, sbeg, eps, haveDomainDecomposition);
          }
          if(m>0) {
            for(int ss = s-m ; ss < s ; ss++) {
              int sc;
              if(haveDomainDecomposition) {
                sc = ss + maxShift;
              } else {
                sc = sbeg + modPositive(ss-sbeg,n);
              }
              ek(k,s,i) += scrhVs(BX1s,k,sc,i);
            }
          } else {
            for(int ss = s ; ss < s-m ; ss++) {
              int sc;
              if(haveDomainDecomposition) {
                sc = ss + maxShift;
              } else {
                sc = sbeg + modPositive(ss-sbeg,n);
              }
              ek(k,s,i) -= scrhVs(BX1s,k,sc,i);
            }
          }
        #elif GEOMETRY == SPHERICAL
          if(eps>=ZERO_F) {
            int som1;
            if(haveDomainDecomposition) {
              som1 = so - 1;
            } else {
              som1 = sbeg + modPositive(so-1-sbeg,n);
            }
            ek(s,j,i) = FargoFlux(scrhVs, BX1s, k, j, i, som1,
                                  n, sbeg, eps, haveDomainDecomposition);

          } else {
            ek(s,j,i) = FargoFlux(scrhVs, BX1s, k, j, i, so,
                                  n, sbeg, eps, haveDomainDecomposition);
          }
          if(m>0) {
            for(int ss = s-m ; ss < s ; ss++) {
              int sc;
              if(haveDomainDecomposition) {
                sc = ss + maxShift;
              } else {
                sc = sbeg + modPositive(ss-sbeg,n);
              }
              ek(s,j,i) += scrhVs(BX1s,sc,j,i);
            }
          } else {
            for(int ss = s ; ss < s-m ; ss++) {
              int sc;
              if(haveDomainDecomposition) {
                sc = ss + maxShift;
              } else {
                sc = sbeg + modPositive(ss-sbeg,n);
              }
              ek(s,j,i) -= scrhVs(BX1s,sc,j,i);
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
        int so;
        if(haveDomainDecomposition) {
          so = s-m + maxShift;    // maxshift corresponds to the offset between
                                  // the indices in scrh and in Uc
        } else {
          so = sbeg + modPositive(s-m-sbeg,n);
        }

        // Compute EMF due to the shift via second order reconstruction
        #if GEOMETRY == CARTESIAN || GEOMETRY == POLAR
          if(eps>=ZERO_F) {
            int som1;
            if(haveDomainDecomposition) {
              som1 = so - 1;
            } else {
              som1 = sbeg + modPositive(so-1-sbeg,n);
            }
            ei(k,s,i) = FargoFlux(scrhVs, BX3s, k, j, i, som1,
                                  n, sbeg, eps, haveDomainDecomposition);
          } else {
            ei(k,s,i) = FargoFlux(scrhVs, BX3s, k, j, i, so,
                                  n, sbeg, eps, haveDomainDecomposition);
          }
          if(m>0) {
            for(int ss = s-m ; ss < s ; ss++) {
              int sc;
              if(haveDomainDecomposition) {
                sc = ss + maxShift;
              } else {
                sc = sbeg + modPositive(ss-sbeg,n);
              }
              ei(k,s,i) += scrhVs(BX3s,k,sc,i);
            }
          } else {
            for(int ss = s ; ss < s-m ; ss++) {
              int sc;
              if(haveDomainDecomposition) {
                sc = ss + maxShift;
              } else {
                sc = sbeg + modPositive(ss-sbeg,n);
              }
              ei(k,s,i) -= scrhVs(BX3s,k,sc,i);
            }
          }
        #elif GEOMETRY == SPHERICAL
        if(eps>=ZERO_F) {
          int som1;
          if(haveDomainDecomposition) {
            som1 = so - 1;
          } else {
            som1 = sbeg + modPositive(so-1-sbeg,n);
          }
          ei(s,j,i) = FargoFlux(scrhVs, BX2s, k, j, i, som1,
                                n, sbeg, eps, haveDomainDecomposition);
        } else {
          ei(s,j,i) = FargoFlux(scrhVs, BX2s, k, j, i, so,
                                n, sbeg, eps, haveDomainDecomposition);
        }
        if(m>0) {
          for(int ss = s-m ; ss < s ; ss++) {
            int sc;
            if(haveDomainDecomposition) {
              sc = ss + maxShift;
            } else {
              sc = sbeg + modPositive(ss-sbeg,n);
            }
            ei(s,j,i) += scrhVs(BX2s,sc,j,i);
          }
        } else {
          for(int ss = s ; ss < s-m ; ss++) {
            int sc;
            if(haveDomainDecomposition) {
              sc = ss + maxShift;
            } else {
              sc = sbeg + modPositive(ss-sbeg,n);
            }
            ei(s,j,i) -= scrhVs(BX2s,sc,j,i);
          }
        }

        #endif  // GEOMETRY

        ei(k,j,i) *= dphi;
      });
  #endif

    // Update field components according to the computed EMFS
    #ifndef EVOLVE_VECTOR_POTENTIAL
      IdefixArray4D<real> Vs = hydro->Vs;
      idefix_for("Fargo::EvolvMagField",
                data->beg[KDIR],data->end[KDIR]+KOFFSET,
                data->beg[JDIR],data->end[JDIR]+JOFFSET,
                data->beg[IDIR],data->end[IDIR]+IOFFSET,
        KOKKOS_LAMBDA (int k, int j, int i) {
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

    #else // EVOLVE_VECTOR_POTENTIAL
      // evolve field using vector potential
      IdefixArray4D<real> Ve = hydro->Ve;
      idefix_for("Fargo::EvolvMagField",
                data->beg[KDIR],data->end[KDIR]+KOFFSET,
                data->beg[JDIR],data->end[JDIR]+JOFFSET,
                data->beg[IDIR],data->end[IDIR]+IOFFSET,
        KOKKOS_LAMBDA (int k, int j, int i) {
          #if GEOMETRY == CARTESIAN
            #if DIMENSIONS == 3
              Ve(AX1e,k,j,i) += ex(k,j,i);
            #endif
            Ve(AX3e,k,j,i) += - ez(k,j,i);
          #elif GEOMETRY == POLAR
            #if DIMENSIONS == 3
              Ve(AX1e,k,j,i) += x1(i) * ex(k,j,i);
            #endif
            Ve(AX3e,k,j,i) +=  - x1m(i) * ez(k,j,i);
          #elif GEOMETRY == SPHERICAL
            #if DIMENSIONS == 3
              Ve(AX1e,k,j,i) += - x1(i) * sinx2m(j) * ex(k,j,i);
              Ve(AX2e,k,j,i) +=   x1m(i) * sinx2(j) * ey(k,j,i);
            #endif
          #endif
        });

    #endif // EVOLVE_VECTOR_POTENTIAL
  }
  #endif // GEOMETRY==CYLINDRICAL
  idfx::popRegion();
}

#endif // DATABLOCK_FARGO_HPP_
