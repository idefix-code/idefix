// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2022 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include <string>
#include <vector>

#include "idefix.hpp"
#include "hydro.hpp"
#include "dataBlock.hpp"
#include "fargo.hpp"

// If no high order fargo, then choose the order according to reconstruction
#ifndef HIGH_ORDER_FARGO
  #if ORDER >= 3
  #define HIGH_ORDER_FARGO
  #endif
#endif

#ifdef HIGH_ORDER_FARGO
KOKKOS_FORCEINLINE_FUNCTION real PPMLim(real dvp, real dvm) {
  if(dvp*dvm >0.0) {
    real dqc = 0.5*(dvp+dvm);
    real d2q = 2.0*( fabs(dvp) < fabs(dvm) ? dvp : dvm);
    return( fabs(d2q) < fabs(dqc) ? d2q : dqc);
  }
  return(ZERO_F);
}

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
    real dqlm = PPMLim(dqm1,dqm2);
    real dql0 = PPMLim(dqp1,dqm1);
    real dqlp = PPMLim(dqp2,dqp1);

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


void Fargo::Init(Input &input, DataBlock *data) {
  idfx::pushRegion("Fargo::Init");
  this->data = data;
  this->hydro = &(data->hydro);

  // A bit of arithmetic to get the sizes of the working array
  this->nghost = data->nghost;
  this->beg = data->beg;
  this->end = data->end;

  if(input.CheckBlock("Fargo")) {
    std::string opType = input.Get<std::string>("Fargo","velocity",0);
    if(opType.compare("userdef")==0) {
       #if GEOMETRY != SPHERICAL && GEOMETRY != POLAR
        IDEFIX_ERROR("Fargo+userdef is only compatible with SPHERICAL and POLAR geometries");
      #endif
      this->type=userdef;
    } else if(opType.compare("shearingbox")==0) {
      this->type=shearingbox;
      #if GEOMETRY != CARTESIAN
        // Actually, this has never really been tested, so assumes it doesn't work...
        IDEFIX_ERROR("Fargo+shearingbox is only compatible with cartesian geometry");
      #endif
    } else {
      IDEFIX_ERROR("Unknown fargo velocity in the input file. "
      "Only userdef and shearingbox are allowed");
    }
    this->maxShift = input.GetOrSet<int>("Fargo", "maxShift",0, 10);
  } else {
    // DEPRECATED: initialisation from the [Hydro] block
    if(input.CheckEntry("Hydro","fargo")>=0) {
      std::string opType = input.Get<std::string>("Hydro","fargo",0);
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
  }

  #if GEOMETRY == CARTESIAN || GEOMETRY == POLAR
    // Check if there is a domain decomposition in the intended fargo direction
    if(data->mygrid->nproc[JDIR]>1) {
      haveDomainDecomposition = true;
      this->nghost[JDIR] += this->maxShift;
      this->beg[JDIR] += this->maxShift;
      this->end[JDIR] += this->maxShift;
      if(data->np_int[JDIR] < this->maxShift + data->nghost[JDIR]) {
        IDEFIX_ERROR("Subdomain size < Fargo:maxShift. "
                     "Try reducting the number of processes along X2");
      }
    }
    if(this->type==userdef)
      this->meanVelocity = IdefixArray2D<real>("FargoVelocity",data->np_tot[KDIR],
                                                             data->np_tot[IDIR]);
  #elif GEOMETRY == SPHERICAL
    // Check if there is a domain decomposition in the intended fargo direction
    if(data->mygrid->nproc[KDIR]>1) {
      haveDomainDecomposition = true;
      this->nghost[KDIR] += this->maxShift;
      this->beg[KDIR] += this->maxShift;
      this->end[KDIR] += this->maxShift;
      if(data->np_int[KDIR] < this->maxShift + data->nghost[KDIR]) {
        IDEFIX_ERROR("Subdomain size < Fargo:maxShift. "
                     "Try reducting the number of processes along X3");
      }
    }
    this->meanVelocity = IdefixArray2D<real>("FargoVelocity",data->np_tot[JDIR],
                                                             data->np_tot[IDIR]);
  #else
    IDEFIX_ERROR("Fargo is not compatible with the GEOMETRY you intend to use");
  #endif


  // Initialise our scratch space
  this->scrhUc = IdefixArray4D<real>("FargoVcScratchSpace",NVAR
                                      ,end[KDIR]-beg[KDIR] + 2*nghost[KDIR]
                                      ,end[JDIR]-beg[JDIR] + 2*nghost[JDIR]
                                      ,end[IDIR]-beg[IDIR] + 2*nghost[IDIR]);

  #if MHD == YES
    if(haveDomainDecomposition) {
      this->scrhVs = IdefixArray4D<real>("FargoVsScratchSpace",DIMENSIONS
                                          ,end[KDIR]-beg[KDIR] + 2*nghost[KDIR]+KOFFSET
                                          ,end[JDIR]-beg[JDIR] + 2*nghost[JDIR]+JOFFSET
                                          ,end[IDIR]-beg[IDIR] + 2*nghost[IDIR]+IOFFSET);


    } else {
      // A separate allocation for scrhVs is only needed with domain decomposition, otherwise,
      // we just make a reference to scrhVs
      this->scrhVs = hydro->Vs;
    }
  #endif
  #ifdef WITH_MPI
    if(haveDomainDecomposition) {
      std::vector<int> vars;
      for(int i=0 ; i < NVAR ; i++) {
        vars.push_back(i);
      }
      #if MHD == YES
        this->mpi.Init(data->mygrid, vars, this->nghost.data(), data->np_int.data(), true);
      #else
        this->mpi.Init(data->mygrid, vars, this->nghost.data(), data->np_int.data());
      #endif
    }
  #endif


  idfx::popRegion();
}

void Fargo::ShowConfig() {
  idfx::pushRegion("Fargo::ShowConfig");
  if(type==userdef) {
    idfx::cout << "Fargo: ENABLED with user-defined velocity function." << std::endl;
    if(!fargoVelocityFunc) {
      IDEFIX_ERROR("No Fargo velocity function has been enabled.");
    }
  } else if(type==shearingbox) {
    idfx::cout << "Fargo: ENABLED with shearing-box velocity function." << std::endl;
  } else {
    IDEFIX_ERROR("Something went wrong during Fargo initialisation");
  }
  #ifdef HIGH_ORDER_FARGO
    idfx::cout << "Fargo: using high order PPM advection scheme." << std::endl;
  #else
    idfx::cout << "Fargo: using standard PLM advection scheme." << std::endl;
  #endif
  if(haveDomainDecomposition) {
    idfx::cout << "Fargo: using domain decomposition along the azimuthal direction"
               << " with maxShift=" << this->maxShift << std::endl;
  }
  idfx::popRegion();
}

void Fargo::EnrollVelocity(FargoVelocityFunc myFunc) {
  if(this->type!=userdef) {
    IDEFIX_WARNING("Fargo velocity function enrollment requires Hydro/Fargo "
                 "to be set to userdef in .ini file");
  }
  this->fargoVelocityFunc = myFunc;
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
    if(this->haveDomainDecomposition) {
      CheckMaxDisplacement();
    }
  }
  idfx::popRegion();
}

// This function checks that the velocity provided does not lead to an azimuthal displacement
// larger than the one allowed
void Fargo::CheckMaxDisplacement() {
    IdefixArray2D<real> meanV = this->meanVelocity;
    int ibeg, iend, jbeg, jend;
    IdefixArray1D<real> xi;
    IdefixArray1D<real> xj;
    IdefixArray1D<real> dxk;
    [[maybe_unused]] FargoType fargoType = type;
    [[maybe_unused]] real sbS = hydro->sbS;
    real invDt = 0;

    // Get domain size
    #if GEOMETRY == CARTESIAN || GEOMETRY == POLAR
      ibeg = data->beg[IDIR];
      iend = data->end[IDIR];
      jbeg = data->beg[KDIR];
      jend = data->end[KDIR];
      xi = data->x[IDIR];
      xj = data->x[KDIR];
      dxk = data->dx[JDIR];
    #elif GEOMETRY == SPHERICAL
      ibeg = data->beg[IDIR];
      iend = data->end[IDIR];
      jbeg = data->beg[JDIR];
      jend = data->end[JDIR];
      xi = data->x[IDIR];
      xj = data->x[JDIR];
      dxk = data->dx[KDIR];
    #endif


    idefix_reduce("CheckMaxDt", jbeg, jend, ibeg, iend,
      KOKKOS_LAMBDA(int j, int i, real &invDtLoc) {
        real w,dphi;
        #if GEOMETRY == CARTESIAN
          if(fargoType==userdef) {
            w = meanV(j,i);
          } else if(fargoType==shearingbox) {
            w = sbS*xi(i);
          }
        #elif GEOMETRY == POLAR
          w = meanV(j,i)/xi(i);
        #elif GEOMETRY == SPHERICAL
          w = meanV(j,i)/(xi(i)*sin(xj(j)));
        #endif
        dphi = dxk(0);  // dphi is supposedly constant when using fargo.
        invDtLoc = FMAX(invDtLoc, FABS(w/dphi));
      },
      Kokkos::Max<real>(invDt));
  #ifdef WITH_MPI
    if(idfx::psize>1) {
          MPI_SAFE_CALL(MPI_Allreduce(MPI_IN_PLACE, &invDt, 1, realMPI, MPI_MAX, MPI_COMM_WORLD));
        }
  #endif
  this->dtMax = this->maxShift / invDt;
}

void Fargo::AddVelocity(const real t) {
  idfx::pushRegion("Fargo::AddVelocity");
  if(type==userdef) {
    GetFargoVelocity(t);
  }
  IdefixArray1D<real> x1 = data->x[IDIR];
  IdefixArray4D<real> Vc = hydro->Vc;
  IdefixArray2D<real> meanV = this->meanVelocity;
  [[maybe_unused]] FargoType fargoType = type;
  [[maybe_unused]] real sbS = hydro->sbS;

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
  [[maybe_unused]] IdefixArray2D<real> meanV = this->meanVelocity;
  [[maybe_unused]] FargoType fargoType = type;
  [[maybe_unused]] real sbS = hydro->sbS;

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

void Fargo::StoreToScratch() {
  IdefixArray4D<real> Uc = hydro->Uc;
  IdefixArray4D<real> scrhUc = this->scrhUc;
  bool haveDomainDecomposition = this->haveDomainDecomposition;
  int maxShift = this->maxShift;

  idefix_for("Fargo:StoreUc",
            0,NVAR,
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

  #if MHD == YES
    #ifdef EVOLVE_VECTOR_POTENTIAL
      // Update Vs to its latest
      hydro->emf.ComputeMagFieldFromA(hydro->Ve,hydro->Vs);
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
  #endif
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

void Fargo::ShiftSolution(const real t, const real dt) {
  idfx::pushRegion("Fargo::ShiftSolution");

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
  StoreToScratch();

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
                int so;
                if(haveDomainDecomposition) {
                  so = s-m + maxShift;    // maxshift corresponds to the offset between
                                          // the indices in scrh and in Uc
                } else {
                  so = sbeg + modPositive(s-m-sbeg, ds);
                }

                // Define Left and right fluxes
                // Fluxes are defined from slop-limited interpolation
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

#if MHD == YES
  IdefixArray4D<real> scrhVs = this->scrhVs;
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
          ek(k,s,i) = FargoFlux(scrhVs, BX1s, k, j, i, som1, n, sbeg, eps, haveDomainDecomposition);

        } else {
          ek(k,s,i) = FargoFlux(scrhVs, BX1s, k, j, i, so, n, sbeg, eps, haveDomainDecomposition);
        }
        if(m>0) {
          for(int ss = s-m ; ss < s ; ss++) {
            int sc;
            if(haveDomainDecomposition) {
              sc = ss;
            } else {
              sc = sbeg + modPositive(ss-sbeg,n);
            }
            ek(k,s,i) += scrhVs(BX1s,k,sc,i);
          }
        } else {
          for(int ss = s ; ss < s-m ; ss++) {
            int sc;
            if(haveDomainDecomposition) {
              sc = ss;
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
          ek(s,j,i) = FargoFlux(scrhVs, BX1s, k, j, i, som1, n, sbeg, eps, haveDomainDecomposition);

        } else {
          ek(s,j,i) = FargoFlux(scrhVs, BX1s, k, j, i, so, n, sbeg, eps, haveDomainDecomposition);
        }
        if(m>0) {
          for(int ss = s-m ; ss < s ; ss++) {
            int sc;
            if(haveDomainDecomposition) {
              sc = ss;
            } else {
              sc = sbeg + modPositive(ss-sbeg,n);
            }
            ek(s,j,i) += scrhVs(BX1s,sc,j,i);
          }
        } else {
          for(int ss = s ; ss < s-m ; ss++) {
            int sc;
            if(haveDomainDecomposition) {
              sc = ss;
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
          ei(k,s,i) = FargoFlux(scrhVs, BX3s, k, j, i, som1, n, sbeg, eps, haveDomainDecomposition);
        } else {
          ei(k,s,i) = FargoFlux(scrhVs, BX3s, k, j, i, so, n, sbeg, eps, haveDomainDecomposition);
        }
        if(m>0) {
          for(int ss = s-m ; ss < s ; ss++) {
            int sc;
            if(haveDomainDecomposition) {
              sc = ss;
            } else {
              sc = sbeg + modPositive(ss-sbeg,n);
            }
            ei(k,s,i) += scrhVs(BX3s,k,sc,i);
          }
        } else {
          for(int ss = s ; ss < s-m ; ss++) {
            int sc;
            if(haveDomainDecomposition) {
              sc = ss;
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
        ei(s,j,i) = FargoFlux(scrhVs, BX2s, k, j, i, som1, n, sbeg, eps, haveDomainDecomposition);
      } else {
        ei(s,j,i) = FargoFlux(scrhVs, BX2s, k, j, i, so, n, sbeg, eps, haveDomainDecomposition);
      }
      if(m>0) {
        for(int ss = s-m ; ss < s ; ss++) {
          int sc;
          if(haveDomainDecomposition) {
            sc = ss;
          } else {
            sc = sbeg + modPositive(ss-sbeg,n);
          }
          ei(s,j,i) += scrhVs(BX2s,sc,j,i);
        }
      } else {
        for(int ss = s ; ss < s-m ; ss++) {
          int sc;
          if(haveDomainDecomposition) {
            sc = ss;
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


#endif // MHD
#endif // GEOMETRY==CYLINDRICAL

  idfx::popRegion();
}
