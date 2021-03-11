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

void Fargo::Init(Input &input, Grid &grid, Hydro *hydro) {
  idfx::pushRegion("Fargo::Init");
  // Todo(lesurg): should work on the shearing box version of this.
  this->hydro = hydro;
  this->data = hydro->data;
  this->type=userdef;
  #if GEOMETRY == CARTESIAN || GEOMETRY == POLAR
    // Check there is no domain decomposition in the intended fargo direction
    if(data->mygrid->nproc[JDIR]>1) {
      IDEFIX_ERROR("Fargo is not yet compatible with MPI decomposition along the X2 direction");
    }
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
  idfx::cout << "Fargo: Enabled with user-defined velocity function" << std::endl;
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
  GetFargoVelocity(t);
  IdefixArray4D<real> Vc = hydro->Vc;
  IdefixArray2D<real> meanV = this->meanVelocity;
  idefix_for("FargoAddVelocity",
              0,data->np_tot[KDIR],
              0,data->np_tot[JDIR],
              0,data->np_tot[IDIR],
              KOKKOS_LAMBDA(int k, int j, int i) {
                #if GEOMETRY == CARTESIAN || GEOMETRY == POLAR
                  Vc(VX2,k,j,i) += meanV(k,i);
                #elif GEOMETRY == SPHERICAL
                  Vc(VX3,k,j,i) += meanV(j,i);
                #endif
              });

  idfx::popRegion();
}

void Fargo::SubstractVelocity(const real t) {
  idfx::pushRegion("Fargo::AddVelocity");

  GetFargoVelocity(t);
  IdefixArray4D<real> Vc = hydro->Vc;
  IdefixArray2D<real> meanV = this->meanVelocity;
  idefix_for("FargoAddVelocity",
              0,data->np_tot[KDIR],
              0,data->np_tot[JDIR],
              0,data->np_tot[IDIR],
              KOKKOS_LAMBDA(int k, int j, int i) {
                #if GEOMETRY == CARTESIAN || GEOMETRY == POLAR
                  Vc(VX2,k,j,i) -= meanV(k,i);
                #elif GEOMETRY == SPHERICAL
                  Vc(VX3,k,j,i) -= meanV(j,i);
                #endif
              });
  idfx::popRegion();
}


void Fargo::ShiftSolution(const real t, const real dt) {
  idfx::pushRegion("Fargo::ShiftSolution");
  // Refresh the fargo velocity function
  GetFargoVelocity(t);

  IdefixArray4D<real> Uc = hydro->Uc;
  IdefixArray4D<real> scrh = this->scratch;
  IdefixArray2D<real> meanV = this->meanVelocity;
  IdefixArray1D<real> x1 = data->x[IDIR];
  IdefixArray1D<real> x2 = data->x[JDIR];
  IdefixArray1D<real> dx2 = data->dx[JDIR];
  IdefixArray1D<real> dx3 = data->dx[KDIR];

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
                int s,sp1,sm1,sp2,sm2;
                #if GEOMETRY == CARTESIAN
                 w = meanV(k,i);
                 dphi = dx2(j);
                 s = j;
                #elif GEOMETRY == POLAR
                 w = meanV(k,i)/x1(i);
                 dphi = dx2(j);
                 s = j;
                #elif GEOMETRY == SPHERICAL
                 w = meanV(j,i)/(x1(i)*SIN(x2(j)));
                 dphi = dx3(k);
                 s = k;
                #endif
                // compute shifted indices, taking into account the fact that we're periodic
                sp1 = sbeg + (s+1-sbeg)%(send-sbeg);
                sp2 = sbeg + (s+2-sbeg)%(send-sbeg);
                sm1 = sbeg + (s-1-sbeg)%(send-sbeg);
                sm2 = sbeg + (s-2-sbeg)%(send-sbeg);
                // Compute the offset in phi, modulo the full domain size
                real dL = std::fmod(w*dt, Lphi);

                // Translate this into # of cells
                int m = static_cast<int> (std::floor(dL/dphi+HALF_F));

                // get the remainding shift
                real eps = dL/dphi - m;

                // target index after the shift
                int starget = sbeg+ (s+m-sbeg)%(send-sbeg);

                // Define Left and right fluxes
                // Fluxes are defined from slop-limited interpolation
                // Using Van-leer slope limiter (consistently with the main advection scheme)
                real Fl,Fr;
                real dqm, dqp, dq;

                #if GEOMETRY == CARTESIAN || GEOMETRY == POLAR
                  if(eps>=ZERO_F) {
                    // Compute Fl
                    dqm = Uc(n,k,sm1,i) - Uc(n,k,sm2,i);
                    dqp = Uc(n,k,s,i) - Uc(n,k,sm1,i);
                    dq = (dqp*dqm > ZERO_F ? TWO_F*dqp*dqm/(dqp + dqm) : ZERO_F);
                    Fl = Uc(n,k,sm1,i) + 0.5*dq*(1.0-eps);
                    //Compute Fr
                    dqm=dqp;
                    dqp = Uc(n,k,sp1,i) - Uc(n,k,s,i);
                    dq = (dqp*dqm > ZERO_F ? TWO_F*dqp*dqm/(dqp + dqm) : ZERO_F);
                    Fr = Uc(n,k,s,i) + 0.5*dq*(1.0-eps);
                  } else {
                    //Compute Fl
                    dqm = Uc(n,k,s,i) - Uc(n,k,sm1,i);
                    dqp = Uc(n,k,sp1,i) - Uc(n,k,s,i);
                    dq = (dqp*dqm > ZERO_F ? TWO_F*dqp*dqm/(dqp + dqm) : ZERO_F);
                    Fl = Uc(n,k,s,i) - 0.5*dq*(1.0+eps);
                    // Compute Fr
                    dqm=dqp;
                    dqp = Uc(n,k,sp2,i) - Uc(n,k,sp1,i);
                    dq = (dqp*dqm > ZERO_F ? TWO_F*dqp*dqm/(dqp + dqm) : ZERO_F);
                    Fr = Uc(n,k,sp1,i) - 0.5*dq*(1.0+eps);
                  }
                  scrh(n,k,starget,i) = Uc(n,k,s,i) - eps*(Fr - Fl);
                #elif GEOMETRY == SPHERICAL
                if(eps>=ZERO_F) {
                    // Compute Fl
                    dqm = Uc(n,sm1,j,i) - Uc(n,sm2,j,i);
                    dqp = Uc(n,s,j,i) - Uc(n,sm1,j,i);
                    dq = (dqp*dqm > ZERO_F ? TWO_F*dqp*dqm/(dqp + dqm) : ZERO_F);
                    Fl = Uc(n,sm1,j,i) + 0.5*dq*(1.0-eps);
                    //Compute Fr
                    dqm=dqp;
                    dqp = Uc(n,sp1,j,i) - Uc(n,s,j,i);
                    dq = (dqp*dqm > ZERO_F ? TWO_F*dqp*dqm/(dqp + dqm) : ZERO_F);
                    Fr = Uc(n,s,j,i) + 0.5*dq*(1.0-eps);
                  } else {
                    //Compute Fl
                    dqm = Uc(n,s,j,i) - Uc(n,sm1,j,i);
                    dqp = Uc(n,sp1,j,i) - Uc(n,s,j,i);
                    dq = (dqp*dqm > ZERO_F ? TWO_F*dqp*dqm/(dqp + dqm) : ZERO_F);
                    Fl = Uc(n,s,j,i) - 0.5*dq*(1.0+eps);
                    // Compute Fr
                    dqm=dqp;
                    dqp = Uc(n,sp2,j,i) - Uc(n,sp1,j,i);
                    dq = (dqp*dqm > ZERO_F ? TWO_F*dqp*dqm/(dqp + dqm) : ZERO_F);
                    Fr = Uc(n,sp1,j,i) - 0.5*dq*(1.0+eps);
                  }
                  scrh(n,starget,j,i) = Uc(n,s,j,i) - eps*(Fr - Fl);
                #endif
              });

  // move back scratch space into Uc
  Kokkos::deep_copy(Uc,scrh);

  idfx::popRegion();
}
