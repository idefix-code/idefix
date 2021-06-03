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
  IdefixArray1D<real> sinx2 = data->sinx2;
  IdefixArray1D<real> sinx2m = data->sinx2m;

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
                 w = meanV(j,i)/(x1(i)*sinx2(j));
                 dphi = dx3(k);
                 s = k;
                #endif
                // compute shifted indices, taking into account the fact that we're periodic
                int ds = send-sbeg;
                sp1 = sbeg + ((s+1-sbeg)%ds+ds)%ds;
                sp2 = sbeg + ((s+2-sbeg)%ds+ds)%ds;
                sm1 = sbeg + ((s-1-sbeg)%ds+ds)%ds;
                sm2 = sbeg + ((s-2-sbeg)%ds+ds)%ds;
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
        w = 0.5*(meanV(k,i-1)+meanV(k,i));
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

      // compute shifted indices, taking into account the fact that we're periodic
      int sop1 = sbeg + ((so+1-sbeg)%n+n)%n;
      int som1 = sbeg + ((so-1-sbeg)%n+n)%n;
      int som2 = sbeg + ((so-2-sbeg)%n+n)%n;

      // Compute EMF due to the shift via second order reconstruction
      real dqm, dqp, dq;

      #if GEOMETRY == CARTESIAN || GEOMETRY == POLAR
        if(eps>=ZERO_F) {
          // Compute extrapolated ek
          dqm = Vs(BX1s,k,som1,i) - Vs(BX1s,k,som2,i);
          dqp = Vs(BX1s,k,so,i) - Vs(BX1s,k,som1,i);
          dq = (dqp*dqm > ZERO_F ? TWO_F*dqp*dqm/(dqp + dqm) : ZERO_F);
          ek(k,s,i) = eps*(Vs(BX1s,k,som1,i) + 0.5*dq*(1.0-eps));
        } else {
          dqm = Vs(BX1s,k,so,i) - Vs(BX1s,k,som1,i);
          dqp = Vs(BX1s,k,sop1,i) - Vs(BX1s,k,so,i);
          dq = (dqp*dqm > ZERO_F ? TWO_F*dqp*dqm/(dqp + dqm) : ZERO_F);
          ek(k,s,i) = eps*(Vs(BX1s,k,so,i) - 0.5*dq*(1.0+eps));
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
        // Compute extrapolated ek
        dqm = Vs(BX1s,som1,j,i) - Vs(BX1s,som2,j,i);
        dqp = Vs(BX1s,so,j,i) - Vs(BX1s,som1,j,i);
        dq = (dqp*dqm > ZERO_F ? TWO_F*dqp*dqm/(dqp + dqm) : ZERO_F);
        ek(s,j,i) = eps*(Vs(BX1s,som1,j,i) + 0.5*dq*(1.0-eps));
      } else {
        dqm = Vs(BX1s,so,j,i) - Vs(BX1s,som1,j,i);
        dqp = Vs(BX1s,sop1,j,i) - Vs(BX1s,so,j,i);
        dq = (dqp*dqm > ZERO_F ? TWO_F*dqp*dqm/(dqp + dqm) : ZERO_F);
        ek(s,j,i) = eps*(Vs(BX1s,so,j,i) - 0.5*dq*(1.0+eps));
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
        w = 0.5*(meanV(k,i)+meanV(k-1,i));
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
          // Compute extrapolated ek
          dqm = Vs(BX3s,k,som1,i) - Vs(BX3s,k,som2,i);
          dqp = Vs(BX3s,k,so,i) - Vs(BX3s,k,som1,i);
          dq = (dqp*dqm > ZERO_F ? TWO_F*dqp*dqm/(dqp + dqm) : ZERO_F);
          ei(k,s,i) = eps*(Vs(BX3s,k,som1,i) + 0.5*dq*(1.0-eps));
        } else {
          dqm = Vs(BX3s,k,so,i) - Vs(BX3s,k,som1,i);
          dqp = Vs(BX3s,k,sop1,i) - Vs(BX3s,k,so,i);
          dq = (dqp*dqm > ZERO_F ? TWO_F*dqp*dqm/(dqp + dqm) : ZERO_F);
          ei(k,s,i) = eps*(Vs(BX3s,k,so,i) - 0.5*dq*(1.0+eps));
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
        // Compute extrapolated ek
        dqm = Vs(BX2s,som1,j,i) - Vs(BX2s,som2,j,i);
        dqp = Vs(BX2s,so,j,i) - Vs(BX2s,som1,j,i);
        dq = (dqp*dqm > ZERO_F ? TWO_F*dqp*dqm/(dqp + dqm) : ZERO_F);
        ei(s,j,i) = eps*(Vs(BX2s,som1,j,i) + 0.5*dq*(1.0-eps));
      } else {
        dqm = Vs(BX2s,so,j,i) - Vs(BX2s,som1,j,i);
        dqp = Vs(BX2s,sop1,j,i) - Vs(BX2s,so,j,i);
        dq = (dqp*dqm > ZERO_F ? TWO_F*dqp*dqm/(dqp + dqm) : ZERO_F);
        ei(s,j,i) = eps*(Vs(BX2s,so,j,i) - 0.5*dq*(1.0+eps));
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
  this->hydro->ReconstructVcField(Uc);

#endif // MHD

  idfx::popRegion();
}
