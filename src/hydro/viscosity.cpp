// ***********************************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************************

// This source code is largely inspired from the viscous_flux of Pluto4.2 
// ((c) P. Tzeferacos & A. Mignone)

#include <string>

#include "../idefix.hpp"
#include "viscosity.hpp"


#define D_DX_I(q,n)  (q(n,k,j,i) - q(n,k,j,i - 1))
#define D_DY_J(q,n)  (q(n,k,j,i) - q(n,k,j - 1,i))
#define D_DZ_K(q,n)  (q(n,k,j,i) - q(n,k - 1,j,i))

#define D_DY_I(q,n)  (  0.25*(q(n,k,j + 1,i) + q(n,k,j + 1,i - 1)) \
                    - 0.25*(q(n,k,j - 1,i) + q(n,k,j - 1,i - 1)))

#define D_DZ_I(q,n)  (  0.25*(q(n,k + 1,j,i) + q(n,k + 1,j,i - 1))  \
                    - 0.25*(q(n,k - 1,j,i) + q(n,k - 1,j,i - 1)))

#define D_DX_J(q,n)  (  0.25*(q(n,k,j,i + 1) + q(n,k,j - 1,i + 1)) \
                    - 0.25*(q(n,k,j,i - 1) + q(n,k,j - 1,i - 1)))

#define D_DZ_J(q,n)  (  0.25*(q(n,k + 1,j,i) + q(n,k + 1,j - 1,i)) \
                    - 0.25*(q(n,k - 1,j,i) + q(n,k - 1,j - 1,i)))

#define D_DX_K(q,n)  (  0.25*(q(n,k,j,i + 1) + q(n,k - 1,j,i + 1)) \
                    - 0.25*(q(n,k,j,i - 1) + q(n,k - 1,j,i - 1)))

#define D_DY_K(q,n)  (  0.25*(q(n,k,j + 1,i) + q(n,k - 1,j + 1,i)) \
                    - 0.25*(q(n,k,j - 1,i) + q(n,k - 1,j - 1,i)))                


Viscosity::Viscosity() {
  // Default constructor
}

void Viscosity::Init(Input &input, Grid &grid, Hydro *hydroin) {
  idfx::pushRegion("Viscosity::Init");
  // Save the parent hydro object
  this->hydro = hydroin;
  idfx::cout << "Hydro:" << hydroin << std::endl;
  if(input.CheckEntry("Hydro","Viscosity")>=0) {
    if(input.GetString("Hydro","Viscosity",0).compare("constant") == 0) {
        idfx::cout << "Viscosity: Enabling constant viscosity function." << std::endl;
        this->eta1 = input.GetReal("Hydro","Viscosity",1);
        // second viscosity?
        if(input.CheckEntry("Hydro","Viscosity")>2) {
          this->eta2 = input.GetReal("Hydro","Viscosity",2);
        } else {
          this->eta2 = 0.0;
          idfx::cout << "Viscosity: Second viscosity not provided. Assuming it is 0." << std::endl;
        }
        this->haveViscosity = Constant;
      } else if(input.GetString("Hydro","Viscosity",0).compare("userdef") == 0) {
        idfx::cout << "Viscosity: Enabling user-defined viscosity function."
                   << std::endl;
        this->haveViscosity = UserDefFunction;
      } else {
        IDEFIX_ERROR("Unknown viscosity definition in idefix.ini. "
                     "Can only be constant or userdef.");
      }
  } else {
    IDEFIX_ERROR("I cannot create a Viscosity object without Viscosity defined"
                   "in the .ini file");
  }

  // Allocate and fill arrays when needed
  #if GEOMETRY != CARTESIAN
    one_dmu = IDefixArray1D<real>("Viscosity_1dmu", hydro->data->np_tot[JDIR]);
    IdefixArray1D<real> th = hydro->data->x[JDIR];
    idefix_for("ViscousInitGeometry",1,hydro->data->np_tot[JDIR],
      KOKKOS_LAMBDA(int j) {
        real scrch = 1.0-cos(th(j))-(1.0-cos(th(j-1))) * (th(j-1) > 0.0 ? 1.0:-1.0);
        one_dmu(j) = 1.0/scrch;

      });
  #endif
  viscSrc = IdefixArray4D<real>("Viscosity_source", COMPONENTS, hydro->data->np_tot[KDIR], 
                                                                hydro->data->np_tot[JDIR], 
                                                                hydro->data->np_tot[IDIR]);
  idfx::popRegion();
}


void Viscosity::EnrollViscousDiffusivity(DiffusivityFunc myFunc) {
  if(this->haveViscosity < UserDefFunction) {
    IDEFIX_ERROR("Viscous diffusivity enrollment requires Hydro/Viscosity "
                 "to be set to userdef in .ini file");
  }
  this->viscousDiffusivityFunc = myFunc;
  idfx::cout << "Viscosity: User-defined viscous diffusion has been enrolled." << std::endl;
}

// This function computes the viscous flux and stores it in hydro->fluxRiemann (this avoids an extra array)
// Associated source terms, present in non-cartesian geometry are also computed
// and stored in this->viscSrc for later use (in calcRhs).
void Viscosity::AddViscousFlux(int dir, const real t) {
  idfx::pushRegion("Viscosity::AddViscousFlux");
  IdefixArray4D<real> Vc = this->hydro->Vc;
  IdefixArray4D<real> Flux = this->hydro->FluxRiemann;
  IdefixArray3D<real> dMax = this->hydro->dMax;
  IdefixArray1D<real> x1 = this->hydro->data->x[IDIR];
  IdefixArray1D<real> x2 = this->hydro->data->x[JDIR];
  IdefixArray1D<real> x3 = this->hydro->data->x[KDIR];
  IdefixArray1D<real> x1l = this->hydro->data->xl[IDIR];
  IdefixArray1D<real> x2l = this->hydro->data->xl[JDIR];
  IdefixArray1D<real> x3l = this->hydro->data->xl[KDIR];
  IdefixArray1D<real> dx1 = this->hydro->data->dx[IDIR];
  IdefixArray1D<real> dx2 = this->hydro->data->dx[JDIR];
  IdefixArray1D<real> dx3 = this->hydro->data->dx[KDIR];


  int ibeg, iend, jbeg, jend, kbeg, kend;
  ibeg = this->hydro->data->beg[IDIR];
  iend = this->hydro->data->end[IDIR];
  jbeg = this->hydro->data->beg[JDIR];
  jend = this->hydro->data->end[JDIR];
  kbeg = this->hydro->data->beg[KDIR];
  kend = this->hydro->data->end[KDIR];

  // Determine the offset along which we do the extrapolation
  if(dir==IDIR) iend++;
  if(dir==JDIR) jend++;
  if(dir==KDIR) kend++;

  real eta1 = this->eta1;
  real eta2 = this->eta2; 

  idefix_for("ViscousFlux",kbeg, kend, jbeg, jend, ibeg, iend,
    KOKKOS_LAMBDA (int k, int j, int i) {
      real tau_xx, tau_xy, tau_xz;
      real tau_yy, tau_yz;
      real tau_zz;

      real dVxi, dVxj, dVxk;
      real dVyi, dVyj, dVyk;
      real dVzi, dVzj, dVzk;

      real divV;

      tau_xx = tau_xy = tau_xz = ZERO_F;
      tau_yy = tau_yz = ZERO_F;
      tau_zz = ZERO_F;

      dVxi = dVxj = dVxk = ZERO_F;
      dVyi = dVyj = dVyk = ZERO_F;
      dVzi = dVzj = dVzk = ZERO_F;

      if(dir == IDIR) {
        EXPAND(  dVxi = D_DX_I(Vc,VX1)/dx1(i); , 
                 dVyi = D_DX_I(Vc,VX2)/dx1(i); , 
                 dVzi = D_DX_I(Vc,VX3)/dx1(i); )

        #if DIMENSIONS >= 2
          EXPAND(  dVxj = D_DY_I(Vc,VX1)/dx2(j); , 
                   dVyj = D_DY_I(Vc,VX2)/dx2(j); , 
                   dVzj = D_DY_I(Vc,VX3)/dx2(j); )
          #if DIMENSIONS == 3
            dVxk = D_DZ_I(Vc,VX1)/dx3(k); 
            dVyk = D_DZ_I(Vc,VX2)/dx3(k); 
            dVzk = D_DZ_I(Vc,VX3)/dx3(k);
          #endif
        #endif

        #if GEOMETRY == CARTESIAN
          divV = D_EXPAND(dVxi, + dVyj, + dVzk);

          tau_xx = 2.0*eta1*dVxi + (eta2 - (2.0/3.0)*eta1)*divV;
          tau_xy = eta1*(dVxj + dVyi);
          tau_xz = eta1*(dVxk + dVzi);
        #endif

        // Update flux with the stress tensor
        EXPAND( Flux(MX1, k, j, i) -= tau_xx; ,
                Flux(MX2, k, j, i) -= tau_xy; ,
                Flux(MX3, k, j, i) -= tau_xz; )
        
        #if HAVE_ENERGY
          Flux(ENG, k, j, i) -= EXPAND(   0.5*(Vc(VX1,k,j,i) + Vc(VX1,k,j,i-1))*tau_xx  ,  
                                        + 0.5*(Vc(VX2,k,j,i) + Vc(VX2,k,j,i-1))*tau_xy  ,
                                        + 0.5*(Vc(VX3,k,j,i) + Vc(VX3,k,j,i-1))*tau_xz) ;
        #endif
        

        dMax(k,j,i) += (FMAX(eta1,eta2))/(0.5*(Vc(RHO,k,j,i)+Vc(RHO,k,j,i-1)));
      }
      if(dir == JDIR) {
        EXPAND(  dVxi = D_DX_J(Vc,VX1)/dx1(i);, 
                 dVyi = D_DX_J(Vc,VX2)/dx1(i);, 
                 dVzi = D_DX_J(Vc,VX3)/dx1(i); )

        EXPAND(  dVxj = D_DY_J(Vc,VX1)/dx2(j);, 
                 dVyj = D_DY_J(Vc,VX2)/dx2(j);, 
                 dVzj = D_DY_J(Vc,VX3)/dx2(j); )

        #if DIMENSIONS == 3
          dVxk = D_DZ_J(Vc,VX1)/dx3(k); 
          dVyk = D_DZ_J(Vc,VX2)/dx3(k); 
          dVzk = D_DZ_J(Vc,VX3)/dx3(k);
        #endif
    
        #if GEOMETRY == CARTESIAN
          divV = D_EXPAND(dVxi, + dVyj, + dVzk);

          tau_xy = eta1*(dVxj+dVyi);
          tau_yy = 2.0*eta1*dVyj + (eta2 - (2.0/3.0)*eta1)*divV;
          tau_yz = eta1*(dVzj + dVyk);
        #endif

        // Update flux with the stress tensor
        EXPAND( Flux(MX1, k, j, i) -= tau_xy; ,
                Flux(MX2, k, j, i) -= tau_yy; ,
                Flux(MX3, k, j, i) -= tau_yz; )
        
        #if HAVE_ENERGY
          Flux(ENG, k, j, i) -= EXPAND(   0.5*(Vc(VX1,k,j,i) + Vc(VX1,k,j-1,i))*tau_xy  ,  
                                        + 0.5*(Vc(VX2,k,j,i) + Vc(VX2,k,j-1,i))*tau_yy  ,
                                        + 0.5*(Vc(VX3,k,j,i) + Vc(VX3,k,j-1,i))*tau_yz) ;
        #endif

        dMax(k,j,i) += (FMAX(eta1,eta2))/(0.5*(Vc(RHO,k,j,i)+Vc(RHO,k,j-1,i)));
      }
      if(dir == KDIR) {
        dVxi = D_DX_K(Vc,VX1)/dx1(i); 
        dVyi = D_DX_K(Vc,VX2)/dx1(i); 
        dVzi = D_DX_K(Vc,VX3)/dx1(i);
        dVxj = D_DY_K(Vc,VX1)/dx2(j);
        dVyj = D_DY_K(Vc,VX2)/dx2(j);
        dVzj = D_DY_K(Vc,VX3)/dx2(j); 
        dVxk = D_DZ_K(Vc,VX1)/dx3(k);
        dVyk = D_DZ_K(Vc,VX2)/dx3(k); 
        dVzk = D_DZ_K(Vc,VX3)/dx3(k);

        #if GEOMETRY == CARTESIAN
          divV = dVxi + dVyj + dVzk;

          tau_xz = eta1*(dVxk+dVzi);
          tau_yz = eta1*(dVyk+dVzj);
          tau_zz = 2.0*eta1*dVzk + (eta2 - (2.0/3.0)*eta1)*divV;
        #endif

        // Update flux with the stress tensor
        EXPAND( Flux(MX1, k, j, i) -= tau_xz; ,
                Flux(MX2, k, j, i) -= tau_yz; ,
                Flux(MX3, k, j, i) -= tau_zz; )
        
        #if HAVE_ENERGY
          Flux(ENG, k, j, i) -= EXPAND(   0.5*(Vc(VX1,k,j,i) + Vc(VX1,k-1,j,i))*tau_xz  ,  
                                        + 0.5*(Vc(VX2,k,j,i) + Vc(VX2,k-1,j,i))*tau_yz  ,
                                        + 0.5*(Vc(VX3,k,j,i) + Vc(VX3,k-1,j,i))*tau_zz) ;
        #endif

        dMax(k,j,i) += (FMAX(eta1,eta2))/(0.5*(Vc(RHO,k,j,i)+Vc(RHO,k-1,j,i)));
      }
    });


  idfx::popRegion();
}