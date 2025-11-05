// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include <string>

#include "forcing.hpp"
#include "dataBlock.hpp"
#include "dataBlockHost.hpp"
#include "input.hpp"

KOKKOS_INLINE_FUNCTION real K_aff_11(real x, real x0, real x1) {
  return 2./(x1-x0)*x - (x1+x0)/(x1-x0);
}

KOKKOS_INLINE_FUNCTION real K_aff_02pi(real x, real x0, real x1) {
  return 2.*M_PI*(x - x0)/(x1 - x0);
}

KOKKOS_INLINE_FUNCTION real K_cheby_fst(int n, real x) {
  if (x >= -1. and x <= 1.) {
    return cos(n*acos(x));
  } else {
    return ZERO_F;
  }
}

KOKKOS_INLINE_FUNCTION real K_cheby_fst_both_hom_dir(int n, real x) {
  if (x >= -1. and x <= 1.) {
    return cos(n*acos(x)) - 1.*((n+1)%2) - x*(n%2);
  } else {
    return ZERO_F;
  }
}

KOKKOS_INLINE_FUNCTION real K_cheby_fst_both_hom_neu(int n, real x) {
  if (x >= -1. and x <= 1.) {
    return cos(n*acos(x)) - pow(n/(n+2),2.)*cos((n+2)*acos(x));
//    return cos(n*acos(x)) - pow(n/(n+2),2.)*cos((n+2)*acos(x))*((n+1)%2) - pow(n/(n+2),2.)*cos((n+2)*acos(x))*(n%2);
  } else {
    return ZERO_F;
  }
}

KOKKOS_INLINE_FUNCTION real K_fourier(int n, real x) {
  if (x >= 0. and x <= 2.*M_PI) {
    return cos(n*x) + sin(n*x);
  } else {
    return ZERO_F;
  }
}

KOKKOS_INLINE_FUNCTION real K_sin(int n, real x) {
  if (x >= 0. and x <= 2.*M_PI) {
    return sin(n*x);
  } else {
    return ZERO_F;
  }
}

KOKKOS_INLINE_FUNCTION real K_cos(int n, real x) {
  if (x >= 0. and x <= 2.*M_PI) {
    return cos(n*x);
  } else {
    return ZERO_F;
  }
}

Forcing::Forcing(Input &input, DataBlock *datain) {
  idfx::pushRegion("Forcing::Forcing");
  this->data = datain;
  this->seed = input.GetOrSet<int>("Forcing","seed",0,0);

  this->write = input.GetOrSet<int>("Forcing","write",0, 0);
//  std::string folder = input.GetOrSet<std::string>("Forcing","filename",0,"testOU");
  this->folder = input.GetOrSet<std::string>("Output","folder",0,"output");

  this->stillHaveForcing = true;
  this->stopTime = input.GetOrSet<real>("Forcing","stoptime",0,std::numeric_limits<double>::infinity());

  this->nForcingModes = 0;

  this->kmin = -1.;
  this->kmax = -1.;
  this->ellmin = -1;
  this->ellmax = -1;
  this->mmin = -1;
  this->mmax = -1;
  this->haveSolenoidalForcing = false;

  this->xbeg = data->mygrid->xbeg[IDIR];
  this->xend = data->mygrid->xend[IDIR];
  this->ybeg = data->mygrid->xbeg[JDIR];
  this->yend = data->mygrid->xend[JDIR];
  this->zbeg = data->mygrid->xbeg[KDIR];
  this->zend = data->mygrid->xend[KDIR];

  this->kx0 = 2.*M_PI/(xend - xbeg);
  this->ky0 = 2.*M_PI/(yend - ybeg);
  this->kz0 = 2.*M_PI/(zend - zbeg);

  if (input.CheckEntry("Forcing", "iso3D")>=0) {
    this->forcingType = iso3D;
    #if GEOMETRY != CARTESIAN
      IDEFIX_ERROR("You cannot have 3D isotropic forcing in another geometry than Cartesian");
    #endif //GEOMETRY != CARTESIAN
    #if COMPONENTS < 3 or DIMENSIONS < 3
std::cout << COMPONENTS << DIMENSIONS << std::endl;
      IDEFIX_ERROR("You cannot have 3D isotropic forcing with less than 3 components and dimensions.");
    #endif //COMPONENTS < 3 or DIMENSIONS < 3
// if () IDEFIX_ERROR("You should not have 3D isotropic forcing with other boundary conditions than fully periodic.");
    #if GEOMETRY == SPHERICAL
      IDEFIX_ERROR("Comparing radial with angular wavenumbers is somewhat shady.");
    #endif //GEOMETRY == SPHERICAL
    this->kmin = input.Get<real>("Forcing","iso3D", 0);
    this->kmax = input.Get<real>("Forcing","iso3D", 1);

    std::vector<std::vector<real>> k3Disovec;

    int nxmax = kmax/kx0 + 1;
    int nymax = kmax/ky0 + 1;
    int nzmax = kmax/kz0 + 1;
    for (int nx=1; nx<nxmax; nx++) {
      for (int ny=1; ny<nymax; ny++) {
        for (int nz=1; nz<nzmax; nz++) {
          real kx = kx0*nx;
          real ky = ky0*ny;
          real kz = kz0*nz;
          real k_2 = kx*kx + ky*ky + kz*kz;
          if (k_2 >= kmin*kmin and k_2 <= kmax*kmax) {
            nForcingModes ++;
            k3Disovec.push_back({kx, ky, kz});
            modeNames.push_back({"I" + std::to_string(nx) + std::to_string(ny) + std::to_string(nz), "J" + std::to_string(nx) + std::to_string(ny) + std::to_string(nz), "K" + std::to_string(nx) + std::to_string(ny) + std::to_string(nz)});
            nForcingModes ++;
            k3Disovec.push_back({kx, -ky, kz});
            modeNames.push_back({"I" + std::to_string(nx) + std::to_string(-ny) + std::to_string(nz), "J" + std::to_string(nx) + std::to_string(-ny) + std::to_string(nz), "K" + std::to_string(nx) + std::to_string(-ny) + std::to_string(nz)});
            nForcingModes ++;
            k3Disovec.push_back({kx, ky, -kz});
            modeNames.push_back({"I" + std::to_string(nx) + std::to_string(ny) + std::to_string(-nz), "J" + std::to_string(nx) + std::to_string(ny) + std::to_string(-nz), "K" + std::to_string(nx) + std::to_string(ny) + std::to_string(-nz)});
            nForcingModes ++;
            k3Disovec.push_back({kx, -ky, -kz});
            modeNames.push_back({"I" + std::to_string(nx) + std::to_string(-ny) + std::to_string(-nz), "J" + std::to_string(nx) + std::to_string(-ny) + std::to_string(-nz), "K" + std::to_string(nx) + std::to_string(-ny) + std::to_string(-nz)});
          }
        }
      }
    }
    k3DisoHost = IdefixHostArray2D<real>("k3DisoHost", nForcingModes, 3);
    k3Diso = IdefixArray2D<real>("k3Diso", nForcingModes, 3);
    for (int l=0; l<nForcingModes; l++) {
      k3DisoHost(l,0) = k3Disovec[l][0];
      k3DisoHost(l,1) = k3Disovec[l][1];
      k3DisoHost(l,2) = k3Disovec[l][2];
    }
    Kokkos::deep_copy(k3Diso, k3DisoHost);
    EXPAND(
    this->forcingModesIdir = IdefixArray4D<Kokkos::complex<real>>("forcingModesIdir", nForcingModes,
                              data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);          ,
    this->forcingModesJdir = IdefixArray4D<Kokkos::complex<real>>("forcingModesJdir", nForcingModes,
                              data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);          ,
    this->forcingModesKdir = IdefixArray4D<Kokkos::complex<real>>("forcingModesKdir", nForcingModes,
                              data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);          )
  }

  else if (input.CheckEntry("Forcing", "iso2D")>=0) {
    this->forcingType = iso2D;
    this->normal2DisoStr = input.Get<std::string>("Forcing","iso2D", 0);
    if (this->normal2DisoStr == "IDIR") this->normal2Diso = IDIR;
    else if (this->normal2DisoStr == "JDIR") this->normal2Diso = JDIR;
    else if (this->normal2DisoStr == "KDIR") this->normal2Diso = KDIR;
    else IDEFIX_ERROR("The normal component for 2D isotropic forcing cannot not be something else than IDIR, JDIR or KDIR");
    #if COMPONENTS < 2 or DIMENSIONS < 2
      IDEFIX_ERROR("You cannot have 2D isotropic forcing with less than 2 components and dimensions.");
    #endif //COMPONENTS < 2 or DIMENSIONS < 2
    #if GEOMETRY == SPHERICAL
      if (this->normal2DisoStr == "JDIR" or this->normal2DisoStr == "KDIR") IDEFIX_ERROR("Comparing radial with angular wavenumbers is somewhat shady.");
    #endif //GEOMETRY == SPHERICAL

    this->kmin = input.Get<real>("Forcing","iso2D", 1);
    this->kmax = input.Get<real>("Forcing","iso2D", 2);

    std::vector<std::vector<real>> k2Disovec;

    kx0 = (this->normal2Diso==IDIR) ? ZERO_F : kx0;
    ky0 = (this->normal2Diso==JDIR) ? ZERO_F : ky0;
    kz0 = (this->normal2Diso==KDIR) ? ZERO_F : kz0;
    int nxmax = (this->normal2Diso==IDIR) ? 1 : kmax/kx0 + 1;
    int nymax = (this->normal2Diso==JDIR) ? 1 : kmax/ky0 + 1;
    int nzmax = (this->normal2Diso==KDIR) ? 1 : kmax/kz0 + 1;
    int xsign, ysign, zsign;
    if (this->normal2Diso==IDIR) {
      xsign = 0;
      ysign = 1;
      zsign = -1;
    } else if (this->normal2Diso==JDIR) {
      xsign = 1;
      ysign = 0;
      zsign = -1;
    } else {
      xsign = 1;
      ysign = -1;
      zsign = 0;
    }
    for (int nx=0; nx<nxmax; nx++) {
      for (int ny=0; ny<nymax; ny++) {
        for (int nz=0; nz<nzmax; nz++) {
          real kx = kx0*nx;
          real ky = ky0*ny;
          real kz = kz0*nz;
          real k_2 = kx*kx + ky*ky + kz*kz;
          if (k_2 >= kmin*kmin and k_2 <= kmax*kmax) {
            nForcingModes ++;
            k2Disovec.push_back({kx, ky, kz});
            modeNames.push_back({"I" + std::to_string(nx) + std::to_string(ny) + std::to_string(nz), "J" + std::to_string(nx) + std::to_string(ny) + std::to_string(nz), "K" + std::to_string(nx) + std::to_string(ny) + std::to_string(nz)});
            nForcingModes ++;
            k2Disovec.push_back({xsign*kx, ysign*ky, zsign*kz});
            modeNames.push_back({"I" + std::to_string(xsign*nx) + std::to_string(ysign*ny) + std::to_string(zsign*nz), "J" + std::to_string(xsign*nx) + std::to_string(ysign*ny) + std::to_string(zsign*nz), "K" + std::to_string(xsign*nx) + std::to_string(ysign*ny) + std::to_string(zsign*nz)});
          }
        }
      }
    }
    k2DisoHost = IdefixHostArray2D<real>("k2DisoHost", nForcingModes, 3);
    k2Diso = IdefixArray2D<real>("k2Diso", nForcingModes, 3);
    for (int l=0; l<nForcingModes; l++) {
      k2DisoHost(l,0) = k2Disovec[l][0];
      k2DisoHost(l,1) = k2Disovec[l][1];
      k2DisoHost(l,2) = k2Disovec[l][2];
    }
    Kokkos::deep_copy(k2Diso, k2DisoHost);
    EXPAND(
    this->forcingModesIdir = IdefixArray4D<Kokkos::complex<real>>("forcingModesIdir", nForcingModes,
                              data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);          ,
    this->forcingModesJdir = IdefixArray4D<Kokkos::complex<real>>("forcingModesJdir", nForcingModes,
                              data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);          ,
    this->forcingModesKdir = IdefixArray4D<Kokkos::complex<real>>("forcingModesKdir", nForcingModes,
                              data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);          )
  }

  else if (input.CheckEntry("Forcing", "ani3D")>=0) {
    #if COMPONENTS < 3 or DIMENSIONS < 3
      IDEFIX_ERROR("You cannot have 3D anisotropic forcing with less than 3 components and dimensions.");
    #endif //COMPONENTS < 3 or DIMENSIONS < 3
    this->forcingType = ani3D;
    this->normal3DaniStr = input.Get<std::string>("Forcing","ani3D", 0);
    if (this->normal3DaniStr == "IDIR") this->normal3Dani = IDIR;
    else if (this->normal3DaniStr == "JDIR") this->normal3Dani = JDIR;
    else if (this->normal3DaniStr == "KDIR") this->normal3Dani = KDIR;
    else IDEFIX_ERROR("The normal component for 3D anisotropic forcing cannot not be something else than IDIR, JDIR or KDIR");
//    this->normal3DaniBoundStr = input.Get<std::string>("Forcing","ani3D", 1);
//    if (this->normal3DaniBoundStr == "bothfree") this->normal3DaniBound = bothFree;
//    else if (this->normal3DaniBoundStr == "righthomdir") this->normal3DaniBound = rightHomDir;
//    else if (this->normal3DaniBoundStr == "bothhomdir") this->normal3DaniBound = bothHomDir;
//    else if (this->normal3DaniBoundStr == "bothhomneu") this->normal3DaniBound = bothHomNeu;
//    else IDEFIX_ERROR("The boundaries for the normal component of the forcing can only be bothfree, righthomdir, bothhomdir or bothhomneu for now.");
//    this->normal3DaniBasisStr = input.Get<std::string>("Forcing","ani3D", 2);
//    if (this->normal3DaniBasisStr == "chebyshev") this->normal3DaniBasis = chebyshev;
////    else if (this->normal3DaniBasisStr == "legendre") this->normal3DaniBasis = legendre;
//    else if (this->normal3DaniBasisStr == "fourier") this->normal3DaniBasis = fourier;
//    else IDEFIX_ERROR("The basis for the normal component of the forcing can only be chebyshev and fourier for now (legendre coming soon).");
////    if (this->normal3DaniBasis = chebyshev and this->normal3DaniBound = bothHomNeu) IDEFIX_ERROR("Cannot have a chebyshev basis which has both homogeneous Neumann boundaries.");
////    else if (normal3DaniBasis = fourier and normal3DaniBound = rightHomDir) IDEFIX_ERROR("Cannot have a fourier basis which has only right homogeneous Dirichlet boundary.");

    this->normal3DaniBasisStr = input.Get<std::string>("Forcing","ani3D", 1);
    if (this->normal3DaniBasisStr == "chebyshev") this->normal3DaniBasis = chebyshev;
//    else if (this->normal3DaniBasisStr == "legendre") this->normal3DaniBasis = legendre;
    else if (this->normal3DaniBasisStr == "fourier") this->normal3DaniBasis = fourier;
    else IDEFIX_ERROR("The basis for the normal component of the forcing can only be chebyshev and fourier for now (legendre coming soon).");

    std::vector<NormalBoundType> vecNormal3DaniBound;
    std::string dirbegStr = "X" + std::to_string(this->normal3Dani + 1) + "-beg";
    std::string direndStr = "X" + std::to_string(this->normal3Dani + 1) + "-end";
    if (input.Get<std::string>("Boundary",dirbegStr, 0) != input.Get<std::string>("Boundary",direndStr, 0)) {
      IDEFIX_ERROR("Chebyshev basis with different BCs at the two sides are not yet implemented");
    } else if (input.Get<std::string>("Boundary",dirbegStr, 0) == "reflective") {
      (this->normal3Dani == IDIR) ? vecNormal3DaniBound.push_back(bothHomDir) : vecNormal3DaniBound.push_back(bothHomNeu);
      (this->normal3Dani == JDIR) ? vecNormal3DaniBound.push_back(bothHomDir) : vecNormal3DaniBound.push_back(bothHomNeu);
      (this->normal3Dani == KDIR) ? vecNormal3DaniBound.push_back(bothHomDir) : vecNormal3DaniBound.push_back(bothHomNeu);
    } else if (input.Get<std::string>("Boundary",dirbegStr, 0) == "outflow") {
      vecNormal3DaniBound.push_back(bothHomNeu);
      vecNormal3DaniBound.push_back(bothHomNeu);
      vecNormal3DaniBound.push_back(bothHomNeu);
    } else if (input.Get<std::string>("Boundary",dirbegStr, 0) == "userdef") {
      for (int dir = IDIR; dir<COMPONENTS; dir++) {
        std::string userdefDirbegStr = "begBoundTypeVX" + std::to_string(dir + 1);
        std::string userdefDirendStr = "endBoundTypeVX" + std::to_string(dir + 1);
        if (input.Get<std::string>("Boundary",userdefDirbegStr, 0) != input.Get<std::string>("Boundary",userdefDirendStr, 0)) {
          IDEFIX_ERROR("Chebyshev basis with different BCs at the two sides are not yet implemented");
        } else if (input.Get<std::string>("Boundary",userdefDirbegStr, 0) == "dirichletZero") {
          vecNormal3DaniBound.push_back(bothHomDir);
        } else if (input.Get<std::string>("Boundary",userdefDirbegStr, 0) == "neumann") {
          vecNormal3DaniBound.push_back(bothHomNeu);
        } else {
          IDEFIX_ERROR("userdef BCs not recognised, they are needed choose the suitable forcing basis.");
        }
      }
    } else if (input.Get<std::string>("Boundary",dirbegStr, 0) == "periodic") {
      IDEFIX_ERROR("periodic BCs recognised, you should use 3D isotropic forcing in this direction.");
    } else {
      IDEFIX_ERROR("BCs not recognised, they are needed choose the suitable forcing basis.");
    }
    normal3DaniBoundHost = IdefixHostArray1D<NormalBoundType>("normal3DaniBoundHost", COMPONENTS);
    normal3DaniBound = IdefixArray1D<NormalBoundType>("normal3DaniBound", COMPONENTS);
    for (int l=IDIR; l<COMPONENTS; l++) {
      normal3DaniBoundHost(l) = vecNormal3DaniBound[l];
      normal3DaniBoundHost(l) = vecNormal3DaniBound[l];
      normal3DaniBoundHost(l) = vecNormal3DaniBound[l];
    }
    Kokkos::deep_copy(normal3DaniBound, normal3DaniBoundHost);

    this->kmin = input.Get<real>("Forcing","ani3D", 2);
    this->kmax = input.Get<real>("Forcing","ani3D", 3);

    if (input.GetOrSet<std::string>("Forcing","ani3D", 4, "n") == "write") WriteNormalBasis(folder);

    std::vector<std::vector<real>> k3Danivec;

    int nxmax = kmax/kx0 + 1;
    int nymax = kmax/ky0 + 1;
    int nzmax = kmax/kz0 + 1;
    int xsign, ysign, zsign;
    if (this->normal3Dani==IDIR) {
      xsign = 1;
      ysign = 1;
      zsign = -1;
    } else if (this->normal3Dani==JDIR) {
      xsign = 1;
      ysign = 1;
      zsign = -1;
    } else if (this->normal3Dani==KDIR) {
      xsign = 1;
      ysign = -1;
      zsign = 1;
    } else {
      IDEFIX_ERROR("Forcing: normal direction not known");
    }
    for (int nx=1; nx<nxmax; nx++) {
      for (int ny=1; ny<nymax; ny++) {
        for (int nz=1; nz<nzmax; nz++) {
          real kx = kx0*nx;
          real ky = ky0*ny;
          real kz = kz0*nz;
          real factor = 1.;
          #if GEOMETRY == SPHERICAL
            factor *= 4./pow(xbeg+xend, 2.);
          #endif //GEOMETRY == SPHERICAL
          //undimensionalise kx so that kmin and kmax corresponds to angular ktheta and kphi
          real k_2 = kx*kx/factor + (ky*ky + kz*kz);
//          real k_2 = kx*kx + factor*(ky*ky + kz*kz);
          if (k_2 >= kmin*kmin and k_2 <= kmax*kmax) {
            nForcingModes ++;
            k3Danivec.push_back({kx, ky, kz});
            modeNames.push_back({"I" + std::to_string(nx) + std::to_string(ny) + std::to_string(nz), "J" + std::to_string(nx) + std::to_string(ny) + std::to_string(nz), "K" + std::to_string(nx) + std::to_string(ny) + std::to_string(nz)});
            nForcingModes ++;
            k3Danivec.push_back({xsign*kx, ysign*ky, zsign*kz});
            modeNames.push_back({"I" + std::to_string(xsign*nx) + std::to_string(ysign*ny) + std::to_string(zsign*nz), "J" + std::to_string(xsign*nx) + std::to_string(ysign*ny) + std::to_string(zsign*nz), "K" + std::to_string(xsign*nx) + std::to_string(ysign*ny) + std::to_string(zsign*nz)});
          }
        }
      }
    }
    k3DaniHost = IdefixHostArray2D<real>("k3DaniHost", nForcingModes, 3);
    k3Dani = IdefixArray2D<real>("k3Dani", nForcingModes, 3);
    for (int l=0; l<nForcingModes; l++) {
      k3DaniHost(l,0) = k3Danivec[l][0];
      k3DaniHost(l,1) = k3Danivec[l][1];
      k3DaniHost(l,2) = k3Danivec[l][2];
    }
    Kokkos::deep_copy(k3Dani, k3DaniHost);
    EXPAND(
    this->forcingModesIdir = IdefixArray4D<Kokkos::complex<real>>("forcingModesIdir", nForcingModes,
                              data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);          ,
    this->forcingModesJdir = IdefixArray4D<Kokkos::complex<real>>("forcingModesJdir", nForcingModes,
                              data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);          ,
    this->forcingModesKdir = IdefixArray4D<Kokkos::complex<real>>("forcingModesKdir", nForcingModes,
                              data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);          )
  }

  else if (input.CheckEntry("Forcing", "vsh")>=0) {
    this->forcingType = vsh;
    #if GEOMETRY != SPHERICAL
      IDEFIX_ERROR("You cannot have vsh forcing in another geometry than spherical.");
    #endif //GEOMETRY != SPHERICAL
    #if VSH == NO
      IDEFIX_ERROR("You cannot have vsh forcing without the VSH module.");
    #endif //VSH == NO
// if () IDEFIX_ERROR("You cannot have vsh forcing without the full sphere.");
    this->ellmin = input.Get<int>("Forcing","vsh", 0);
    this->ellmax = input.Get<int>("Forcing","vsh", 1);
    this->mmin = input.GetOrSet<int>("Forcing","vsh", 2, 0);
    this->mmax = input.GetOrSet<int>("Forcing","vsh", 3, this->ellmax);

    std::vector<std::vector<int>> ellmVshvec;
    for (int ell=ellmin; ell<ellmax; ell++) {
      for (int m=mmin; m<mmax & m<ell+1 ; m++) {
        this->nForcingModes++;
        ellmVshvec.push_back({ell,m});
        modeNames.push_back({"Y" + std::to_string(ell) + std::to_string(m), "S" + std::to_string(ell) + std::to_string(m), "T" + std::to_string(ell) + std::to_string(m)});
      }
    }
    ellmVshHost = IdefixHostArray2D<int>("ellmVshHost", nForcingModes, 2);
    ellmVsh = IdefixArray2D<int>("ellmVsh", nForcingModes, 2);
    for (int l=0; l<nForcingModes; l++) {
      ellmVshHost(l,0) = ellmVshvec[l][0];
      ellmVshHost(l,1) = ellmVshvec[l][1];
    }
    Kokkos::deep_copy(ellmVsh, ellmVshHost);

    EXPAND(
    this->forcingModesIdir = IdefixArray4D<Kokkos::complex<real>>("forcingModesIdir", nForcingModes,
                              data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);            ,
    // WARNING, vsh case needs a special treatment since the forcing modes between
    // the theta and phi directions are coupled to each other...
    this->forcingModesJdir = IdefixArray4D<Kokkos::complex<real>>("forcingModesJdir", 2*nForcingModes,
                              data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);            ,
    this->forcingModesKdir = IdefixArray4D<Kokkos::complex<real>>("forcingModesKdir", 2*nForcingModes,
                              data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);            )
  }

//  else { this->haveUserDef = true;}
  else { this->forcingType = userDef;} // TO BE CODED THOUGH

  if (input.CheckEntry("Forcing", "solenoidal")>=0) {
    #if GEOMETRY == POLAR or GEOMETRY == CYLINDRICAL
      IDEFIX_ERROR("Cannot have solenoidal forcing in polar in cylindrical geometry.");
    #endif //GEOMETRY == POLAR or GEOMETRY == CYLINDRICAL
    if (this->forcingType == vsh) IDEFIX_ERROR("Cannot have solenoidal forcing with vsh forcing. Select yourself the toroidal vector spherical harmonics in this case.");
    else if (this->forcingType == userDef) IDEFIX_ERROR("Cannot have solenoidal forcing with userdef forcing.");
    else this->haveSolenoidalForcing = true;
  }

  // Allocate required arrays
  this->forcingTerm = IdefixArray4D<real>("forcingTerm", COMPONENTS,
                              data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
  this->pristineForcingTerm = IdefixArray4D<real>("pristineForcingTerm", COMPONENTS,
                              data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
  this->solenoidalForcingTerm = IdefixArray4D<real>("solenoidalForcingTerm", COMPONENTS,
                              data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
//  this->compressiveForcingTerm = IdefixArray4D<real>("compressiveForcingTerm", COMPONENTS,
//                              data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);

  this->means = IdefixArray2D<real>("means", nForcingModes, COMPONENTS);
  this->tcorrs = IdefixArray2D<real>("tcorrs", nForcingModes, COMPONENTS);
  this->epsilons = IdefixArray2D<real>("epsilons", nForcingModes, COMPONENTS);

// WARNING DO SOMETHING WITH MEANS TCORRS AND EPSILONS
  this->hostMeans = IdefixHostArray2D<real>("hostMeans", nForcingModes, COMPONENTS);
  this->hostTcorrs = IdefixHostArray2D<real>("hostTcorrs", nForcingModes, COMPONENTS);
  this->hostEpsilons = IdefixHostArray2D<real>("hostEpsilons", nForcingModes, COMPONENTS);
  this->machNumber = input.GetOrSet<real>("Forcing", "mach", 0, -1.);
  if (this->machNumber >= 0) { // in this case they will be reset later
    tcorr = 1.;
    epsilon = -1.;
  } else{
    tcorr = input.Get<real>("Forcing", "t_corr", 0);
    epsilon = input.Get<real>("Forcing", "epsilon", 0);
  }
  #ifdef ISOTHERMAL
    this->cs = input.Get<real>("Hydro", "csiso", 1); //WARNING DOESN'T WORK FOR USERDEF csiso
  #endif

  this->oUprocesses.InitProcesses(this->folder, this->seed, this->nForcingModes, this->modeNames);
  idfx::popRegion();
}

void Forcing::ShowConfig() {
  idfx::cout << "Forcing: ENABLED with seed " << seed << "." << std::endl;
  if (stopTime < std::numeric_limits<double>::infinity()) {
    idfx::cout << "Forcing: will be stopped at t=" << stopTime << " ." << std::endl;
  }
  switch(forcingType) {
    case ForcingType::iso3D:
      idfx::cout << "Forcing: 3D isotropic." << std::endl;
      idfx::cout << "Forcing: kmin=" << kmin << " and kmax=" << kmax << " ." << std::endl;
      idfx::cout << "Forcing: There are " << nForcingModes << " different forcing modes." << std::endl;
      if (haveSolenoidalForcing) idfx::cout << "Forcing: solenoidal." << std::endl;
      break;
    case ForcingType::iso2D:
      idfx::cout << "Forcing: 2D isotropic with normal " << normal2DisoStr << "." << std::endl;
      idfx::cout << "Forcing: kmin=" << kmin << " and kmax=" << kmax << " ." << std::endl;
      idfx::cout << "Forcing: There are " << nForcingModes << " different forcing modes." << std::endl;
      if (haveSolenoidalForcing) idfx::cout << "Forcing: solenoidal." << std::endl;
      break;
    case ForcingType::ani3D:
      idfx::cout << "Forcing: 3D anisotropic with normal " << normal3DaniStr << " and " << normal3DaniBasisStr << " basis." << std::endl;
      idfx::cout << "Forcing: kmin=" << kmin << " and kmax=" << kmax << " ." << std::endl;
      idfx::cout << "Forcing: There are " << nForcingModes << " different forcing modes." << std::endl;
      if (haveSolenoidalForcing) idfx::cout << "Forcing: solenoidal." << std::endl;
      break;
    case ForcingType::vsh:
      idfx::cout << "Forcing: vector spherical harmonics." << std::endl;
      idfx::cout << "Forcing: ellmin=" << ellmin << " and ellmax=" << ellmax << " ." << std::endl;
      idfx::cout << "Forcing: There are " << nForcingModes << " different forcing modes." << std::endl;
      break;
    case ForcingType::userDef:
      idfx::cout << "Forcing: userdef." << std::endl;
      break;
  }

  if (machNumber >= ZERO_F) {
    idfx::cout << "Forcing: Mach number=" << machNumber << " ." << std::endl;
  } else {
    idfx::cout << "Forcing: epsilon=" << epsilon << " and tcorr=" << tcorr << " ." << std::endl;
  }

//    if(skipGravity>1) {
//      idfx::cout << "Gravity: gravity field will be updated every " << skipGravity
//                 << " cycles." << std::endl;
//    }
}

void Forcing::InitForcingParameters() {
  idfx::pushRegion("Forcing::InitForcingParameters");

  if (this->machNumber >= 0) {
    real kf = HALF_F*(kmin+kmax);
    #if HAVE_ENERGY 
      ComputeAverageSoundSpeed();
    #endif //ISOTHERMAL
    this->epsilon = pow(this->machNumber*cs,3.)*kf/(2.*M_PI);
    this->tcorr = 2.*M_PI/(this->machNumber*cs*kf);
  }
  for (int l=0; l<nForcingModes; l++) {
    for (int dir=IDIR; dir<COMPONENTS; dir++) {
      this->hostMeans(l,dir) = ZERO_F;
      this->hostTcorrs(l,dir) = tcorr;
      this->hostEpsilons(l,dir) = epsilon;
    }
  }
  Kokkos::deep_copy(this->means, this->hostMeans);
  Kokkos::deep_copy(this->tcorrs, this->hostTcorrs);
  Kokkos::deep_copy(this->epsilons, this->hostEpsilons);
  this->oUprocesses.SetProcesses(this->means, this->tcorrs, this->epsilons);

  idfx::popRegion();
}

void Forcing::InitForcingModes() {
  idfx::pushRegion("Forcing::InitForcingModes");

  int normal3Dani = this->normal3Dani;
//  IDEFIX_ERROR("DSTOP");
  IdefixArray1D<NormalBoundType> normal3DaniBound = this->normal3DaniBound;
  int normal3DaniBasis = this->normal3DaniBasis;
  real kx0 = this->kx0;
  real ky0 = this->ky0;
  real kz0 = this->kz0;
  real xbeg = this->xbeg;
  real ybeg = this->ybeg;
  real zbeg = this->zbeg;
  real xend = this->xend;
  real yend = this->yend;
  real zend = this->zend;
  EXPAND(
  IdefixArray4D<Kokkos::complex<real>> forcingModesIdir = this->forcingModesIdir; ,
  IdefixArray4D<Kokkos::complex<real>> forcingModesJdir = this->forcingModesJdir; ,
  IdefixArray4D<Kokkos::complex<real>> forcingModesKdir = this->forcingModesKdir; )
  IdefixArray1D<real> x1 = data->x[IDIR];
  IdefixArray1D<real> x2 = data->x[JDIR];
  IdefixArray1D<real> x3 = data->x[KDIR];
  IdefixArray2D<real> k3Diso = this->k3Diso;
  IdefixArray2D<real> k2Diso = this->k2Diso;
  IdefixArray2D<real> k3Dani = this->k3Dani;
  IdefixArray2D<int> ellmVsh = this->ellmVsh;
  Kokkos::complex unit_j(0.,1.);
  int nForcingModes = this->nForcingModes;
  ForcingType forcingType = this->forcingType;
  idefix_for("Forcing::InitForcingModes",
              0, nForcingModes,
              0, data->np_tot[KDIR],
              0, data->np_tot[JDIR],
              0, data->np_tot[IDIR],
              KOKKOS_LAMBDA(int l, int k, int j, int i) {
                EXPAND(
                forcingModesIdir(l,k,j,i) = ZERO_F; ,
                forcingModesJdir(l,k,j,i) = ZERO_F; ,
                forcingModesKdir(l,k,j,i) = ZERO_F; )
                if (forcingType == vsh) {
                  EXPAND(                                           ,
                  forcingModesJdir(l+nForcingModes,k,j,i) = ZERO_F; ,
                  forcingModesKdir(l+nForcingModes,k,j,i) = ZERO_F; )
                }
              });
  switch(forcingType) {
    case ForcingType::iso3D:
      idefix_for("iso3D", 0, nForcingModes, 0, data->np_tot[KDIR], 0, data->np_tot[JDIR], 0, data->np_tot[IDIR],
                  KOKKOS_LAMBDA (int l, int k, int j, int i) {
                    real kx = k3Diso(l, IDIR);
                    real ky = k3Diso(l, JDIR);
                    real kz = k3Diso(l, KDIR);
                    real kdotx = kx*x1(i) + ky*x2(j) + kz*x3(k);
                    EXPAND(
                    forcingModesIdir(l,k,j,i) = exp(unit_j * kdotx); ,
                    forcingModesJdir(l,k,j,i) = exp(unit_j * kdotx); ,
                    forcingModesKdir(l,k,j,i) = exp(unit_j * kdotx); )
      });
      break;
    case ForcingType::iso2D:
      idefix_for("iso2D", 0, nForcingModes, 0, data->np_tot[KDIR], 0, data->np_tot[JDIR], 0, data->np_tot[IDIR],
                  KOKKOS_LAMBDA (int l, int k, int j, int i) {
                    real kx = k2Diso(l, IDIR);
                    real ky = k2Diso(l, JDIR);
                    real kz = k2Diso(l, KDIR);
                    real kdotx = kx*x1(i) + ky*x2(j) + kz*x3(k);
                    EXPAND(
                    forcingModesIdir(l,k,j,i) = exp(unit_j * kdotx); ,
                    forcingModesJdir(l,k,j,i) = exp(unit_j * kdotx); ,
                    forcingModesKdir(l,k,j,i) = exp(unit_j * kdotx); )
      });
      break;
    case ForcingType::ani3D:
      idefix_for("ani3D", 0, nForcingModes, 0, data->np_tot[KDIR], 0, data->np_tot[JDIR], 0, data->np_tot[IDIR],
                  KOKKOS_LAMBDA (int l, int k, int j, int i) {
                    real kx = k3Dani(l, IDIR);
                    real ky = k3Dani(l, JDIR);
                    real kz = k3Dani(l, KDIR);
                    int offx = (normal3Dani==IDIR) ? 0 : 1;
                    int offy = (normal3Dani==JDIR) ? 0 : 1;
                    int offz = (normal3Dani==KDIR) ? 0 : 1;
                    real kdotx = offx*kx*x1(i) + offy*ky*x2(j) + offz*kz*x3(k);
                    int oppx = (normal3Dani==IDIR) ? 1 : 0;
                    int oppy = (normal3Dani==JDIR) ? 1 : 0;
                    int oppz = (normal3Dani==KDIR) ? 1 : 0;
                    int order = 0;
                    order = (normal3Dani==IDIR) ? kx/kx0 : order;
                    order = (normal3Dani==JDIR) ? ky/ky0 : kz/kz0;
                    real rightx = oppx*x1(i) + oppy*x2(j) + oppz*x3(k);
                    real rightx0 = oppx*xbeg + oppy*ybeg + oppz*zbeg;
                    real rightx1 = oppx*xend + oppy*yend + oppz*zend;
                    switch(normal3DaniBasis) {
                      case NormalBasis::chebyshev:
                        for (int dir=IDIR; dir<COMPONENTS; dir++) {
                          int bound = normal3DaniBound(dir);
                          switch(bound) {
                            case NormalBoundType::bothHomDir:
                              if (dir==IDIR) forcingModesIdir(dir,k,j,i) = K_cheby_fst_both_hom_dir(order, K_aff_11(rightx, rightx0, rightx1)) * exp(unit_j*kdotx);
                              else if (dir==JDIR) forcingModesJdir(dir,k,j,i) = K_cheby_fst_both_hom_dir(order, K_aff_11(rightx, rightx0, rightx1)) * exp(unit_j*kdotx);
                              else if (dir==KDIR) forcingModesKdir(dir,k,j,i) = K_cheby_fst_both_hom_dir(order, K_aff_11(rightx, rightx0, rightx1)) * exp(unit_j*kdotx);
                            break;
                            case NormalBoundType::bothHomNeu:
                              if (dir==IDIR) forcingModesIdir(dir,k,j,i) = K_cheby_fst_both_hom_neu(order, K_aff_11(rightx, rightx0, rightx1)) * exp(unit_j*kdotx);
                              else if (dir==JDIR) forcingModesJdir(dir,k,j,i) = K_cheby_fst_both_hom_neu(order, K_aff_11(rightx, rightx0, rightx1)) * exp(unit_j*kdotx);
                              else if (dir==KDIR) forcingModesKdir(dir,k,j,i) = K_cheby_fst_both_hom_neu(order, K_aff_11(rightx, rightx0, rightx1)) * exp(unit_j*kdotx);
                            break;
                          }
                        }
                      break;
                      case NormalBasis::fourier:
                        for (int dir=IDIR; dir<COMPONENTS; dir++) {
                          int bound = normal3DaniBound(dir);
                          switch(bound) {
                            case NormalBoundType::bothHomDir:
                              if (dir==IDIR) forcingModesIdir(dir,k,j,i) = K_sin(order, K_aff_02pi(rightx, rightx0, rightx1)) * exp(unit_j*kdotx);
                              else if (dir==JDIR) forcingModesJdir(dir,k,j,i) = K_sin(order, K_aff_02pi(rightx, rightx0, rightx1)) * exp(unit_j*kdotx);
                              else if (dir==KDIR) forcingModesKdir(dir,k,j,i) = K_sin(order, K_aff_02pi(rightx, rightx0, rightx1)) * exp(unit_j*kdotx);
                            break;
                            case NormalBoundType::bothHomNeu:
                              if (dir==IDIR) forcingModesIdir(dir,k,j,i) = K_cos(order, K_aff_02pi(rightx, rightx0, rightx1)) * exp(unit_j*kdotx);
                              else if (dir==JDIR) forcingModesJdir(dir,k,j,i) = K_cos(order, K_aff_02pi(rightx, rightx0, rightx1)) * exp(unit_j*kdotx);
                              else if (dir==KDIR) forcingModesKdir(dir,k,j,i) = K_cos(order, K_aff_02pi(rightx, rightx0, rightx1)) * exp(unit_j*kdotx);
                            break;
                          }
                        }
                    }
      });
      break;
    #if VSH == YES
      case ForcingType::vsh:
        IdefixArray4D<Kokkos::complex<real>> Ylm_r = this->data->Ylm_r;
        IdefixArray4D<Kokkos::complex<real>> Slm_th = this->data->Slm_th;
        IdefixArray4D<Kokkos::complex<real>> Slm_phi = this->data->Slm_phi;
        IdefixArray4D<Kokkos::complex<real>> Tlm_th = this->data->Tlm_th;
        IdefixArray4D<Kokkos::complex<real>> Tlm_phi = this->data->Tlm_phi;
        int nForcingModes = this->nForcingModes;
        idefix_for("vsh", 0, nForcingModes, 0, data->np_tot[KDIR], 0, data->np_tot[JDIR], 0, data->np_tot[IDIR],
                    KOKKOS_LAMBDA (int l, int k, int j, int i) {
                      int ell = ellmVsh(l, 0);
                      int m = ellmVsh(l, 1);
                      forcingModesIdir(l,k,j,i) = Ylm_r(ell,m,k,j);
                      forcingModesJdir(l,k,j,i) = Slm_th(ell,m,k,j);
                      forcingModesJdir(l+nForcingModes,k,j,i) = Tlm_th(ell,m,k,j);
                      forcingModesKdir(l,k,j,i) = Slm_phi(ell,m,k,j);
                      forcingModesKdir(l+nForcingModes,k,j,i) = Tlm_phi(ell,m,k,j);
        });
        break;
    #endif //VSH == YES
//    case ForcingType::userDef:
//      break;
  }
  idfx::popRegion();
}

// This function writes the normal basis in a txt file to later plot them with matplotlib
// Working in mono-domain and in the X1 direction only so far
void Forcing::WriteNormalBasis(std::string folder) {
  IdefixArray1D<real> x1 = data->x[IDIR];
  int order_max = kmax/kx0 + 1;
  int order_min = kmin/kx0 - 1;
  order_min = std::max(1, order_min);
  real xbeg = this->xbeg;
  real xend = this->xend;
  IdefixArray1D<NormalBoundType> normal3DaniBound = this->normal3DaniBound;
  int normal3DaniBasis = this->normal3DaniBasis;
  IdefixHostArray3D<real> normalBasisHost("normalBasisHost", order_max, COMPONENTS, data->np_tot[IDIR]);
  IdefixArray3D<real> normalBasis("normalBasis", order_max, COMPONENTS, data->np_tot[IDIR]);
  idefix_for("write", 0, order_max, 0, data->np_tot[IDIR],
              KOKKOS_LAMBDA (int order, int i) {
                switch(normal3DaniBasis) {
                  case NormalBasis::chebyshev:
                    for (int dir = IDIR; dir < COMPONENTS; dir++) {
                      int bound = normal3DaniBound(dir);
                      switch(bound) {
                        case NormalBoundType::bothHomDir:
                          normalBasis(order,dir,i) = K_cheby_fst_both_hom_dir(order, K_aff_11(x1(i), xbeg, xend));
                        break;
                        case NormalBoundType::bothHomNeu:
                          normalBasis(order,dir,i) = K_cheby_fst_both_hom_neu(order, K_aff_11(x1(i), xbeg, xend));
                        break;
                      }
                    }
                    break;
                  case NormalBasis::fourier:
                    for (int dir = IDIR; dir < COMPONENTS; dir++) {
                      int bound = normal3DaniBound(dir);
                      switch(bound) {
                        case NormalBoundType::bothHomDir:
                          normalBasis(order,dir,i) = K_sin(order, K_aff_02pi(x1(i), xbeg, xend));
                        break;
                        case NormalBoundType::bothHomNeu:
                          normalBasis(order,dir,i) = K_cos(order, K_aff_02pi(x1(i), xbeg, xend));
                        break;
                      }
                    }
                    break;
                  }
  });
  IdefixHostArray1D<real> x1host("x1host", data->np_tot[IDIR]);
  Kokkos::deep_copy(normalBasisHost, normalBasis);
  Kokkos::deep_copy(x1host, x1);
  if(idfx::prank==0) {
    int precision = 10;
    int col_width = precision + 10;
    std::string filename = folder + "/normalBasis.dat";
    std::ofstream file;
    file.open(filename, std::ios::trunc);
    file << std::setw(col_width) << "x";
    for (int dir=IDIR; dir<COMPONENTS; dir++) {
      for (int order=order_min; order<order_max; order++) {
        std::string current_name;
        if (dir == IDIR) current_name = "fI"+std::to_string(order);
        else if (dir == JDIR) current_name = "fJ"+std::to_string(order);
        else if (dir == KDIR) current_name = "fK"+std::to_string(order);
        else IDEFIX_ERROR("Fucking error");
        file << std::setw(col_width) << current_name;
      }
    }
    file << std::endl;
    file.precision(precision);
    for (int i=0; i<data->np_tot[IDIR]; i++) {
      file << std::setw(col_width) << x1host(i);
      for (int dir=IDIR; dir<COMPONENTS; dir++) {
        for (int order=order_min; order<order_max; order++) {
          file << std::scientific << std::setw(col_width) << normalBasisHost(order,dir,i);
        }
      }
      file << std::endl;
    }
    file << std::endl;
    file.close();
  }
}

// This function compute the required forcing field
void Forcing::ComputeForcing(real dt) {
  idfx::pushRegion("Forcing::ComputeForcing");
  
  ResetForcingTerms();
  if (haveSolenoidalForcing) {
    ComputePristineForcing(dt);
    ComputeSolenoidalForcing(dt);
    forcingTerm = solenoidalForcingTerm;
  } else {
    ComputePristineForcing(dt);
    forcingTerm = pristineForcingTerm;
  }
  oUprocesses.UpdateProcessesValues(dt);
  idfx::popRegion();
}

// This function compute the pristine forcing field
void Forcing::ComputePristineForcing(real dt) {
  idfx::pushRegion("Forcing::ComputePristineForcing");
  
  EXPAND(
  IdefixArray4D<Kokkos::complex<real>> forcingModesIdir = this->forcingModesIdir; ,
  IdefixArray4D<Kokkos::complex<real>> forcingModesJdir = this->forcingModesJdir; ,
  IdefixArray4D<Kokkos::complex<real>> forcingModesKdir = this->forcingModesKdir; )
  IdefixArray4D<real> pristineForcingTerm = this->pristineForcingTerm;
  IdefixArray2D<Kokkos::complex<real>> ouValues = this->oUprocesses.ouValues;
  int nForcingModes = this->nForcingModes;
  idefix_for("ComputePristineForcing", 0, data->np_tot[KDIR], 0, data->np_tot[JDIR], 0, data->np_tot[IDIR],
              KOKKOS_LAMBDA (int k, int j, int i) {
                  Kokkos::complex<real> forcingX1=0.;
                  Kokkos::complex<real> forcingX2=0.;
                  Kokkos::complex<real> forcingX3=0.;
                  for (int l=0; l<nForcingModes; l++) {
                    forcingX1 += ouValues(l,IDIR)*forcingModesIdir(l,k,j,i);
                    #if COMPONENTS >= 2
                      forcingX2 += ouValues(l,JDIR)*forcingModesJdir(l,k,j,i);
                      #if VSH == YES
                        //in this case (l,JDIR) -> (l,Slm) and (l,KDIR) -> (l,Tlm)
                        forcingX2 += ouValues(l,KDIR)*forcingModesJdir(l+nForcingModes,k,j,i);
                      #endif //VSH == YES
                      #if COMPONENTS == 3
                        #if VSH == YES
                          //don't forget (l,JDIR) -> (l,Slm) and (l,KDIR) -> (l,Tlm)
                          forcingX3 += ouValues(l,JDIR)*forcingModesKdir(l,k,j,i);
                          forcingX3 += ouValues(l,KDIR)*forcingModesJdir(l+nForcingModes,k,j,i);
                        #else
                          forcingX3 += ouValues(l,KDIR)*forcingModesKdir(l,k,j,i);
                        #endif //VSH == YES
                      #endif //COMPONENTS == 3
                    #endif //COMPONENTS >= 2
                  }
                  pristineForcingTerm(IDIR,k,j,i) += forcingX1.real();
                  pristineForcingTerm(JDIR,k,j,i) += forcingX2.real();
                  pristineForcingTerm(KDIR,k,j,i) += forcingX3.real();
  });

  idfx::popRegion();
}

// This function compute the solenoidal part of the forcing field
void Forcing::ComputeSolenoidalForcing(real dt) {
  idfx::pushRegion("Forcing::ComputeSolenoidalForcing");
  
  IdefixArray2D<real> kiso = (this->forcingType == iso3D) ? this->k3Diso : this->k2Diso;
  EXPAND(
  IdefixArray4D<Kokkos::complex<real>> forcingModesIdir = this->forcingModesIdir; ,
  IdefixArray4D<Kokkos::complex<real>> forcingModesJdir = this->forcingModesJdir; ,
  IdefixArray4D<Kokkos::complex<real>> forcingModesKdir = this->forcingModesKdir; )
  IdefixArray4D<real> pristineForcingTerm = this->pristineForcingTerm;
  IdefixArray4D<real> solenoidalForcingTerm = this->solenoidalForcingTerm;

  IdefixArray1D<real> dx1 = data->dx[IDIR];
  IdefixArray1D<real> dx2 = data->dx[JDIR];
  IdefixArray1D<real> dx3 = data->dx[KDIR];

  IdefixArray1D<real> x1m = data->xl[IDIR];
  IdefixArray1D<real> x2m = data->xl[JDIR];
  IdefixArray1D<real> x1 = data->x[IDIR];
  IdefixArray1D<real> x2 = data->x[JDIR];

  IdefixArray1D<real> sinx2 = data->sinx2;

  real kmean = HALF_F*(kmax + kmin);
  idefix_for("ComputeSolenoidalForcing", 1, data->np_tot[KDIR] - 1, 1, data->np_tot[JDIR] - 1, 1, data->np_tot[IDIR] - 1,
              KOKKOS_LAMBDA (int k, int j, int i) {
                #if GEOMETRY == CARTESIAN
                      solenoidalForcingTerm(IDIR,k,j,i) = D_EXPAND( ZERO_F                                               ,
                                                 + 1/dx2(j) * (pristineForcingTerm(KDIR,k,j+1,i) - pristineForcingTerm(KDIR,k,j-1,i) )  ,
                                                 - 1/dx3(k) * (pristineForcingTerm(JDIR,k+1,j,i) - pristineForcingTerm(JDIR,k-1,j,i) )  );
                  #if DIMENSIONS >= 2
                      solenoidalForcingTerm(JDIR,k,j,i) = D_EXPAND( - 1/dx1(i) * (pristineForcingTerm(KDIR,k,j,i+1) - pristineForcingTerm(KDIR,k,j,i-1) )  ,
                                                                                                      ,
                                                 + 1/dx3(k) * (pristineForcingTerm(IDIR,k+1,j,i) - pristineForcingTerm(IDIR,k-1,j,i) )  );
                  #endif
                  #if DIMENSIONS == 3
                      solenoidalForcingTerm(KDIR,k,j,i) = 1/dx1(i) * (pristineForcingTerm(JDIR,k,j,i+1) - pristineForcingTerm(JDIR,k,j,i-1) )
                                       - 1/dx2(j) * (pristineForcingTerm(IDIR,k,j+1,i) - pristineForcingTerm(IDIR,k,j-1,i) );
                  #endif
                #endif
                #if GEOMETRY == CYLINDRICAL
                      IDEFIX_ERROR("Not yet defined");
                #endif
                #if GEOMETRY == POLAR
                      solenoidalForcingTerm(IDIR,k,j,i) = D_EXPAND( ZERO_F                                                       ,
                                                 + 1/(x1(i)*dx2(j)) * (pristineForcingTerm(KDIR,k,j+1,i) - pristineForcingTerm(KDIR,k,j-1,i) ) ,
                                                 - 1/dx3(k) * (pristineForcingTerm(JDIR,k+1,j,i) - pristineForcingTerm(JDIR,k-1,j,i) )          );

                      solenoidalForcingTerm(JDIR,k,j,i) = D_EXPAND( - 1/dx1(i) * (pristineForcingTerm(KDIR,k,j,i+1) - pristineForcingTerm(KDIR,k,j,i-1) )  ,
                                                                                                      ,
                                                 + 1/dx3(k) * (pristineForcingTerm(IDIR,k+1,j,i) - pristineForcingTerm(IDIR,k-1,j,i) )  );

                  #if DIMENSIONS == 3
                      solenoidalForcingTerm(KDIR,k,j,i) = 1/(x1(i)*dx1(i)) * (x1(i+1)*pristineForcingTerm(JDIR,k,j,i+1) - x1(i-1)*pristineForcingTerm(JDIR,k,j,i-1) )
                                       - 1/(x1(i)*dx2(j)) * (pristineForcingTerm(IDIR,k,j+1,i) - pristineForcingTerm(IDIR,k,j-1,i) );
                  #endif
                #endif
                #if GEOMETRY == SPHERICAL
                      real Ax2 = fabs(sinx2(j));
//                      // Regularisation along the axis
//                      if(FABS(Ax2)<1e-12) Ax2 = ONE_F;
                      solenoidalForcingTerm(IDIR,k,j,i) = D_EXPAND( ZERO_F                                                ,
                                                 + 1/(x1m(i)*Ax2) * (sinx2(j+1)*pristineForcingTerm(KDIR,k,j+1,i)
                                                 - sinx2(j-1)*pristineForcingTerm(KDIR,k,j-1,i) )                       ,
                                                 - 1/(x1(i)*Ax2*dx3(k)) * (pristineForcingTerm(JDIR,k+1,j,i)
                                                 - pristineForcingTerm(JDIR,k-1,j,i) )                                   );

                      solenoidalForcingTerm(JDIR,k,j,i) = D_EXPAND( - 1/(x1(i)*dx1(i)) * (x1(i+1)*pristineForcingTerm(KDIR,k,j,i+1)
                                                 - x1(i-1)*pristineForcingTerm(KDIR,k,j,i-1) )                           ,
                                                                                                      ,
                                                 + 1/(x1(i)*Ax2*dx3(k)) * (pristineForcingTerm(IDIR,k+1,j,i)
                                                 - pristineForcingTerm(IDIR,k-1,j,i) )                                  );

                  #if DIMENSIONS == 3
                      solenoidalForcingTerm(KDIR,k,j,i) = 1/(x1(i)*dx1(i)) * (x1(i+1)*pristineForcingTerm(JDIR,k,j,i+1) - x1(i-1)*pristineForcingTerm(JDIR,k,j,i-1) )
                                      - 1/(x1(i)*dx2(j)) * (pristineForcingTerm(IDIR,k,j+1,i) - pristineForcingTerm(IDIR,k,j-1,i) );
                  #endif
                #endif

                solenoidalForcingTerm(IDIR,k,j,i) /= 2.*kmean;
                solenoidalForcingTerm(JDIR,k,j,i) /= 2.*kmean;
                solenoidalForcingTerm(KDIR,k,j,i) /= 2.*kmean;
  });

  idfx::popRegion();
}

//// This function compute the compressive part of the forcing field
//void Forcing::ComputeCompressiveForcing(real dt) {
//  idfx::pushRegion("Forcing::ComputeCompressiveForcing");
//  
//  ResetForcingTerms();
//
//  IdefixArray2D<real> kiso = (this->forcingType == iso3D) ? this->k3Diso : this->k2Diso;
//  IdefixArray4D<Kokkos::complex<real>> forcingModesIdir = this->forcingModesIdir;
//  IdefixArray4D<Kokkos::complex<real>> forcingModesJdir = this->forcingModesJdir;
//  IdefixArray4D<Kokkos::complex<real>> forcingModesKdir = this->forcingModesKdir;
//  IdefixArray4D<real> compressiveForcingTerm = this->compressiveForcingTerm;
//  IdefixArray2D<Kokkos::complex<real>> ouValues = this->oUprocesses.ouValues;
//  int nForcingModes = this->nForcingModes;
//  Kokkos::complex unit_j(0.,1.);
//  idefix_for("ComputeCompressiveForcing", 0, data->np_tot[KDIR], 0, data->np_tot[JDIR], 0, data->np_tot[IDIR],
//              KOKKOS_LAMBDA (int k, int j, int i) {
//                  Kokkos::complex<real> forcingX1 = 0.;
//                  Kokkos::complex<real> forcingX2 = 0.;
//                  Kokkos::complex<real> forcingX3 = 0.;
//                  real kx1, kx2, kx3;
//                  for (int l=0; l<nForcingModes; l++) {
//                    kx1 = kiso(l,0);
//                    kx2 = kiso(l,1);
//                    kx3 = kiso(l,2);
//                    kk = kx1*kx1 + kx2*kx2 + kx3*kx3;
//                    //Here we are computing unit k
//                    kx1 /= kk;
//                    kx2 /= kk;
//                    kx3 /= kk;
//                    //Here we are coding, for each mode k, a_k . ik/k
//                    forcingX1 += unit_j*kx1*ouValues(l,IDIR)*forcingModesIdir(l,k,j,i);
//                    #if COMPONENTS >= 2
//                      forcingX2 += unit_j*(ouValues(l,KDIR)*kx1 - ouValues(l,IDIR)*kx3)*forcingModesJdir(l,k,j,i);
//                      #if COMPONENTS == 3
//                        forcingX3 += unit_j*(ouValues(l,IDIR)*kx2 - ouValues(l,JDIR)*kx1)*forcingModesKdir(l,k,j,i);
//                      #endif //COMPONENTS == 3
//                    #endif //COMPONENTS >= 2
//                  }
//                  compressiveForcingTerm(IDIR,k,j,i) += forcingX1.real();
//                  compressiveForcingTerm(JDIR,k,j,i) += forcingX2.real();
//                  compressiveForcingTerm(KDIR,k,j,i) += forcingX3.real();
//  });
//
//  idfx::popRegion();
//}

// Fill the forcing term with zeros
void Forcing::ResetForcingTerms() {
  idfx::pushRegion("Forcing::ResetForcingTerms");
  IdefixArray4D<real> forcingTerm = this->forcingTerm;
  IdefixArray4D<real> pristineForcingTerm = this->pristineForcingTerm;
  IdefixArray4D<real> solenoidalForcingTerm = this->solenoidalForcingTerm;
//  IdefixArray4D<real> compressiveForcingTerm = this->compressiveForcingTerm;
  idefix_for("Forcing::ResetForcingTerms",
              0, data->np_tot[KDIR],
              0, data->np_tot[JDIR],
              0, data->np_tot[IDIR],
              KOKKOS_LAMBDA(int k, int j, int i) {
                forcingTerm(IDIR,k,j,i) = ZERO_F;
                forcingTerm(JDIR,k,j,i) = ZERO_F;
                forcingTerm(KDIR,k,j,i) = ZERO_F;
                pristineForcingTerm(IDIR,k,j,i) = ZERO_F;
                pristineForcingTerm(JDIR,k,j,i) = ZERO_F;
                pristineForcingTerm(KDIR,k,j,i) = ZERO_F;
                solenoidalForcingTerm(IDIR,k,j,i) = ZERO_F;
                solenoidalForcingTerm(JDIR,k,j,i) = ZERO_F;
                solenoidalForcingTerm(KDIR,k,j,i) = ZERO_F;
//                compressiveForcingTerm(IDIR,k,j,i) = ZERO_F;
//                compressiveForcingTerm(JDIR,k,j,i) = ZERO_F;
//                compressiveForcingTerm(KDIR,k,j,i) = ZERO_F;
              });
  idfx::popRegion();
}

void Forcing::ComputeAverageSoundSpeed() {
  idfx::pushRegion("Forcing::ComputeAverageSoundSpeed");
  this->cs = 1.;

  DataBlockHost d(*this->data);
  d.SyncFromDevice();
  IdefixHostArray1D<real> x1=d.x[IDIR];
  IdefixHostArray1D<real> x2=d.x[JDIR];
  IdefixHostArray1D<real> x3=d.x[KDIR];
  IdefixHostArray4D<real> Vc=d.Vc;

  #if GEOMETRY == CARTESIAN || GEOMETRY == SPHERICAL
    real currentCs = ZERO_F;
    real currentVol = ZERO_F;
    #if GEOMETRY == CARTESIAN
      for(int i = d.beg[IDIR]; i < d.end[IDIR] ; i++) {
        for(int j = d.beg[JDIR]; j < d.end[JDIR] ; j++) {
          for(int k = d.beg[KDIR]; k < d.end[KDIR] ; k++) {
                currentCs += sqrt(d.Vc(PRS,k,j,i)/d.Vc(RHO,k,j,i));
                currentVol += d.dV(k,j,i);
           }
         }
       }
    #endif //GEOMETRY == CARTESIAN
    #if GEOMETRY == SPHERICAL
//      real r2, sint;
      for(int i = d.beg[IDIR]; i < d.end[IDIR] ; i++) {
//        r2 = x1(i)*x1(i);
        for(int j = d.beg[JDIR]; j < d.end[JDIR] ; j++) {
//          sint = SIN(d.x[JDIR](j));
          for(int k = d.beg[KDIR]; k < d.end[KDIR] ; k++) {
//                currentCs += r2*sint*sqrt(d.Vc(PRS,k,j,i)/d.Vc(RHO,k,j,i));
                currentCs += d.dV(k,j,i)*sqrt(d.Vc(PRS,k,j,i)/d.Vc(RHO,k,j,i));
                currentVol += d.dV(k,j,i);
          }
        }
      }
    #endif //GEOMETRY == SPHERICAL
    #ifdef WITH_MPI
      MPI_Allreduce(MPI_IN_PLACE, &currentCs, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &currentVol, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    #endif
    this->cs = currentCs/currentVol;
  #endif //GEOMETRY == CARTESIAN || GEOMETRY == SPHERICAL
  idfx::popRegion();
}

