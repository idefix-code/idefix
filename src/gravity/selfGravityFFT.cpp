// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifdef WITH_FFT

#include <cmath>
#include <memory>
#include <string>
#include <utility>
#include <vector>
#include <KokkosFFT.hpp>

#include "idefix.hpp"
#include "selfGravityFFT.hpp"
#include "dataBlock.hpp"
#include "fluid.hpp"
#include "vector.hpp"
#include "grid.hpp"

void SelfGravityFFT::Init(Input &input, DataBlock *datain) {
  idfx::pushRegion("SelfGravityFFT::Init");
  this->data = datain;
  Grid *grid = data->mygrid;
  CheckCompatibility();

  // FFT path currently requires periodic BC in all active dimensions
  for (int dir = 0; dir < DIMENSIONS; dir++) {
    if(grid->lbound[dir] != periodic || grid->rbound[dir] != periodic) {
      IDEFIX_ERROR("[SelfGravityFFT]: FFT solver requires periodic boundary "
                    "conditions in all active dimensions.");
    }
  }

  // storage with laplacian local array
  rho = IdefixArray3D<real>("Density", data->np_int[KDIR],
                                       data->np_tot[JDIR],
                                       data->np_tot[IDIR]);
  phi = IdefixArray3D<real>("Potential", data->np_tot[KDIR],
                                         data->np_tot[JDIR],
                                         data->np_tot[IDIR]);

  npr_glob = {grid->np_int[IDIR], grid->np_int[JDIR], grid->np_int[KDIR]};
  npr = {grid->np_int[IDIR], grid->np_int[JDIR], grid->np_int[KDIR]/idfx::psize};
  npf = {grid->np_int[IDIR]/2+1, grid->np_int[JDIR], grid->np_int[KDIR]/idfx::psize};
  npf_t = {grid->np_int[IDIR]/2+1, grid->np_int[KDIR], grid->np_int[JDIR]/idfx::psize};

  for(int dir = 0; dir < 3; dir++) {
    real d = (grid->xend[dir]-grid->xbeg[dir])/(2.0*M_PI*static_cast<real>(npr_glob[dir]));
    kx_glob[dir] = KokkosFFT::fftfreq(Device(), npr_glob[dir], d);
    int p = idfx::prank;
    if(dir != JDIR) kx[dir] = kx_glob[dir];
    else            kx[dir] = Kokkos::subview(kx_glob[dir],
                                  std::pair<int,int>(p * npf_t[KDIR], (p+1) * npf_t[KDIR]));
  }

  #ifdef WITH_MPI
    rhoF = IdefixArray3D<complex>("rhoHatFFT", npf_t[KDIR], npf_t[JDIR], npf_t[IDIR]);
    phiF = IdefixArray3D<complex>("phiHatFFT", npf_t[KDIR], npf_t[JDIR], npf_t[IDIR]);
  #else
    rhoF = IdefixArray3D<complex>("rhoHatFFT", npf[KDIR], npf[JDIR], npf[IDIR]);
    phiF = IdefixArray3D<complex>("phiHatFFT", npf[KDIR], npf[JDIR], npf[IDIR]);
  #endif

  std::array<int,3> nfft_real = {npr_glob[KDIR], npr_glob[JDIR], npr_glob[IDIR]};
  std::array<int,3> nfft_complex = {npr_glob[KDIR], npr_glob[JDIR], npr_glob[IDIR]/2+1};

  this->fft = std::make_unique<FFT>(nfft_real, nfft_complex);

  #ifdef WITH_MPI

    int ntarget = 0;
    std::vector<int> mapVars;
    mapVars.push_back(ntarget);

    this->mpi.Init(data->mygrid, mapVars, data->nghost,
                   data->np_int, data->lbound, data->rbound, false);
  #endif
  idfx::popRegion();
}

void SelfGravityFFT::CheckCompatibility() {
    #if GEOMETRY != CARTESIAN
    IDEFIX_ERROR("SelfGravityFFT supports CARTESIAN geometry only.");
    #endif
    Grid *grid = data->mygrid;
    if(grid->np_int[IDIR] != data->np_int[IDIR] || grid->np_int[JDIR] != data->np_int[JDIR]) {
      IDEFIX_ERROR("SelfGravityFFT requires no domain decomposition in the X1 and X2 directions.");
    }
    if(grid->np_int[KDIR] % idfx::psize != 0 || grid->np_int[JDIR] % idfx::psize != 0) {
      IDEFIX_ERROR("SelfGravityFFT requires that the number of grid points in the X3 and "
                   "X2 directions are divisible by the number of MPI processes.");
    }
    if(grid->np_int[IDIR] % 2 != 0) {
      IDEFIX_ERROR("SelfGravityFFT requires an even number of grid points in the X1 direction.");
    }
}

void SelfGravityFFT::ShowConfig() {
  idfx::cout << "SelfGravity: Using FFT Poisson solver (periodic Cartesian)." << std::endl;
  if(skipSelfGravity>1) {
    idfx::cout << "SelfGravity: self-gravity updated every " << skipSelfGravity
               << " cycles." << std::endl;
  }
}

void SelfGravityFFT::EnforcePeriodic(int dir, BoundarySide side, IdefixArray3D<real> &array) {
  idfx::pushRegion("Laplacian::EnforceBoundary");

  IdefixArray3D<real> localVar = array;

  // Number of active cells
  const int nxi = data->np_int[IDIR];
  const int nxj = data->np_int[JDIR];
  const int nxk = data->np_int[KDIR];

  // Number of ghost cells
  const int ighost = data->nghost[IDIR];
  const int jghost = data->nghost[JDIR];
  const int kghost = data->nghost[KDIR];

  // Boundaries of the loop
  const int ibeg = (dir == IDIR) ? side*(ighost+nxi) : 0;
  const int iend = (dir == IDIR) ? ighost + side*(ighost+nxi) : data->np_tot[IDIR];
  const int jbeg = (dir == JDIR) ? side*(jghost+nxj) : 0;
  const int jend = (dir == JDIR) ? jghost + side*(jghost+nxj) : data->np_tot[JDIR];
  const int kbeg = (dir == KDIR) ? side*(kghost+nxk) : 0;
  const int kend = (dir == KDIR) ? kghost + side*(kghost+nxk) : data->np_tot[KDIR];

  // Periodicity already enforced by MPI calls
  if(data->mygrid->nproc[dir] == 1) {
    idefix_for("BoundaryPeriodic", kbeg, kend, jbeg, jend, ibeg, iend,
        KOKKOS_LAMBDA (int k, int j, int i) {
            int iref, jref, kref;
            // This hack takes care of cases where we have more ghost zones than active zones
            if(dir==IDIR)
            iref = ighost + (i+ighost*(nxi-1))%nxi;
            else
            iref = i;
            if(dir==JDIR)
            jref = jghost + (j+jghost*(nxj-1))%nxj;
            else
            jref = j;
            if(dir==KDIR)
            kref = kghost + (k+kghost*(nxk-1))%nxk;
            else
            kref = k;

            localVar(k,j,i) = localVar(kref,jref,iref);
    });
  }
  idfx::popRegion();
}

void SelfGravityFFT::SetBoundaries(IdefixArray3D<real> &arr) {
  idfx::pushRegion("SelfGravityFFT::SetBoundaries");

  #ifdef WITH_MPI
  IdefixArray4D<real> arr4D = IdefixArray4D<real> (arr.data(), 1, data->np_tot[KDIR],
                                                    data->np_tot[JDIR],
                                                    data->np_tot[IDIR]);
  #endif

  for(int dir = 0 ; dir < DIMENSIONS ; dir++) {
    // MPI Exchange data when needed
    #ifdef WITH_MPI
    if(data->mygrid->nproc[dir]>1) {
      switch(dir) {
        case 0:
          this->mpi.ExchangeX1(arr4D);
          break;
        case 1:
          this->mpi.ExchangeX2(arr4D);
          break;
        case 2:
          this->mpi.ExchangeX3(arr4D);
          break;
      }
    }
    #endif

    EnforcePeriodic(dir, left, arr);
    EnforcePeriodic(dir, right, arr);
  }

  idfx::popRegion();
}

void SelfGravityFFT::SolvePoisson() {
  idfx::pushRegion("SelfGravityFFT::SolvePoisson");
  Kokkos::Timer timer;
  elapsedTime -= timer.seconds();

  // Make a view omitting the ghost cells
  auto rhoReal = Kokkos::subview(data->hydro->Vc,
                                RHO,
                                std::pair<int,int>(data->beg[KDIR],data->end[KDIR]),
                                std::pair<int,int>(data->beg[JDIR],data->end[JDIR]),
                                std::pair<int,int>(data->beg[IDIR],data->end[IDIR]));

  auto rhoF = this->rhoF;
  auto phiF = this->phiF;
  // Forward transform of the density field
  fft->R2C(rhoReal, rhoF, false);

  // Inverse poisson in Fourier space
  auto kx1 = kx[IDIR];
  auto kx2 = kx[JDIR];
  auto kx3 = kx[KDIR];
  #ifdef WITH_MPI
    // Work with transposed arrays
    idefix_for("PoissonFFT", 0, npf_t[KDIR], 0, npf_t[JDIR], 0, npf_t[IDIR],
        KOKKOS_LAMBDA(int j, int k, int i) {
        const real k2 = kx1(i)*kx1(i) + kx2(j)*kx2(j) + kx3(k)*kx3(k);
        real inv_k2 = (k2 > 0.0) ? -1.0/k2 : 0.0;
        phiF(j,k,i) = rhoF(j,k,i) * inv_k2;
        }
    );
  #else
    idefix_for("PoissonFFT", 0, npf[KDIR], 0, npf[JDIR], 0, npf[IDIR],
        KOKKOS_LAMBDA(int k, int j, int i) {
        const real k2 = kx1(i)*kx1(i) + kx2(j)*kx2(j) + kx3(k)*kx3(k);
        real inv_k2 = (k2 > 0.0) ? -1.0/k2 : 0.0;
        phiF(k,j,i) = rhoF(k,j,i) * inv_k2;
        }
    );
  #endif

  // make a ghost-cell free view of the potential field to apply the inverse transform
  auto phiReal = Kokkos::subview(phi,
                                std::pair<int,int>(data->beg[KDIR],data->end[KDIR]),
                                std::pair<int,int>(data->beg[JDIR],data->end[JDIR]),
                                std::pair<int,int>(data->beg[IDIR],data->end[IDIR]));

  fft->C2R(phiF, phiReal, false);
  // Need to apply the boundary conditions to the potential field to fill the ghost cells,
  SetBoundaries(phi);
  elapsedTime += timer.seconds();

  idfx::popRegion();
}

void SelfGravityFFT::AddSelfGravityPotential(IdefixArray3D<real> &phiP) {
  IdefixArray3D<real> localPot = phiP;
  IdefixArray3D<real> pot = this->phi;
  real gravCst = this->data->gravity->gravCst;


  idefix_for("AddSelfGravityPotentialFFT", 0, data->np_tot[KDIR],
                                          0, data->np_tot[JDIR],
                                          0, data->np_tot[IDIR],
    KOKKOS_LAMBDA (int k, int j, int i) {
      localPot(k,j,i) += 4.*M_PI*gravCst * pot(k, j, i);
    });
}

void SelfGravityFFT::EnrollUserDefBoundary(Laplacian::UserDefBoundaryFunc myFunc) {
  (void) myFunc;
  IDEFIX_ERROR("SelfGravityFFT only supports periodic boundaries; "
               "userdef boundaries are not supported.");
}

#endif // WITH_FFT
