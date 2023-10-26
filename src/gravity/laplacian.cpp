// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************



#include <string>
#include <vector>
#include "laplacian.hpp"
#include "selfGravity.hpp"
#include "dataBlock.hpp"


Laplacian::Laplacian(DataBlock *datain, std::array<LaplacianBoundaryType,3> leftBound,
                                        std::array<LaplacianBoundaryType,3> rightBound,
                                        bool havePreconditionnerIn) {
  idfx::pushRegion("Laplacian::Laplacian");
  // Save the parents data objects
  this->data = datain;
  this->havePreconditioner = havePreconditionnerIn;

  // init the grid elements using the parent datablock
  this->np_tot = data->np_tot;
  this->np_int = data->np_int;
  this->nghost = data->nghost;
  this->beg = data->beg;
  this->end = data->end;
  this->x = data->x;
  this->dx = data->dx;
  this->sinx2 = data->sinx2;
  this->dV = data->dV;
  this->A = data->A;
  this->loffset = {0,0,0};
  this->roffset = {0,0,0};

  this->lbound = leftBound;
  this->rbound = rightBound;

  isPeriodic = true;
  for(int dir = 0 ; dir < 3 ; dir++) {
    if(lbound[dir] != LaplacianBoundaryType::periodic) isPeriodic = false;
    if(rbound[dir] != LaplacianBoundaryType::periodic) isPeriodic = false;
  }

  #ifdef WITH_MPI
    if(lbound[IDIR] == origin) {
      // create communicator for spherical radius
      int remainDims[3] = {false, true, true};
      MPI_SAFE_CALL(MPI_Cart_sub(data->mygrid->CartComm, remainDims, &originComm));
    }

  // Update internal boundaries in case of domain decomposition
    for(int dir = 0 ; dir < DIMENSIONS ; dir++) {
      if(data->mygrid->nproc[dir]>1) {
        if(this->data->lbound[dir]==internal ) {
          this->lbound[dir] = Laplacian::internalgrav;
        }
        if(this->data->rbound[dir]==internal) {
          this->rbound[dir] = Laplacian::internalgrav;
        }
      }
    }
  #endif

  #if GEOMETRY == SPHERICAL
    if ((this->rbound[JDIR]==axis) || (this->lbound[JDIR]==axis)) {
      // Check wether the x3 spherical axis is full two pi
      if(fabs((data->mygrid->xend[KDIR] - data->mygrid->xbeg[KDIR] -2.0*M_PI)) < 1e-10) {
        this->isTwoPi = true;
      }

    #ifdef WITH_MPI
      // Check that there is no domain decomposition in phi
      if(data->mygrid->nproc[KDIR]>1) {
        IDEFIX_ERROR("Laplacian:: Axis boundaries are not compatible with "
                     "MPI domain decomposition in X3");
      }
    #endif
    }
  #endif

  if(this->lbound[IDIR] == origin) {
    InitInternalGrid();
  }

  // Init preconditionner if needed
  if(havePreconditioner) {
    InitPreconditionner();
  }
  // Initialise the Laplacian coefficients
  PreComputeLaplacian();

  // Init MPI stack when needed
  #ifdef WITH_MPI
    this->arr4D = IdefixArray4D<real> ("WorkingArrayMpi", 1, this->np_tot[KDIR],
                                                            this->np_tot[JDIR],
                                                            this->np_tot[IDIR]);

    int ntarget = 0;
    std::vector<int> mapVars;
    mapVars.push_back(ntarget);

    this->mpi.Init(data->mygrid, mapVars, this->nghost.data(), this->np_int.data());
  #endif

  idfx::popRegion();
}

void Laplacian::InitInternalGrid() {
  idfx::pushRegion("Laplacian::InitInternalGrid");
  // Extend the grid so that the inner radius will be 1/10 of the initial inner radius
  GridHost gh(*data->mygrid);
  gh.SyncFromDevice();

  real dx0 = gh.dx[IDIR](0);
  real rin = gh.xbeg[IDIR]/10.0;

  loffset[IDIR] = (gh.xl[IDIR](0) - rin) / dx0;

  this->np_tot[IDIR] += loffset[IDIR];
  this->np_int[IDIR] += loffset[IDIR];
  this->end[IDIR] += loffset[IDIR];

  // Allocate required arrays
  this->x[IDIR] = IdefixArray1D<real>("SG_x1", this->np_tot[IDIR]);
  this->dx[IDIR] = IdefixArray1D<real>("SG_dx1", this->np_tot[IDIR]);
  this->dV = IdefixArray3D<real>("SG_dV", this->np_tot[KDIR],
                                          this->np_tot[JDIR],
                                          this->np_tot[IDIR]);
  for(int dir = 0 ; dir < 3 ; dir ++) {
    this->A[dir] = IdefixArray3D<real>("SG_A", this->np_tot[KDIR]+KOFFSET,
                                          this->np_tot[JDIR]+JOFFSET,
                                          this->np_tot[IDIR]+IOFFSET);
  }
  // xl and xr are used only temporarily in this subroutine.
  auto x1l = IdefixArray1D<real>("SG_x1l", this->np_tot[IDIR]);
  auto x1r = IdefixArray1D<real>("SG_x1r", this->np_tot[IDIR]);
  // incidentally, sinx2 is not affected by a grid extension in x1, so no need to
  // allocate a new array


  // Fill selfgravity with the existing datablock grid
  int ioffset = loffset[IDIR];

  auto x1in = data->x[IDIR];
  auto x1lin = data->xl[IDIR];
  auto x1rin = data->xr[IDIR];
  auto dx1in = data->dx[IDIR];
  auto x1 = this->x[IDIR];
  auto dx1 = this->dx[IDIR];
  idefix_for("InternalGridCopy",0, data->np_tot[IDIR],
     KOKKOS_LAMBDA(int i) {
       x1(i+ioffset) = x1in(i);
       dx1(i+ioffset) = dx1in(i);
       x1l(i+ioffset) = x1lin(i);
       x1r(i+ioffset) = x1rin(i);
     });

  // extend with a uniform grid and spacing equal to first of the datablock
  idefix_for("InternalGridFill",0, ioffset,
     KOKKOS_LAMBDA(int i) {
       const real dx0 = dx1(ioffset);
       dx1(i) = dx0;
       x1(i) = x1(ioffset) + (i-ioffset)*dx0;
       x1l(i) = x1l(ioffset) + (i-ioffset)*dx0;
       x1r(i) = x1l(ioffset) + (i-ioffset+1)*dx0;
     });

  // Check that all is well
  IdefixArray1D<real>::HostMirror xH = Kokkos::create_mirror_view(x1l);
  Kokkos::deep_copy(xH, x1l);

  if(xH(0)<0.0) {
    IDEFIX_ERROR("SelfGravity: Your grid extension goes beyond the origin!");
  }

  IdefixArray3D<real> dV  = this->dV;
  //IdefixArray1D<real> dx1 = this->dx[IDIR];
  IdefixArray1D<real> dx2 = this->dx[JDIR];
  IdefixArray1D<real> dx3 = this->dx[KDIR];
  //IdefixArray1D<real> x1  = this->x[IDIR];
  IdefixArray1D<real> x2  = this->x[JDIR];
  IdefixArray1D<real> x3  = this->x[KDIR];
  IdefixArray1D<real> x1p = x1r;
  IdefixArray1D<real> x2p = data->xr[JDIR];
  IdefixArray1D<real> x3p = data->xr[KDIR];
  IdefixArray1D<real> x1m = x1l;
  IdefixArray1D<real> x2m = data->xl[JDIR];
  IdefixArray1D<real> x3m = data->xl[KDIR];

  // Fill the volume and area arrays (assume spherical geometry)
  idefix_for("Volumes",0,this->np_tot[KDIR],0,this->np_tot[JDIR],0,this->np_tot[IDIR],
    KOKKOS_LAMBDA (int k, int j, int i) {
      real dVr = FABS(x1p(i)*x1p(i)*x1p(i) - x1m(i)*x1m(i)*x1m(i))/3.0;
      real dmu = FABS(cos(x2m(j)) - cos(x2p(j)));
      dV(k,j,i) = D_EXPAND( dVr, *dmu, *dx3(k));
    }
  );

  // Compute Areas
  IdefixArray3D<real> Ax1 = this->A[IDIR];
  IdefixArray3D<real> Ax2 = this->A[JDIR];
  IdefixArray3D<real> Ax3 = this->A[KDIR];

  // X1 direction
  int end = this->np_tot[IDIR];
  idefix_for("AreaX1",0,this->np_tot[KDIR],0,this->np_tot[JDIR],0,this->np_tot[IDIR]+IOFFSET,
    KOKKOS_LAMBDA (int k, int j, int i) {
      real dmu = FABS(cos(x2m(j)) - cos(x2p(j)));
      if(i == end) {
        Ax1(k,j,i) = D_EXPAND(x1p(i-1)*x1p(i-1), *dmu, *dx3(k));
      } else {
        Ax1(k,j,i) = D_EXPAND(x1m(i)*x1m(i), *dmu, *dx3(k)); // r^2*dmu*dphi
      }
    });

  // X2 direction
  end = this->np_tot[JDIR];
  idefix_for("AreaX2",0,this->np_tot[KDIR],0,this->np_tot[JDIR]+JOFFSET,0,this->np_tot[IDIR],
    KOKKOS_LAMBDA (int k, int j, int i) {
      if (j == end) {
        Ax2(k,j,i) = D_EXPAND(x1(i)*dx1(i), *FABS(sin(x2p(j-1))), *dx3(k));
      } else {
        Ax2(k,j,i) = D_EXPAND(x1(i)*dx1(i), *FABS(sin(x2m(j))), *dx3(k)); // = r*dr*sin(thp)*dphi
      }
    });

  // X3 direction
  end = this->np_tot[KDIR];
  idefix_for("AreaX3",0,this->np_tot[KDIR]+KOFFSET,0,this->np_tot[JDIR],0,this->np_tot[IDIR],
    KOKKOS_LAMBDA (int k, int j, int i) {
      Ax3(k,j,i) = D_EXPAND(x1(i)*dx1(i), *dx2(j), *ONE_F);   // = r*dr*dth
    }
  );

  idfx::popRegion();
}

void Laplacian::InitPreconditionner() {
  idfx::pushRegion("Laplacian::InitPreconditioner");
  IdefixArray3D<real> P = IdefixArray3D<real> ("Preconditionner", this->np_tot[KDIR],
                                                                  this->np_tot[JDIR],
                                                                  this->np_tot[IDIR]);

  int ibeg, iend, jbeg, jend, kbeg, kend;
  ibeg = this->beg[IDIR];
  iend = this->end[IDIR];
  jbeg = this->beg[JDIR];
  jend = this->end[JDIR];
  kbeg = this->beg[KDIR];
  kend = this->end[KDIR];

  // Set array to zero
  idefix_for("ResetPrecond", kbeg, kend, jbeg, jend, ibeg, iend,
    KOKKOS_LAMBDA (int k, int j, int i) {
      P(k,j,i) = 0;
    });

  IdefixArray3D<real> dV = this->dV;

  #if (GEOMETRY==CYLINDRICAL) || (GEOMETRY==POLAR)
  IdefixArray1D<real> r = this->x[IDIR];
  #elif GEOMETRY==SPHERICAL
  IdefixArray1D<real> r = this->x[IDIR];
  IdefixArray1D<real> sinth = this->sinx2;
  #endif

  // Loop on dimensions to get the diagonal elements
  for(int dir = 0 ; dir < DIMENSIONS ; dir++) {
    IdefixArray1D<real> dx = this->dx[dir];
    IdefixArray3D<real> Ax = this->A[dir];
    int ioffset = (dir == IDIR) ? 1 : 0;
    int joffset = (dir == JDIR) ? 1 : 0;
    int koffset = (dir == KDIR) ? 1 : 0;

    idefix_for("InitPrecond", kbeg, kend, jbeg, jend, ibeg, iend,
    KOKKOS_LAMBDA (int k, int j, int i) {
      int index = i*ioffset + j*joffset + k*koffset;
      real dx0 = dx(index);
      real dxm = dx(index-1);
      real dxp = dx(index+1);

      real Ap = Ax(k+koffset, j+joffset, i+ioffset);
      real Am = Ax(k,j,i);

      // Compute the diagonal elements coming from the laplacian
      real d = -(2.0*Ap/(dx0+dxp) + 2.0*Am/(dx0+dxm)) / dV(k,j,i);

      // Handling curvature terms
      real h1, h2, h3;
      #if GEOMETRY==CARTESIAN
      h1=1.;
      h2=1.;
      h3=1.;
      #elif GEOMETRY==CYLINDRICAL
      h1=1.;
      h2=1.;
      h3=r(i); // Should never be used (cf. Idefix docu)
      #elif GEOMETRY==POLAR
      h1=1.;
      h2=r(i);
      h3=1.;
      #else
      h1=1.;
      h2=r(i);
      h3=r(i)*sinth(j);
      #endif

      // Full diagonal element (including curvature terms)
      d = d/(h1*ioffset+h2*joffset+h3*koffset);

      // We store in P the sum of all these elements
      P(k,j,i) += d;
    });
  }

  // Store the array for future use
  this->precond = P;

  idfx::popRegion();
}

void Laplacian::PreComputeLaplacian() {
  idfx::pushRegion("Laplacian::PreComputeLaplacian");
  // Precompute Laplacian Factor

  // Allocate Laplacian factors
  this->Lx1 = IdefixArray4D<real>("SelfGravity_Lx1",2,
                                                    this->np_tot[KDIR],
                                                    this->np_tot[JDIR],
                                                    this->np_tot[IDIR]);
  #if DIMENSIONS > 1
    this->Lx2 = IdefixArray4D<real>("SelfGravity_Lx2",2,
                                                      this->np_tot[KDIR],
                                                      this->np_tot[JDIR],
                                                      this->np_tot[IDIR]);

    #if DIMENSIONS > 2
      this->Lx3 = IdefixArray4D<real>("SelfGravity_Lx3",2,
                                                        this->np_tot[KDIR],
                                                        this->np_tot[JDIR],
                                                        this->np_tot[IDIR]);
    #endif
  #endif


  // Local copy of Laplacian arrays
  IdefixArray4D<real> Lx1 = this->Lx1; // X stride
  IdefixArray4D<real> Lx2 = this->Lx2; // Y stride
  IdefixArray4D<real> Lx3 = this->Lx3; // Z stride


  IdefixArray3D<real> P = this->precond;
  bool havePreconditioner = this->havePreconditioner;
  IdefixArray1D<real> dx1 = this->dx[IDIR];
  IdefixArray1D<real> dx2 = this->dx[JDIR];
  IdefixArray1D<real> dx3 = this->dx[KDIR];
  IdefixArray3D<real> Ax1 = this->A[IDIR];
  IdefixArray3D<real> Ax2 = this->A[JDIR];
  IdefixArray3D<real> Ax3 = this->A[KDIR];
  IdefixArray3D<real> dV = this->dV;
  #if GEOMETRY == POLAR
  IdefixArray1D<real> r = this->x[IDIR];
  #elif GEOMETRY == SPHERICAL
  IdefixArray1D<real> r = this->x[IDIR];
  IdefixArray1D<real> sinth = this->sinx2;
  #endif

  idefix_for("L_Factor",KOFFSET,this->np_tot[KDIR]-KOFFSET,
                        JOFFSET,this->np_tot[JDIR]-JOFFSET,
                        IOFFSET,this->np_tot[IDIR]-IOFFSET,
    KOKKOS_LAMBDA(int k, int j,int i) {
      [[maybe_unused]] real h1, h2, h3;
      #if GEOMETRY == CARTESIAN
      h1 = 1.;
      h2 = 1.;
      h3 = 1.;
      #elif GEOMETRY == POLAR
      h1 = 1.;
      h2 = r(i);
      h3 = 1.;
      #else
      h1 = 1.;
      h2 = r(i);
      h3 = r(i) * sinth(j);
      #endif

      // i-1 coefficient
      Lx1(0,k,j,i) = 2.0 * Ax1(k, j, i) / h1 / (dx1(i) + dx1(i-1)) / dV(k,j,i);
      // i+1 coefficient
      Lx1(1,k,j,i) = 2.0 * Ax1(k, j, i+1) / h1 / (dx1(i+1) + dx1(i)) / dV(k,j,i);
      if(havePreconditioner) {
        Lx1(0,k,j,i) /= P(k,j,i);
        Lx1(1,k,j,i) /= P(k,j,i);
      }
      #if DIMENSIONS > 1
        Lx2(0,k,j,i) = 2.0 * Ax2(k, j, i) / h2 / (dx2(j) + dx2(j-1)) / dV(k,j,i);
        Lx2(1,k,j,i) = 2.0 * Ax2(k, j+1, i) / h2 / (dx2(j+1) + dx2(j)) / dV(k,j,i);
        if(havePreconditioner) {
          Lx2(0,k,j,i) /= P(k,j,i);
          Lx2(1,k,j,i) /= P(k,j,i);
        }
        #if DIMENSIONS > 2
          Lx3(0,k,j,i) = 2.0 * Ax3(k, j, i) / h3 / (dx3(k) + dx3(k-1)) / dV(k,j,i);
          Lx3(1,k,j,i) = 2.0 * Ax3(k+1, j, i) / h3 / (dx3(k+1) + dx3(k)) / dV(k,j,i);
          if(havePreconditioner) {
            Lx3(0,k,j,i) /= P(k,j,i);
            Lx3(1,k,j,i) /= P(k,j,i);
          }
        #endif
      #endif
    });
  idfx::popRegion();
}


void Laplacian::operator()(IdefixArray3D<real> array, IdefixArray3D<real> laplacian) {
  idfx::pushRegion("Laplacian::ComputeLaplacian");

  int ibeg, iend, jbeg, jend, kbeg, kend;
  ibeg = this->beg[IDIR];
  iend = this->end[IDIR];
  jbeg = this->beg[JDIR];
  jend = this->end[JDIR];
  kbeg = this->beg[KDIR];
  kend = this->end[KDIR];
  IdefixArray4D<real> Lx1 = this->Lx1;
  #if DIMENSIONS > 1
    IdefixArray4D<real> Lx2 = this->Lx2;
    #if DIMENSIONS > 2
      IdefixArray4D<real> Lx3 = this->Lx3;
    #endif
  #endif

  // Handling boundaries before laplacian calculation
  this->SetBoundaries(array);

  idefix_for("FiniteDifference", kbeg, kend, jbeg, jend, ibeg, iend,
    KOKKOS_LAMBDA (int k, int j, int i) {
      real gc = 0;
      real Delta = 0;
      real Lm, Lr;
      #if DIMENSIONS > 2
      Lm = Lx3(0,k,j,i);
      Lr = Lx3(1,k,j,i);
      gc += Lm + Lr;
      Delta += array(k-1,j,i)*Lm + array(k+1,j,i)*Lr;
      #endif
      #if DIMENSIONS > 1
        Lm = Lx2(0,k,j,i);
        Lr = Lx2(1,k,j,i);
        gc += Lm + Lr;
        Delta += array(k,j-1,i)*Lm + array(k,j+1,i)*Lr;
      #endif
      Lm = Lx1(0,k,j,i);
      Lr = Lx1(1,k,j,i);
      gc += Lm + Lr;
      Delta += array(k,j,i-1)*Lm + array(k,j,i+1)*Lr;

      Delta -= gc * array(k,j,i);

      laplacian(k, j, i) = Delta;
    });

  idfx::popRegion();
}

void Laplacian::EnforceBoundary(int dir, BoundarySide side, LaplacianBoundaryType type,
                                  IdefixArray3D<real> &arr) {
  idfx::pushRegion("Laplacian::EnforceBoundary");

  IdefixArray3D<real> localVar = arr;

  // Number of active cells
  const int nxi = this->np_int[IDIR];
  const int nxj = this->np_int[JDIR];
  const int nxk = this->np_int[KDIR];

  // Number of ghost cells
  const int ighost = this->nghost[IDIR];
  const int jghost = this->nghost[JDIR];
  const int kghost = this->nghost[KDIR];

  // Boundaries of the loop
  const int ibeg = (dir == IDIR) ? side*(ighost+nxi) : 0;
  const int iend = (dir == IDIR) ? ighost + side*(ighost+nxi) : this->np_tot[IDIR];
  const int jbeg = (dir == JDIR) ? side*(jghost+nxj) : 0;
  const int jend = (dir == JDIR) ? jghost + side*(jghost+nxj) : this->np_tot[JDIR];
  const int kbeg = (dir == KDIR) ? side*(kghost+nxk) : 0;
  const int kend = (dir == KDIR) ? kghost + side*(kghost+nxk) : this->np_tot[KDIR];

  switch(type) {
    case internalgrav:
      // internal is used for MPI-enforced boundary conditions. Nothing to be done here.
      break;

    case periodic: {
      if(data->mygrid->nproc[dir] > 1) break; // Periodicity already enforced by MPI calls

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
      break;
    }

    case userdef: {
      if(this->haveUserDefBoundary) {
        // Warning: unlike hydro userdef boundary functions, the selfGravity
        // userdef boundary functions take an additional argument arr which
        // specifies the array for which boundaries are to be handled
        this->userDefBoundaryFunc(*data, dir, side, data->t, arr);
      } else {
        IDEFIX_ERROR("Laplacian:: No function enrolled to define your own boundary conditions");
      }
      break;
    }

    case nullpot: {
      idefix_for("BoundaryNullPot", kbeg, kend, jbeg, jend, ibeg, iend,
            KOKKOS_LAMBDA (int k, int j, int i) {
              localVar(k,j,i) = 0.0;
      });
      break;
    }

    case nullgrad: {
      idefix_for("BoundaryNullGrad", kbeg, kend, jbeg, jend, ibeg, iend,
            KOKKOS_LAMBDA (int k, int j, int i) {
              const int iref = (dir==IDIR) ? ighost + side*(nxi-1) : i;
              const int jref = (dir==JDIR) ? jghost + side*(nxj-1) : j;
              const int kref = (dir==KDIR) ? kghost + side*(nxk-1) : k;

              localVar(k,j,i) = localVar(kref,jref,iref);
      });
      break;
    }

    case axis: {
      // Handling specific loop boundaries
      int jref, offset;
      if(side == left) {
        jref = this->beg[JDIR];
        offset = -1;
      }
      if(side == right) {
        jref = this->end[JDIR]-1;
        offset = 1;
      }

      // NB: we assume no domain decomposition along phi here

      int np_int_k = this->np_int[KDIR];
      int nghost_k = this->nghost[KDIR];

      if(this->isTwoPi) {
        idefix_for("BoundaryAxis",kbeg,kend,jbeg,jend,ibeg,iend,
                KOKKOS_LAMBDA (int k, int j, int i) {
                  int kcomp = nghost_k + (( k - nghost_k + np_int_k/2) % np_int_k);
                  // Assuming sVc=1 for a scalar
                  localVar(k,j,i) = localVar(kcomp, 2*jref-j+offset,i);
                });
      } else { // not 2pi
        idefix_for("BoundaryAxis",kbeg,kend,jbeg,jend,ibeg,iend,
                KOKKOS_LAMBDA (int k, int j, int i) {
                  // kcomp = k by construction since we're doing a fraction of twopi
                  // Assuming sVc=1 for a scalar
                  localVar(k,j,i) = localVar(k, 2*jref-j+offset,i);
                });
      }
      break;
    }

    case origin: {
      // Assume the grid is extended inwards close to the origin. Hence the inner ghost
      // all have the same value.
      int iref = this->beg[IDIR];

      real psiIn = 0.0;

      idefix_reduce("meanPsiIn",
                    beg[KDIR], end[KDIR],
                    beg[JDIR], end[JDIR],
                    KOKKOS_LAMBDA(int k, int j, real &psi) {
                      psi += localVar(k,j,iref);
                    },Kokkos::Sum<real> (psiIn));

      #ifdef WITH_MPI
        MPI_Allreduce(MPI_IN_PLACE, &psiIn, 1, realMPI, MPI_SUM, originComm);
      #endif
      // Do a mean by dividing by the number of points
      psiIn = psiIn/(data->mygrid->np_int[JDIR]*data->mygrid->np_int[KDIR]);

      // put this in the ghost cells
      idefix_for("BoundaryOrigin",kbeg,kend,jbeg,jend,ibeg,iend,
                KOKKOS_LAMBDA (int k, int j, int i) {
                  localVar(k,j,i) = psiIn;
                });
      break;
    }

    default: {
      std::stringstream msg ("Laplacian:: Boundary condition type is not yet implemented");
      IDEFIX_ERROR(msg);
    }
  }

  idfx::popRegion();
}


void Laplacian::EnrollUserDefBoundary(UserDefBoundaryFunc myFunc) {
  this->userDefBoundaryFunc = myFunc;
  this->haveUserDefBoundary = true;
}

void Laplacian::SetBoundaries(IdefixArray3D<real> &arr) {
  idfx::pushRegion("Laplacian::SetBoundaries");

  #ifdef WITH_MPI
  this->arr4D = IdefixArray4D<real> (arr.data(), 1, this->np_tot[KDIR],
                                                    this->np_tot[JDIR],
                                                    this->np_tot[IDIR]);
  #endif

  for(int dir = 0 ; dir < DIMENSIONS ; dir++) {
    // MPI Exchange data when needed
    #ifdef WITH_MPI
    if(data->mygrid->nproc[dir]>1) {
      switch(dir) {
        case 0:
          this->mpi.ExchangeX1(this->arr4D);
          break;
        case 1:
          this->mpi.ExchangeX2(this->arr4D);
          break;
        case 2:
          this->mpi.ExchangeX3(this->arr4D);
          break;
      }
    }
    #endif

    EnforceBoundary(dir, left, this->lbound[dir], arr);
    EnforceBoundary(dir, right, this->rbound[dir], arr);
  }

  idfx::popRegion();
}

real Laplacian::ComputeCFL() {
  idfx::pushRegion("Laplacian::ComputeCFL");

  int ibeg = this->beg[IDIR];
  int iend = this->end[IDIR];
  int jbeg = this->beg[JDIR];
  int jend = this->end[JDIR];
  int kbeg = this->beg[KDIR];
  int kend = this->end[KDIR];

  #if DIMENSIONS == 1
  IdefixArray1D<real> dx1 = this->dx[IDIR];
  real dx1min2;

  idefix_reduce("GetMin1",
                kbeg, kend,
                jbeg, jend,
                ibeg, iend,
                KOKKOS_LAMBDA (int k, int j, int i, real &localMin) {
                  localMin = std::fmin(dx1(i) * dx1(i), localMin);
                },
                Kokkos::Min<real>(dx1min2));

    // Reduction on the whole grid
    #ifdef WITH_MPI
    MPI_Allreduce(MPI_IN_PLACE, &dx1min2, 1, realMPI, MPI_MIN, MPI_COMM_WORLD);
    #endif

  real dtmax = 1. / 2. * dx1min2;

  #elif DIMENSIONS == 2
  IdefixArray1D<real> dx1 = this->dx[IDIR];
  IdefixArray1D<real> dx2 = this->dx[JDIR];
  real dx1min2, dx2min2;
    #if (GEOMETRY == POLAR || GEOMETRY == SPHERICAL)
    IdefixArray1D<real> r = this->x[IDIR];
    #endif

  idefix_reduce("GetMin1",
                kbeg, kend,
                jbeg, jend,
                ibeg, iend,
                KOKKOS_LAMBDA (int k, int j, int i, real &localMin) {
                  localMin = std::fmin(dx1(i) * dx1(i), localMin);
                },
                Kokkos::Min<real>(dx1min2));

  idefix_reduce("GetMin2",
                kbeg, kend,
                jbeg, jend,
                ibeg, iend,
                KOKKOS_LAMBDA (int k, int j, int i, real &localMin) {
                  real dl = dx2(j);
                  #if (GEOMETRY == POLAR || GEOMETRY == SPHERICAL)
                  dl = dl * r(i);
                  #endif
                  localMin = std::fmin(dl * dl, localMin);
                },
                Kokkos::Min<real>(dx2min2));

    // Reduction on the whole grid
    #ifdef WITH_MPI
    MPI_Allreduce(MPI_IN_PLACE, &dx1min2, 1, realMPI, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &dx2min2, 1, realMPI, MPI_MIN, MPI_COMM_WORLD);
    #endif

  real dtmax = 1. / 2. * 1. / ( 1. / dx1min2 + 1. / dx2min2);

  #else
  IdefixArray1D<real> dx1 = this->dx[IDIR];
  IdefixArray1D<real> dx2 = this->dx[JDIR];
  IdefixArray1D<real> dx3 = this->dx[KDIR];
  real dx1min2, dx2min2, dx3min2;
    #if GEOMETRY == POLAR
    IdefixArray1D<real> r = x[IDIR];
    #elif GEOMETRY == SPHERICAL
    IdefixArray1D<real> r = x[IDIR];
    IdefixArray1D<real> sinth = sinx2;
    #endif

  idefix_reduce("GetMin1",
                kbeg, kend,
                jbeg, jend,
                ibeg, iend,
                KOKKOS_LAMBDA (int k, int j, int i, real &localMin) {
                  localMin = std::fmin(dx1(i) * dx1(i), localMin);
                },
                Kokkos::Min<real>(dx1min2));

  idefix_reduce("GetMin2",
                kbeg, kend,
                jbeg, jend,
                ibeg, iend,
                KOKKOS_LAMBDA (int k, int j, int i, real &localMin) {
                  real dl = dx2(j);
                  #if (GEOMETRY == POLAR || GEOMETRY == SPHERICAL)
                  dl = dl * r(i);
                  #endif
                  localMin = std::fmin(dl * dl, localMin);
                },
                Kokkos::Min<real>(dx2min2));

  idefix_reduce("GetMin3",  // Cylindrical not taken into account as it shouldn't be used in 3D
                kbeg, kend,
                jbeg, jend,
                ibeg, iend,
                KOKKOS_LAMBDA (int k, int j, int i, real &localMin) {
                  real dl = dx3(k);
                  #if GEOMETRY == SPHERICAL
                  dl = dl * r(i) * sinth(j);
                  #endif
                  localMin = std::fmin(dl * dl, localMin);
                },
                Kokkos::Min<real>(dx3min2));

    // Reduction on the whole grid
    #ifdef WITH_MPI
    MPI_Allreduce(MPI_IN_PLACE, &dx1min2, 1, realMPI, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &dx2min2, 1, realMPI, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &dx3min2, 1, realMPI, MPI_MIN, MPI_COMM_WORLD);
    #endif

  real dtmax = 1. / 2. * 1. / ( 1. / dx1min2 + 1. / dx2min2 + 1. / dx3min2);
  #endif

  idfx::popRegion();

  return(0.95 * dtmax); // Taking a percentage to avoid dt=dtmax leading to a breakup
}
