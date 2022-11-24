// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2022 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************



#include <string>
#include <vector>

#include "selfGravity.hpp"
#include "dataBlock.hpp"
#include "hydro.hpp"
#include "vector.hpp"
#include "bicgstab.hpp"
#include "cg.hpp"
#include "minres.hpp"
#include "jacobi.hpp"




SelfGravity::SelfGravity() {
  // Default constructor
}

void SelfGravity::Init(Input &input, DataBlock *datain) {
  idfx::pushRegion("SelfGravity::Init");

  // Save the parents data objects
  this->data = datain;

  // Initialise (default) solver parameters
  this->dt = 0.;
  this->isPeriodic = true;

  // Update targetError when provided
  real targetError = input.GetOrSet<real>("SelfGravity","targetError",0,1e-2);

  // Get maxiter when provided
  real maxiter = input.GetOrSet<int>("SelfGravity","maxIter",0,1000);

  // Get the number of skipped cycles when provided and check consistency
  this->skipSelfGravity = input.GetOrSet<int>("SelfGravity","skip",0,1);
  if(skipSelfGravity<1) {
    IDEFIX_ERROR("[SelfGravity]:skip should be a strictly positive integer");
  }

  // Get the gravity-related boundary conditions
  for(int dir = 0 ; dir < 3 ; dir++) {
    std::string label = std::string("boundary-X")+std::to_string(dir+1)+std::string("-beg");
    std::string boundary = input.Get<std::string>("SelfGravity",label,0);

    if(boundary.compare("nullpot") == 0) {
      this->lbound[dir] = nullpot;
      this->isPeriodic = false;
    } else if(boundary.compare("periodic") == 0) {
      this->lbound[dir] = periodic;
    } else if(boundary.compare("nullgrad") == 0) {
      this->lbound[dir] = nullgrad;
      this->isPeriodic = false;
    } else if(boundary.compare("internalgrav") == 0) {
      this->lbound[dir] = internalgrav;
      this->isPeriodic = false;
    } else if(boundary.compare("userdef") == 0) {
      this->lbound[dir] = userdef;
      this->isPeriodic = false;
    } else if(boundary.compare("axis") == 0) {
      this->lbound[dir] = axis;
      this->isPeriodic = false;
    } else if(boundary.compare("origin") == 0) {
      this->lbound[dir] = origin;
      this->isPeriodic = false;
      // origin only compatible with spherical & axis=IDIR
      #if GEOMETRY != SPHERICAL
        IDEFIX_ERROR("Origin boundary conditions are working in spherical coordinates");
      #endif
      if(dir != IDIR) {
        IDEFIX_ERROR("Origin boundary conditions are meaningful only on the X1 direction");
      }
      #ifdef WITH_MPI
        // create communicator for spherical radius
        int remainDims[3] = {false, true, true};
        MPI_SAFE_CALL(MPI_Cart_sub(data->mygrid->CartComm, remainDims, &originComm));
      #endif

    } else {
      std::stringstream msg;
      msg << "SelfGravity:: Unknown boundary type " << boundary;
      IDEFIX_ERROR(msg);
    }

    label = std::string("boundary-X")+std::to_string(dir+1)+std::string("-end");
    boundary = input.Get<std::string>("SelfGravity",label,0);
    if(boundary.compare("nullpot") == 0) {
      this->rbound[dir] = nullpot;
      this->isPeriodic = false;
    } else if(boundary.compare("periodic") == 0) {
      this->rbound[dir] = periodic;
    } else if(boundary.compare("nullgrad") == 0) {
      this->rbound[dir] = nullgrad;
      this->isPeriodic = false;
    } else if(boundary.compare("internalgrav") == 0) {
      this->rbound[dir] = internalgrav;
      this->isPeriodic = false;
    } else if(boundary.compare("userdef") == 0) {
      this->rbound[dir] = userdef;
      this->isPeriodic = false;
    } else if(boundary.compare("axis") == 0) {
      this->rbound[dir] = axis;
      this->isPeriodic = false;
    } else {
      std::stringstream msg;
      msg << "SelfGravity:: Unknown boundary type " << boundary;
      IDEFIX_ERROR(msg);
    }
  }

  // Update internal boundaries in case of domain decomposition
  #ifdef WITH_MPI
  for(int dir = 0 ; dir < DIMENSIONS ; dir++) {
    if(data->mygrid->nproc[dir]>1) {
      if(this->data->lbound[dir]==internal) {
        this->lbound[dir] = internalgrav;
      }
      if(this->data->rbound[dir]==internal) {
        this->rbound[dir] = internalgrav;
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
        IDEFIX_ERROR("SelfGravity:: Axis boundaries are not compatible with "
                     "MPI domain decomposition in X3");
      }
    #endif
    }
  #endif

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
  this->loffset = std::vector<int>(3,0);
  this->roffset = std::vector<int>(3,0);

  if(this->lbound[IDIR] == origin) {
    InitInternalGrid();
  }

  // Update solver when provided
  if(input.CheckEntry("SelfGravity","solver") >= 0) {
    std::string strSolver = input.Get<std::string>("SelfGravity","solver",0);
    if(strSolver.compare("Jacobi")==0) {
      solver = JACOBI;
    } else if(strSolver.compare("BICGSTAB")==0) {
      solver = BICGSTAB;
    } else if(strSolver.compare("PBICGSTAB")==0) {
      solver = PBICGSTAB;
    } else if(strSolver.compare("CG")==0) {
      solver = CG;
    } else if(strSolver.compare("PCG")==0) {
      solver = PCG;
    } else if(strSolver.compare("MINRES")==0) {
      solver = MINRES;
    } else if(strSolver.compare("PMINRES")==0) {
      solver = PMINRES;
    } else {
      try {
        // Try to use the old solver definition with integer (deprecated)
        int s = std::stoi(strSolver);
        if(s<0 || s > 2) throw std::runtime_error("Unknown solver number (should be 0,1 or 2)");
        this->solver = static_cast<GravitySolver> (s);
        IDEFIX_DEPRECATED("The use of integer to define self-gravity solver is deprecated.");
      } catch(const std::exception& e) {
        std::stringstream msg;
        msg << "SelfGravity: Unknown solver \"" << strSolver << "\"."
            << "Use \"Jacobi\", \"BICGSTAB\" or \"PBICGSTAB\"."
            << std::endl;
        IDEFIX_ERROR(msg);
      }
    }
  } else {
    this->solver = BICGSTAB;
  }

  // Enable preconditionner
  if(this->solver==PBICGSTAB || this->solver == PCG || this->solver == PMINRES) {
    this->havePreconditioner = true;
  }

  // Instantiate the bicgstab solver
  if(solver == BICGSTAB || solver == PBICGSTAB) {
    iterativeSolver = new Bicgstab<SelfGravity>(this, &SelfGravity::ComputeLaplacian,
                                                targetError, maxiter,
                                                this->np_tot, this->beg, this->end);
  } else if(solver == CG || solver == PCG) {
    iterativeSolver = new Cg<SelfGravity>(this, &SelfGravity::ComputeLaplacian,
                                                targetError, maxiter,
                                                this->np_tot, this->beg, this->end);
  } else if(solver == MINRES || solver == PMINRES) {
    iterativeSolver = new Minres<SelfGravity>(this, &SelfGravity::ComputeLaplacian,
                                                targetError, maxiter,
                                                this->np_tot, this->beg, this->end);
  } else {
      real step = ComputeJacobiCFL();
      iterativeSolver = new Jacobi<SelfGravity>(this, &SelfGravity::ComputeLaplacian,
                                                targetError, maxiter, step,
                                                this->np_tot, this->beg, this->end);
  }

  // Init preconditionner if needed
  if(havePreconditioner) {
    InitPreconditionner();
  }

  // Arrays initialisation
  this->density = IdefixArray3D<real> ("Density", this->np_tot[KDIR],
                                                  this->np_tot[JDIR],
                                                  this->np_tot[IDIR]);
  // Fill density array with 0
  {
    auto d = this->density;
    idefix_for("InitDensity",0,this->np_tot[KDIR],0,this->np_tot[JDIR],0,this->np_tot[IDIR],
      KOKKOS_LAMBDA (int k, int j, int i) {
        d(k,j,i) = 0.0;
      });
  }

  this->potential = IdefixArray3D<real> ("Potential", this->np_tot[KDIR],
                                                      this->np_tot[JDIR],
                                                      this->np_tot[IDIR]);

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

void SelfGravity::InitInternalGrid() {
  idfx::pushRegion("SelfGravity::InitInternalGrid");
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

void SelfGravity::ShowConfig() {
  idfx::cout << "SelfGravity: Using ";
  switch(solver) {
    case JACOBI:
      idfx::cout << "Jacobi";
      break;
    case BICGSTAB:
      idfx::cout << "unpreconditionned BICGSTAB";
      break;
    case PBICGSTAB:
      idfx::cout << "preconditionned BICGSTAB";
      break;
    case PCG:
      idfx::cout << "preconditionned CG";
      break;
    case CG:
      idfx::cout << "unpreconditionned CG";
      break;
    case MINRES:
      idfx::cout << "unpreconditionned MinRes";
      break;
    case PMINRES:
      idfx::cout << "preconditionned MinRes";
      break;
    default:
      IDEFIX_ERROR("SelfGravity:: Unknown solver");
  }
  idfx::cout << " solver." << std::endl;
  // idfx::cout << "SelfGravity: target L2 norm error=" << targetError << "." << std::endl;
  // idfx::cout << "SelfGravity: 4piG=" << gravCst << "." << std::endl;

  // The setup is periodic if it passes the previous boundary loading
  if(this->isPeriodic == true) {
    idfx::cout << "SelfGravity: Setup is periodic, using specific mass"
              << " re-normalisation." << std::endl;
  }

  if(this->lbound[IDIR] == origin) {
    idfx::cout << "SelfGravity: using origin boundary with " << loffset[IDIR] << " additional "
               << "radial points." << std::endl;
  }

  if(this->skipSelfGravity>1) {
    idfx::cout << "SelfGravity: self-gravity field will be updated every " << skipSelfGravity
               << " cycles." << std::endl;
  }
  iterativeSolver->ShowConfig();
}

void SelfGravity::InitPreconditionner() {
  idfx::pushRegion("SelfGravity::InitPreconditioner");
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

void SelfGravity::InitSolver() {
  idfx::pushRegion("SelfGravity::InitSolver");

  // Loading needed attributes
  IdefixArray3D<real> density = this->density;
  IdefixArray4D<real> Vc = data->hydro.Vc;

  // Initialise the density field
  // todo: check bounds
  int ioffset = this->loffset[IDIR];
  int joffset = this->loffset[JDIR];
  int koffset = this->loffset[KDIR];

  idefix_for("InitDensity", data->beg[KDIR], data->end[KDIR],
                            data->beg[JDIR], data->end[JDIR],
                            data->beg[IDIR], data->end[IDIR],
    KOKKOS_LAMBDA (int k, int j, int i) {
      density(k+koffset, j+joffset, i+ioffset) = Vc(RHO, k, j, i);
    });

  // Deal with the mean issue for periodic density distribution
  if(this->isPeriodic == true) {
    SubstractMeanDensity();  // Remove density mean
  }

  // divide density by preconditionner if we're doing the preconditionned version
  if(havePreconditioner) {
    int ibeg, iend, jbeg, jend, kbeg, kend;
    ibeg = this->beg[IDIR];
    iend = this->end[IDIR];
    jbeg = this->beg[JDIR];
    jend = this->end[JDIR];
    kbeg = this->beg[KDIR];
    kend = this->end[KDIR];
    IdefixArray3D<real> P = this->precond;
    idefix_for("Precond density", kbeg, kend, jbeg, jend, ibeg, iend,
    KOKKOS_LAMBDA (int k, int j, int i) {
      density(k, j, i) = density(k,j,i) / P(k,j,i);
    });
  }

  // Look for Nans in the input field
  int nanDensity = 0;
  idefix_reduce("checkNanVc",0, this->np_tot[KDIR], 0, this->np_tot[JDIR], 0, this->np_tot[IDIR],
    KOKKOS_LAMBDA (int k, int j, int i, int &nnan) {
      if(std::isnan(density(k,j,i))) nnan++;
    }, Kokkos::Sum<int>(nanDensity) // reduction variable
  );
  #ifdef WITH_MPI
    MPI_Allreduce(MPI_IN_PLACE, &nanDensity,1,MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  #endif

  if(nanDensity>0) {
    std::stringstream msg;
    msg << "Input density in self-gravity contains "<< nanDensity <<  " NaNs" << std::endl;
    throw std::runtime_error(msg.str());
  }


  #ifdef DEBUG_GRAVITY
  WriteField(rhoFile, density);
  #endif

  idfx::popRegion();
}

real SelfGravity::ComputeJacobiCFL() {
  idfx::pushRegion("SelfGravity::ComputeCFL");

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


void SelfGravity::SubstractMeanDensity() {
  idfx::pushRegion("SelfGravity::SubstractMeanDensity");

  // Loading needed attributes
  IdefixArray3D<real> density = this->density;
  IdefixArray3D<real> dV = this->dV;

  int ibeg, iend, jbeg, jend, kbeg, kend;
  ibeg = this->beg[IDIR];
  iend = this->end[IDIR];
  jbeg = this->beg[JDIR];
  jend = this->end[JDIR];
  kbeg = this->beg[KDIR];
  kend = this->end[KDIR];

  // Do the reduction on a vector
  MyVector meanDensityVector;

  // Sum the density over the grid, weighted by cell volume
  // and compute the total grid volume as a normalisation constant
  // both stored in a 2D reduction vector
  idefix_reduce("SumWeightedRho",
                kbeg, kend,
                jbeg, jend,
                ibeg, iend,
                KOKKOS_LAMBDA (int k, int j, int i, MyVector &localVector) {
                  localVector.v[0] += density(k,j,i) * dV(k,j,i);
                  localVector.v[1] += dV(k,j,i);
                },
                Kokkos::Sum<MyVector>(meanDensityVector));

  // Reduction on the whole grid
  #ifdef WITH_MPI
  MPI_Allreduce(MPI_IN_PLACE, &meanDensityVector.v, 2, realMPI, MPI_SUM, MPI_COMM_WORLD);
  #endif

  real mean = meanDensityVector.v[0] / meanDensityVector.v[1];

  // Remove the mean value of the density field
  idefix_for("SubstractMeanDensity",
             0, this->np_tot[KDIR],
             0, this->np_tot[JDIR],
             0, this->np_tot[IDIR],
             KOKKOS_LAMBDA (int k, int j, int i) {
               density(k, j, i) -= mean;
             });

  #ifdef DEBUG_GRAVITY
  idfx::cout << "SelfGravity:: Average value " << mean << " substracted to density" << std::endl;
  #endif

  idfx::popRegion();
}



void SelfGravity::ComputeLaplacian(IdefixArray3D<real> array, IdefixArray3D<real> laplacian) {
  idfx::pushRegion("SelfGravity::ComputeLaplacian");

  int ibeg, iend, jbeg, jend, kbeg, kend;
  ibeg = this->beg[IDIR];
  iend = this->end[IDIR];
  jbeg = this->beg[JDIR];
  jend = this->end[JDIR];
  kbeg = this->beg[KDIR];
  kend = this->end[KDIR];
  IdefixArray3D<real> P = this->precond;
  bool havePreconditioner = this->havePreconditioner;

  // Handling boundaries before laplacian calculation
  this->SetBoundaries(array);

  #if DIMENSIONS == 1
  IdefixArray1D<real> dx1 = this->dx[IDIR];
  IdefixArray3D<real> Ax1 = this->A[IDIR];
  IdefixArray3D<real> dV = this->dV;

  // No change in 1D regarding the geometry as the curvature terms are always one
  idefix_for("FiniteDifference", kbeg, kend, jbeg, jend, ibeg, iend,
    KOKKOS_LAMBDA (int k, int j, int i) {
      real gxm, gxp;
      gxm =  2. * (array(k, j, i) - array(k, j, i-1)) / (dx1(i) + dx1(i-1));
      gxp =  2. * (array(k, j, i+1) - array(k, j, i)) / (dx1(i) + dx1(i+1));
      real Delta = (gxp * Ax1(k, j, i+1) - gxm * Ax1(k, j, i)) / dV(k, j, i);
      if(havePreconditioner) {
        Delta=Delta/P(k,j,i);
      }
      laplacian(k, j, i) = Delta;
    });

  #elif DIMENSIONS == 2
  IdefixArray1D<real> dx1 = this->dx[IDIR];
  IdefixArray1D<real> dx2 = this->dx[JDIR];
  IdefixArray3D<real> Ax1 = this->A[IDIR];
  IdefixArray3D<real> Ax2 = this->A[JDIR];
  IdefixArray3D<real> dV = this->dV;
    #if (GEOMETRY == POLAR || GEOMETRY == SPHERICAL)
    IdefixArray1D<real> r = this->x[IDIR];
    #endif

  idefix_for("FiniteDifference", kbeg, kend, jbeg, jend, ibeg, iend,
    KOKKOS_LAMBDA (int k, int j, int i) {
      real h1, h2;
      #if (GEOMETRY == CARTESIAN || GEOMETRY == CYLINDRICAL)
      h1 = 1.;
      h2 = 1.;
      #else
      h1 = 1.;
      h2 = r(i);
      #endif
      real gxm, gxp, gym, gyp;
      gxm =  2. / h1 * (array(k, j, i) - array(k, j, i-1)) / (dx1(i) + dx1(i-1));
      gxp =  2. / h1 * (array(k, j, i+1) - array(k, j, i)) / (dx1(i) + dx1(i+1));
      gym =  2. / h2 * (array(k, j, i) - array(k, j-1, i)) / (dx2(j) + dx2(j-1));
      gyp =  2. / h2 * (array(k, j+1, i) - array(k, j, i)) / (dx2(j) + dx2(j+1));
      real Delta = (gxp * Ax1(k, j, i+1) - gxm * Ax1(k, j, i)
                           + gyp * Ax2(k, j+1, i) - gym * Ax2(k, j, i))
                           / dV(k, j, i);
      if(havePreconditioner) {
        Delta=Delta/P(k,j,i);
      }
      laplacian(k, j, i) = Delta;
    });

  #else
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

  idefix_for("FiniteDifference", kbeg, kend, jbeg, jend, ibeg, iend,
    KOKKOS_LAMBDA (int k, int j, int i) {
      real h1, h2, h3;
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
      real gxm, gxp, gym, gyp, gzm, gzp;
      gxm =  2. / h1 * (array(k, j, i) - array(k, j, i-1)) / (dx1(i) + dx1(i-1));
      gxp =  2. / h1 * (array(k, j, i+1) - array(k, j, i)) / (dx1(i) + dx1(i+1));
      gym =  2. / h2 * (array(k, j, i) - array(k, j-1, i)) / (dx2(j) + dx2(j-1));
      gyp =  2. / h2 * (array(k, j+1, i) - array(k, j, i)) / (dx2(j) + dx2(j+1));
      gzm =  2. / h3 * (array(k, j, i) - array(k-1, j, i)) / (dx3(k) + dx3(k-1));
      gzp =  2. / h3 * (array(k+1, j, i) - array(k, j, i)) / (dx3(k) + dx3(k+1));
      real Delta = (gxp * Ax1(k, j, i+1) - gxm * Ax1(k, j, i)
                           + gyp * Ax2(k, j+1, i) - gym * Ax2(k, j, i)
                           + gzp * Ax3(k+1, j, i) - gzm * Ax3(k, j, i))
                           / dV(k, j, i);
      if(havePreconditioner) {
        Delta=Delta/P(k,j,i);
      }

      laplacian(k, j, i) = Delta;
    });
  #endif


  idfx::popRegion();
}
/*


*/

void SelfGravity::EnforceBoundary(int dir, BoundarySide side, GravityBoundaryType type,
                                  IdefixArray3D<real> &arr) {
  idfx::pushRegion("SelfGravity::EnforceBoundary");

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
        IDEFIX_ERROR("SelfGravity:: No function enrolled to define your own boundary conditions");
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
      std::stringstream msg ("SelfGravity:: Boundary condition type is not yet implemented");
      IDEFIX_ERROR(msg);
    }
  }

  idfx::popRegion();
}

void SelfGravity::EnrollUserDefBoundary(UserDefBoundaryFunc myFunc) {
  this->userDefBoundaryFunc = myFunc;
  this->haveUserDefBoundary = true;
  idfx::cout << "SelfGravity:: User-defined boundary condition has been enrolled" << std::endl;
}

void SelfGravity::SetBoundaries(IdefixArray3D<real> &arr) {
  idfx::pushRegion("SelfGravity::SetBoundaries");

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



#ifdef DEBUG_GRAVITY
// This routine is for debugging purpose
void SelfGravity::WriteField(std::ofstream &stream, IdefixArray3D<real> &in, int index) {
  stream.precision(15);
  stream << std::scientific;// << index;
  for(int i = 0 ; i < this->np_tot[IDIR] ; i++)  {
    stream << in(0,0,i) << "\t";
  }
  stream << std::endl;
}

void SelfGravity::WriteField(std::ofstream &stream, IdefixArray1D<real> &in, int index) {
  stream.precision(15);
  stream << std::scientific;// << index;
  for(int i = 0 ; i < this->np_tot[IDIR] ; i++)  {
    stream << in(i) << "\t";
  }
  stream << std::endl;
}
#endif

void SelfGravity::SolvePoisson() {
  idfx::pushRegion("SelfGravity::SolvePoisson");

  Kokkos::Timer timer;

  elapsedTime -= timer.seconds();

  #ifdef DEBUG_GRAVITY
    rhoFile.open("rho.dat",std::ios::trunc);
    potentialFile.open("potential.dat",std::ios::trunc);
    geometryFile.open("geometry.dat",std::ios::trunc);
    WriteField(geometryFile,this->x[IDIR]);
    WriteField(geometryFile,this->dx[IDIR]);
    WriteField(geometryFile,this->A[IDIR]);
    WriteField(geometryFile,this->dV);
    geometryFile.close();
  #endif

  InitSolver(); // (Re)initialise the solver

  this->nsteps = iterativeSolver->Solve(potential, density);
  if (this->nsteps<0) {
    idfx::cout << "SelfGravity:: BICGSTAB failed, resetting potential" << std::endl;

    // Look for Nans to explain the repetitive failing
    if(data->CheckNan()>0) {
      std::stringstream msg;
      msg << "Nan found after BICGSTAB failed at time " << data->t << std::endl;
      throw std::runtime_error(msg.str());
    }

    // Re-initialise potential
    IdefixArray3D<real> potential = this->potential;

    idefix_for("ResetPotential",
                0, this->np_tot[KDIR],
                0, this->np_tot[JDIR],
                0, this->np_tot[IDIR],
      KOKKOS_LAMBDA (int k, int j, int i) {
        potential(k, j, i) = ZERO_F;
      });

    // Try again !
    this->nsteps = iterativeSolver->Solve(this->potential, density);
    if (this->nsteps<0) {
      IDEFIX_ERROR("SelfGravity:: BICGSTAB failed despite restart");
    }
  }

  currentError = iterativeSolver->GetError();

  #ifdef DEBUG_GRAVITY
  WriteField(potentialFile,potential,n);
  if(this->convStatus == true) {
    idfx::cout << "SelfGravity:: Reached convergence after " << n << " iterations" << std::endl;
    rhoFile.close();
    potentialFile.close();
  } else if(n == maxiter) {
    idfx::cout << "SelfGravity:: Reached max iter" << std::endl;
    rhoFile.close();
    potentialFile.close();
    IDEFIX_WARNING("SelfGravity:: Failed to converge before reaching max iter");
  }
  #endif

  elapsedTime += timer.seconds();
  idfx::popRegion();
}

void SelfGravity::AddSelfGravityPotential(IdefixArray3D<real> &phiP) {
  idfx::pushRegion("SelfGravity::AddSelfGravityPotential");

  // Loading needed data
  IdefixArray3D<real> localPot = phiP;
  IdefixArray3D<real> potential = this->potential;
  real gravCst = this->data->gravity.gravCst;

  // Updating ghost cells before to return potential
  this->SetBoundaries(potential);

  // Adding self-gravity contribution
  int ioffset = this->loffset[IDIR];
  int joffset = this->loffset[JDIR];
  int koffset = this->loffset[KDIR];
  idefix_for("AddSelfGravityPotential", 0, data->np_tot[KDIR],
                                        0, data->np_tot[JDIR],
                                        0, data->np_tot[IDIR],
    KOKKOS_LAMBDA (int k, int j, int i) {
      // Takes into account the unit conversion, scaled by the choice of gravCst
      localPot(k, j, i) += 4.*M_PI*gravCst * potential(k+koffset, j+joffset, i+ioffset);
  });

  idfx::popRegion();
}
