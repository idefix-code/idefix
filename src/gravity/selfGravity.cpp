// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include <memory>
#include <string>
#include <vector>

#include "selfGravity.hpp"
#include "dataBlock.hpp"
#include "fluid.hpp"
#include "vector.hpp"
#include "bicgstab.hpp"
#include "cg.hpp"
#include "minres.hpp"
#include "jacobi.hpp"


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
  for (int dir = 0 ; dir < 3 ; dir++) {
    this->lbound[dir] = Laplacian::LaplacianBoundaryType::undefined;
    this->rbound[dir] = Laplacian::LaplacianBoundaryType::undefined;
  }
  for(int dir = 0 ; dir < DIMENSIONS ; dir++) {
    std::string label = std::string("boundary-X")+std::to_string(dir+1)+std::string("-beg");
    std::string boundary = input.Get<std::string>("SelfGravity",label,0);

    if(boundary.compare("nullpot") == 0) {
      this->lbound[dir] = Laplacian::LaplacianBoundaryType::nullpot;
      this->isPeriodic = false;
    } else if(boundary.compare("periodic") == 0) {
      this->lbound[dir] = Laplacian::LaplacianBoundaryType::periodic;
    } else if(boundary.compare("nullgrad") == 0) {
      this->lbound[dir] = Laplacian::LaplacianBoundaryType::nullgrad;
      this->isPeriodic = false;
    } else if(boundary.compare("internalgrav") == 0) {
      this->lbound[dir] = Laplacian::LaplacianBoundaryType::internalgrav;
      this->isPeriodic = false;
    } else if(boundary.compare("userdef") == 0) {
      this->lbound[dir] = Laplacian::LaplacianBoundaryType::userdef;
      this->isPeriodic = false;
    } else if(boundary.compare("axis") == 0) {
      this->lbound[dir] = Laplacian::LaplacianBoundaryType::axis;
      this->isPeriodic = false;
    } else if(boundary.compare("origin") == 0) {
      this->lbound[dir] = Laplacian::LaplacianBoundaryType::origin;
      this->isPeriodic = false;
      // origin only compatible with spherical & axis=IDIR
      #if GEOMETRY != SPHERICAL
        IDEFIX_ERROR("Origin boundary conditions are working in spherical coordinates");
      #endif
      if(dir != IDIR) {
        IDEFIX_ERROR("Origin boundary conditions are meaningful only on the X1 direction");
      }
    } else {
      std::stringstream msg;
      msg << "SelfGravity:: Unknown boundary type " << boundary;
      IDEFIX_ERROR(msg);
    }

    label = std::string("boundary-X")+std::to_string(dir+1)+std::string("-end");
    boundary = input.Get<std::string>("SelfGravity",label,0);
    if(boundary.compare("nullpot") == 0) {
      this->rbound[dir] = Laplacian::LaplacianBoundaryType::nullpot;
      this->isPeriodic = false;
    } else if(boundary.compare("periodic") == 0) {
      this->rbound[dir] = Laplacian::LaplacianBoundaryType::periodic;
    } else if(boundary.compare("nullgrad") == 0) {
      this->rbound[dir] = Laplacian::LaplacianBoundaryType::nullgrad;
      this->isPeriodic = false;
    } else if(boundary.compare("internalgrav") == 0) {
      this->rbound[dir] = Laplacian::LaplacianBoundaryType::internalgrav;
      this->isPeriodic = false;
    } else if(boundary.compare("userdef") == 0) {
      this->rbound[dir] = Laplacian::LaplacianBoundaryType::userdef;
      this->isPeriodic = false;
    } else if(boundary.compare("axis") == 0) {
      this->rbound[dir] = Laplacian::LaplacianBoundaryType::axis;
      this->isPeriodic = false;
    } else {
      std::stringstream msg;
      msg << "SelfGravity:: Unknown boundary type " << boundary;
      IDEFIX_ERROR(msg);
    }
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

  // Make the Laplacian operator
  laplacian = std::make_unique<Laplacian>(data, lbound, rbound, this->havePreconditioner );

  np_tot = laplacian->np_tot;

  // Instantiate the bicgstab solver
  if(solver == BICGSTAB || solver == PBICGSTAB) {
    iterativeSolver = new Bicgstab<Laplacian>(*laplacian.get(), targetError, maxiter,
                                              laplacian->np_tot, laplacian->beg, laplacian->end);
  } else if(solver == CG || solver == PCG) {
    iterativeSolver = new Cg<Laplacian>(*laplacian.get(), targetError, maxiter,
                                        laplacian->np_tot, laplacian->beg, laplacian->end);
  } else if(solver == MINRES || solver == PMINRES) {
    iterativeSolver = new Minres<Laplacian>(*laplacian.get(),
                                  targetError, maxiter,
                                  laplacian->np_tot, laplacian->beg, laplacian->end);
  } else {
      real step = laplacian->ComputeCFL();
      iterativeSolver = new Jacobi<Laplacian>(*laplacian.get(), targetError, maxiter, step,
                                              laplacian->np_tot, laplacian->beg, laplacian->end);
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

  if(this->lbound[IDIR] == Laplacian::LaplacianBoundaryType::origin) {
    idfx::cout << "SelfGravity: using origin boundary with " << laplacian->loffset[IDIR]
               << " additional radial points." << std::endl;
  }

  if(this->skipSelfGravity>1) {
    idfx::cout << "SelfGravity: self-gravity field will be updated every " << skipSelfGravity
               << " cycles." << std::endl;
  }
  iterativeSolver->ShowConfig();
}



void SelfGravity::InitSolver() {
  idfx::pushRegion("SelfGravity::InitSolver");

  // Loading needed attributes
  IdefixArray3D<real> density = this->density;
  IdefixArray4D<real> Vc = data->hydro->Vc;

  // Initialise the density field
  // todo: check bounds
  int ioffset = laplacian->loffset[IDIR];
  int joffset = laplacian->loffset[JDIR];
  int koffset = laplacian->loffset[KDIR];

  idefix_for("InitDensity", data->beg[KDIR], data->end[KDIR],
                            data->beg[JDIR], data->end[JDIR],
                            data->beg[IDIR], data->end[IDIR],
    KOKKOS_LAMBDA (int k, int j, int i) {
      density(k+koffset, j+joffset, i+ioffset) = Vc(RHO, k, j, i);
    });

  // Make sure that dust mass contributes to the self-gravitating field
  if(data->haveDust) {
    for(int i = 0 ; i < data->dust.size() ; i++) {
      IdefixArray4D<real> VcDust = data->dust[i]->Vc;
      idefix_for("InitDustDensity", data->beg[KDIR], data->end[KDIR],
                                    data->beg[JDIR], data->end[JDIR],
                                    data->beg[IDIR], data->end[IDIR],
        KOKKOS_LAMBDA (int k, int j, int i) {
          density(k+koffset, j+joffset, i+ioffset) += VcDust(RHO, k, j, i);
      });
    }
  }

  // Deal with the mean issue for periodic density distribution
  if(this->isPeriodic == true) {
    SubstractMeanDensity();  // Remove density mean
  }

  // divide density by preconditionner if we're doing the preconditionned version
  if(havePreconditioner) {
    int ibeg, iend, jbeg, jend, kbeg, kend;
    ibeg = laplacian->beg[IDIR];
    iend = laplacian->end[IDIR];
    jbeg = laplacian->beg[JDIR];
    jend = laplacian->end[JDIR];
    kbeg = laplacian->beg[KDIR];
    kend = laplacian->end[KDIR];
    IdefixArray3D<real> P = laplacian->precond;
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

  idfx::popRegion();
}


void SelfGravity::SubstractMeanDensity() {
  idfx::pushRegion("SelfGravity::SubstractMeanDensity");

  // Loading needed attributes
  IdefixArray3D<real> density = this->density;
  IdefixArray3D<real> dV = laplacian->dV;

  int ibeg, iend, jbeg, jend, kbeg, kend;
  ibeg = laplacian->beg[IDIR];
  iend = laplacian->end[IDIR];
  jbeg = laplacian->beg[JDIR];
  jend = laplacian->end[JDIR];
  kbeg = laplacian->beg[KDIR];
  kend = laplacian->end[KDIR];

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

  idfx::popRegion();
}

void SelfGravity::EnrollUserDefBoundary(Laplacian::UserDefBoundaryFunc myFunc) {
  laplacian->EnrollUserDefBoundary(myFunc);
}



void SelfGravity::SolvePoisson() {
  idfx::pushRegion("SelfGravity::SolvePoisson");

  Kokkos::Timer timer;

  elapsedTime -= timer.seconds();

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


  elapsedTime += timer.seconds();
  idfx::popRegion();
}

void SelfGravity::AddSelfGravityPotential(IdefixArray3D<real> &phiP) {
  idfx::pushRegion("SelfGravity::AddSelfGravityPotential");

  // Loading needed data
  IdefixArray3D<real> localPot = phiP;
  IdefixArray3D<real> potential = this->potential;
  real gravCst = this->data->gravity->gravCst;

  // Updating ghost cells before to return potential
  laplacian->SetBoundaries(potential);

  // Adding self-gravity contribution
  int ioffset = laplacian->loffset[IDIR];
  int joffset = laplacian->loffset[JDIR];
  int koffset = laplacian->loffset[KDIR];
  idefix_for("AddSelfGravityPotential", 0, data->np_tot[KDIR],
                                        0, data->np_tot[JDIR],
                                        0, data->np_tot[IDIR],
    KOKKOS_LAMBDA (int k, int j, int i) {
      // Takes into account the unit conversion, scaled by the choice of gravCst
      localPot(k, j, i) += 4.*M_PI*gravCst * potential(k+koffset, j+joffset, i+ioffset);
  });

  idfx::popRegion();
}
