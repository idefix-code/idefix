// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef UTILS_ITERATIVESOLVER_ITERATIVESOLVER_HPP_
#define UTILS_ITERATIVESOLVER_ITERATIVESOLVER_HPP_

#include <vector>
#include "idefix.hpp"
#include "vector.hpp"

template <class C>
class IterativeSolver {
 public:
  using LinearFunction = void(C::*) (IdefixArray3D<real> in, IdefixArray3D<real> out);

  IterativeSolver(C *parent, LinearFunction func, real error, int maxIter,
                  std::vector<int> ntot, std::vector<int> beg, std::vector<int> end);

  real GetError();  // return the current error of the solver

  virtual int Solve(IdefixArray3D<real> &guess, IdefixArray3D<real> &rhs) = 0;
  virtual void ShowConfig() = 0;

  // Internal functions (left public for Lambda capture)
  void SetRes();  // Set residual from current guess
  void TestErrorL1();  // Test the convergence status of the current iteration with L1 norm
  void TestErrorL2();  // Test the convergence status of the current iteration with L2 norm
  void TestErrorLINF();  // Test the convergence status of the current iteration with LINF norm
  real ComputeDotProduct(IdefixArray3D<real> mat1, IdefixArray3D<real> mat2);

 protected:
  LinearFunction myFunc;
  real currentError;
  real targetError;
  C *parent;          // parent class
  int maxiter;        // Maximum iteration allowed to achieve convergence
  bool convStatus;    // Convergence status
  bool restart{false};

  std::vector<int> beg;
  std::vector<int> end;
  std::vector<int> ntot;

  IdefixArray3D<real> solution;
  IdefixArray3D<real> rhs;
  IdefixArray3D<real> res; // Residual
};

template <class C>
IterativeSolver<C>::IterativeSolver(C *p, LinearFunction f, real error, int maxiter,
            std::vector<int> ntot, std::vector<int> beg, std::vector<int> end) {
  this->myFunc = f;
  this->targetError = error;
  this->maxiter = maxiter;
  this->beg = beg;
  this->end = end;
  this->ntot = ntot;
  this->restart = false;
  this->parent = p;
  this->currentError = 0;
  this->res = IdefixArray3D<real> ("Residual", this->ntot[KDIR],
                                              this->ntot[JDIR],
                                              this->ntot[IDIR]);
}

template <class C>
void IterativeSolver<C>::TestErrorL1() {
  idfx::pushRegion("IterativeSolver::TestErrorL1");

  // Loading needed attributes
  IdefixArray3D<real> res = this->res;
  IdefixArray3D<real> rhs = this->rhs;

  int ibeg, iend, jbeg, jend, kbeg, kend;
  ibeg = this->beg[IDIR];
  iend = this->end[IDIR];
  jbeg = this->beg[JDIR];
  jend = this->end[JDIR];
  kbeg = this->beg[KDIR];
  kend = this->end[KDIR];

  // Do the reduction on a vector
  MyVector normL1Vector;

  // Sum of absolute residuals over the grid and sum of squared rhs
  // both stored in a 2D reduction vector
  idefix_reduce("absRes",
                kbeg, kend,
                jbeg, jend,
                ibeg, iend,
                KOKKOS_LAMBDA (int k, int j, int i, MyVector &localVector) {
                  localVector.v[0] += FABS(res(k,j,i));
                  localVector.v[1] += rhs(k,j,i) * rhs(k,j,i);
                },
                Kokkos::Sum<MyVector>(normL1Vector));

  // Reduction on the whole grid
  #ifdef WITH_MPI
  MPI_Allreduce(MPI_IN_PLACE, &normL1Vector.v, 2, realMPI, MPI_SUM, MPI_COMM_WORLD);
  #endif

  // Squared error
  this->currentError = sqrt(normL1Vector.v[0] * normL1Vector.v[0] / normL1Vector.v[1]);

  // Checking Nans
  if(std::isnan(currentError)) {
    std::stringstream msg;
    msg << "IterativeSolver:: Current error is Nan." << std::endl;
    throw std::runtime_error(msg.str());
  }

  // Checking convergence
  if(currentError <= this->targetError) {
    this->convStatus = true;
    #ifdef DEBUG_BICGSTAB
    idfx::cout << "IterativeSolver:: Squared, normalized norm L1 is " << currentError
               << " at convergence." << std::endl;
    #endif
  }

  idfx::popRegion();
}

template <class C>
void IterativeSolver<C>::TestErrorL2() {
  idfx::pushRegion("IterativeSolver::TestErrorL2");

  // Loading needed attributes
  IdefixArray3D<real> res = this->res;
  IdefixArray3D<real> rhs = this->rhs;

  int ibeg, iend, jbeg, jend, kbeg, kend;
  ibeg = this->beg[IDIR];
  iend = this->end[IDIR];
  jbeg = this->beg[JDIR];
  jend = this->end[JDIR];
  kbeg = this->beg[KDIR];
  kend = this->end[KDIR];

  // Do the reduction on a vector
  MyVector normL2Vector;

  // Sum of squared residuals over the grid and sum of squared rhs
  // both stored in a 2D reduction vector
  idefix_reduce("SumRes2",
                kbeg, kend,
                jbeg, jend,
                ibeg, iend,
                KOKKOS_LAMBDA (int k, int j, int i, MyVector &localVector) {
                  localVector.v[0] += res(k,j,i) * res(k,j,i);
                  localVector.v[1] += rhs(k,j,i) * rhs(k,j,i);
                },
                Kokkos::Sum<MyVector>(normL2Vector));

  // Reduction on the whole grid
  #ifdef WITH_MPI
  MPI_Allreduce(MPI_IN_PLACE, &normL2Vector.v, 2, realMPI, MPI_SUM, MPI_COMM_WORLD);
  #endif

  // Squared error
  this->currentError = sqrt(normL2Vector.v[0] / normL2Vector.v[1]);
  //idfx::cout << "Error=" << this->currentError << std::endl;

  // Checking Nans
  if(std::isnan(currentError)) {
    std::stringstream msg;
    msg << "IterativeSolver:: Current error is Nan." << std::endl;
    throw std::runtime_error(msg.str());
  }
  // Checking convergence
  if(currentError <= this->targetError) {
    this->convStatus = true;
    #ifdef DEBUG_BICGSTAB
    idfx::cout << "IterativeSolver:: Squared, normalized norm L2 is " << currentError
               << " at convergence." << std::endl;
    #endif
  }

  idfx::popRegion();
}

template <class C>
void IterativeSolver<C>::TestErrorLINF() {
  idfx::pushRegion("IterativeSolver::TestErrorLINF");

  // Loading needed attributes
  IdefixArray3D<real> res = this->res;
  IdefixArray3D<real> rhs = this->rhs;

  int ibeg, iend, jbeg, jend, kbeg, kend;
  ibeg = this->beg[IDIR];
  iend = this->end[IDIR];
  jbeg = this->beg[JDIR];
  jend = this->end[JDIR];
  kbeg = this->beg[KDIR];
  kend = this->end[KDIR];

  real maxRes2, rho2;

  // Searching for the maximum residual over the grid
  idefix_reduce("MaxRes2",
                kbeg, kend,
                jbeg, jend,
                ibeg, iend,
                KOKKOS_LAMBDA (int k, int j, int i, real &localMax) {
                  localMax = std::fmax(res(k,j,i) * res(k,j,i), localMax);
                },
                Kokkos::Max<real>(maxRes2));

  // Sum of squared rhs over the grid
  idefix_reduce("SumDensity2",
                kbeg, kend,
                jbeg, jend,
                ibeg, iend,
                KOKKOS_LAMBDA (int k, int j, int i, real &localSum) {
                  localSum += rhs(k,j,i) * rhs(k,j,i);
                },
                Kokkos::Sum<real>(rho2));

  // Reduction on the whole grid
  #ifdef WITH_MPI
  MPI_Allreduce(MPI_IN_PLACE, &maxRes2, 1, realMPI, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &rho2, 1, realMPI, MPI_SUM, MPI_COMM_WORLD);
  #endif

  // Squared error
  currentError = sqrt(maxRes2/rho2);

  // Checking Nans
  if(std::isnan(currentError)) {
    std::stringstream msg;
    msg << "IterativeSolver:: Current error is Nan." << std::endl;
    throw std::runtime_error(msg.str());
  }
  // Checking convergence
  if(currentError <= this->targetError) {
    this->convStatus = true;
    #ifdef DEBUG_BICGSTAB
    idfx::cout << "IterativeSolver:: Squared, normalized norm LINF is " << currentError
               << " at convergence." << std::endl;
    #endif
  }

  idfx::popRegion();
}

template <class C>
void IterativeSolver<C>::SetRes() {
  idfx::pushRegion("IterativeSolver::SetRes");

  // Loading needed attributes
  IdefixArray3D<real> rhs = this->rhs;
  IdefixArray3D<real> solution = this->solution;
  IdefixArray3D<real> res = this->res;

  // Computing laplacian
  (parent->*myFunc)(solution, res); // We store function output in res to spare workRes array

  int ibeg, iend, jbeg, jend, kbeg, kend;
  ibeg = this->beg[IDIR];
  iend = this->end[IDIR];
  jbeg = this->beg[JDIR];
  jend = this->end[JDIR];
  kbeg = this->beg[KDIR];
  kend = this->end[KDIR];

  // Computing residual
  idefix_for("SetRes", kbeg, kend, jbeg, jend, ibeg, iend,
    KOKKOS_LAMBDA (int k, int j, int i) {
      res(k,j,i) = rhs(k,j,i) - res(k,j,i);
    });

  idfx::popRegion();
}

template <class C>
real IterativeSolver<C>::ComputeDotProduct(IdefixArray3D<real> mat1, IdefixArray3D<real> mat2) {
  idfx::pushRegion("IterativeSolver::ComputeDotProduct");

  int ibeg, iend, jbeg, jend, kbeg, kend;
  ibeg = this->beg[IDIR];
  iend = this->end[IDIR];
  jbeg = this->beg[JDIR];
  jend = this->end[JDIR];
  kbeg = this->beg[KDIR];
  kend = this->end[KDIR];

  real sum;

  idefix_reduce("DotProduct",
                kbeg, kend,
                jbeg, jend,
                ibeg, iend,
                KOKKOS_LAMBDA (int k, int j, int i, real &localSum) {
                  localSum += mat1(k,j,i) * mat2(k,j,i);
                },
                Kokkos::Sum<real>(sum));

  // Reduction on the whole grid
  #ifdef WITH_MPI
  MPI_Allreduce(MPI_IN_PLACE, &sum, 1, realMPI, MPI_SUM, MPI_COMM_WORLD);
  #endif

  idfx::popRegion();

  return sum;
}


template <class C>
real IterativeSolver<C>::GetError() {
  return(currentError);
}

#endif //UTILS_ITERATIVESOLVER_ITERATIVESOLVER_HPP_
