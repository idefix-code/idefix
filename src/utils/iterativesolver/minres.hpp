// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef UTILS_ITERATIVESOLVER_MINRES_HPP_
#define UTILS_ITERATIVESOLVER_MINRES_HPP_
#include <vector>
#include "idefix.hpp"
#include "vector.hpp"
#include "iterativesolver.hpp"

// The conjugate gradient derives from the iterativesolver class
template <class T>
class Minres : public IterativeSolver<T> {
 public:
  Minres(T &op, real error, int maxIter,
           std::array<int,3> ntot, std::array<int,3> beg, std::array<int,3> end);

  int Solve(IdefixArray3D<real> &guess, IdefixArray3D<real> &rhs);

  void PerformIter();
  void InitSolver();
  void ShowConfig();

 private:
  real alpha;         // MINRES parameter
  real beta;         // MINRES parameter
  real previousError;
  bool shouldReset{false};
  bool firstStep{false};

  IdefixArray3D<real> p0; // Search direction for gradient descent
  IdefixArray3D<real> p1; // Search direction for gradient descent
  IdefixArray3D<real> p2; // Search direction for gradient descent
  IdefixArray3D<real> s0; // Search direction for gradient descent
  IdefixArray3D<real> s1; // Search direction for gradient descent
  IdefixArray3D<real> s2; // Search direction for gradient descent
  IdefixArray3D<real> work1; // work array
  IdefixArray3D<real> work2; // work array
  IdefixArray3D<real> work3; // work array
};

template <class T>
Minres<T>::Minres(T &op, real error, int maxiter,
            std::array<int,3> ntot, std::array<int,3> beg, std::array<int,3> end) :
            IterativeSolver<T>(op, error, maxiter, ntot, beg, end) {
  // MINRES scalars initialisation
  this->alpha = 1.0;
  this->beta = 1.0;



  this->p0 = IdefixArray3D<real> ("p0", this->ntot[KDIR],
                                                this->ntot[JDIR],
                                                this->ntot[IDIR]);
  this->p1 = IdefixArray3D<real> ("p0", this->ntot[KDIR],
                                                this->ntot[JDIR],
                                                this->ntot[IDIR]);
  this->p2 = IdefixArray3D<real> ("p0", this->ntot[KDIR],
                                                this->ntot[JDIR],
                                                this->ntot[IDIR]);

  this->s0 = IdefixArray3D<real> ("p0", this->ntot[KDIR],
                                                this->ntot[JDIR],
                                                this->ntot[IDIR]);
  this->s1 = IdefixArray3D<real> ("p0", this->ntot[KDIR],
                                                this->ntot[JDIR],
                                                this->ntot[IDIR]);
  this->s2 = IdefixArray3D<real> ("p0", this->ntot[KDIR],
                                                this->ntot[JDIR],
                                                this->ntot[IDIR]);
}

template <class T>
int Minres<T>::Solve(IdefixArray3D<real> &guess, IdefixArray3D<real> &rhs) {
  idfx::pushRegion("Minres::Solve");
  this->solution = guess;
  this->rhs = rhs;

  // Re-initialise convStatus
  this->convStatus = false;

  int n = 0;

  while(this->convStatus != true && n < this->maxiter) {
    if(n%20==0) {
      // Reset Minres every 20 step to avoid breakdown
      this->firstStep  = true;
      this->InitSolver();
    }
    this->PerformIter();
    n++;
  }

  if(n == this->maxiter) {
    idfx::cout << "Minres:: Reached max iter." << std::endl;
    IDEFIX_ERROR("Minres:: Failed to converge before reaching max iter."
                    "You should consider to use a preconditionner.");
  }

  idfx::popRegion();
  return(n);
}

template <class T>
void Minres<T>::InitSolver() {
  idfx::pushRegion("Minres::InitSolver");
  // Residual initialisation
  this->SetRes();

  Kokkos::deep_copy(this->p0, this->res); // (Re)setting reference residual
  this->linearOperator(this->p0, this->s0); // (Re)setting associated laplacian
  //Kokkos::deep_copy(this->p1, this->p0);
  //Kokkos::deep_copy(this->s1, this->s0);

  idfx::popRegion();
}

template <class T>
void Minres<T>::PerformIter() {
  idfx::pushRegion("Minres::PerformIter");

  // Loading needed attributes
  auto x = this->solution;
  auto r = this->res;
  auto p0 = this->p0;
  auto p1 = this->p1;
  auto p2 = this->p2;
  auto s0 = this->s0;
  auto s1 = this->s1;
  auto s2 = this->s2;

  Kokkos::deep_copy(p2,p1);
  Kokkos::deep_copy(p1,p0);
  Kokkos::deep_copy(s2,s1);
  Kokkos::deep_copy(s1,s0);

  int ibeg, iend, jbeg, jend, kbeg, kend;
  ibeg = this->beg[IDIR];
  iend = this->end[IDIR];
  jbeg = this->beg[JDIR];
  jend = this->end[JDIR];
  kbeg = this->beg[KDIR];
  kend = this->end[KDIR];

  // ***** Step 1.

  real alpha = this->ComputeDotProduct(r,s1) / (this->ComputeDotProduct(s1,s1));

  // Checking for Nans
  if(std::isnan(alpha)) {
    idfx::cout << "Minres:: alpha is nan in step 1." << std::endl;
    this->shouldReset = true;
    idfx::popRegion();
    return;
    // IDEFIX_ERROR("rho is nan in step 1");
  }

  // ******* Step 2
  idefix_for("UpdateDir", kbeg, kend, jbeg, jend, ibeg, iend,
    KOKKOS_LAMBDA (int k, int j, int i) {
      x(k,j,i) = x(k,j,i) + alpha * p1(k,j,i);
      r(k,j,i) = r(k,j,i) - alpha * s1(k,j,i);
    });

  this->TestErrorL2();
  /*
  if(this->currentError/this->previousError>0.999) {
    this->firstStep = true;
    Kokkos::deep_copy(p0,r);
    this->linearOperator(this->p0, this->s0);
    idfx::cout << "Reset" << std::endl;
    idfx::popRegion();
    return;
  }*/

  this->previousError = this->currentError;
  Kokkos::deep_copy(p0, s1);

  this->linearOperator(s1, s0);

  real beta1 = this->ComputeDotProduct(s0,s1) / this->ComputeDotProduct(s1,s1);

  // Checking for Nans
  if(std::isnan(beta1)) {
    idfx::cout << "Minres:: beta1 is nan in step 2." << std::endl;
    idfx::popRegion();
    return;
    // IDEFIX_ERROR("rho is nan in step 1");
  }

  idefix_for("UpdateDir", kbeg, kend, jbeg, jend, ibeg, iend,
    KOKKOS_LAMBDA (int k, int j, int i) {
      p0(k,j,i) -= beta1 * p1(k,j,i);
      s0(k,j,i) -= beta1 * s1(k,j,i);
  });
  if(!this->firstStep) {
    real beta2 = this->ComputeDotProduct(s0,s2) / this->ComputeDotProduct(s2,s2);
    idefix_for("UpdateDir", kbeg, kend, jbeg, jend, ibeg, iend,
      KOKKOS_LAMBDA (int k, int j, int i) {
        p0(k,j,i) -= beta2 * p2(k,j,i);
        s0(k,j,i) -= beta2 * s2(k,j,i);
      });
  } else {
    this->firstStep = false;
  }

  idfx::popRegion();
}



template <class T>
void Minres<T>::ShowConfig() {
  idfx::pushRegion("Minres::ShowConfig");
  idfx::cout << "Minres: TargetError: " << this->targetError << std::endl;
  idfx::cout << "Minres: Maximum iterations: " << this->maxiter << std::endl;
  idfx::popRegion();
  return;
}


#endif // UTILS_ITERATIVESOLVER_MINRES_HPP_
