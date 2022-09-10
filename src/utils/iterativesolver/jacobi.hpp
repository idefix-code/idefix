// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2022 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef UTILS_ITERATIVESOLVER_JACOBI_HPP_
#define UTILS_ITERATIVESOLVER_JACOBI_HPP_
#include <vector>
#include "idefix.hpp"
#include "vector.hpp"
#include "iterativesolver.hpp"

// The jacobi derives from the iterativesolver class
template <class C>
class Jacobi : public IterativeSolver<C> {
  using LinearFunction = void(C::*) (IdefixArray3D<real> in, IdefixArray3D<real> out);

 public:
  Jacobi(C *parent, LinearFunction func, real error, int maxIter, real step,
           std::vector<int> ntot, std::vector<int> beg, std::vector<int> end);

  int Solve(IdefixArray3D<real> &guess, IdefixArray3D<real> &rhs);
  void ShowConfig();
  void PerformIter();

 private:
  real dt;
};

template <class C>
Jacobi<C>::Jacobi(C *p, LinearFunction f, real error, int maxIter, real step,
           std::vector<int> ntot, std::vector<int> beg, std::vector<int> end) :
            IterativeSolver<C>(p, f, error, maxIter, ntot, beg, end), dt(step) {
              // do nothing
            }

template <class C>
int Jacobi<C>::Solve(IdefixArray3D<real> &guess, IdefixArray3D<real> &rhs) {
  idfx::pushRegion("Jacobi::Solve");
  this->solution = guess;
  this->rhs = rhs;
  // Re-initialise convStatus
  this->convStatus = false;
  // Residual initialisation
  this->SetRes();

  int n = 0;
  while(this->convStatus != true && n < this->maxiter) {
    this->PerformIter();
    n++;
  }

  if(n == this->maxiter) {
    IDEFIX_WARNING("Jacobi:: Failed to converge before reaching max iter."
                  "You should consider to use (P)BICGSTAB.");
  }

  idfx::popRegion();
  return(n);
}

template <class C>
void Jacobi<C>::PerformIter() {
  idfx::pushRegion("Jacobi::PerformIter");

  // Loading needed attributes
  IdefixArray3D<real> solution = this->solution;
  IdefixArray3D<real> res = this->res;
  real dt = this->dt;

  int ibeg, iend, jbeg, jend, kbeg, kend;
  ibeg = this->beg[IDIR];
  iend = this->end[IDIR];
  jbeg = this->beg[JDIR];
  jend = this->end[JDIR];
  kbeg = this->beg[KDIR];
  kend = this->end[KDIR];

  // Perform one iteration of the Jacobi method
  idefix_for("JacIter", kbeg, kend, jbeg, jend, ibeg, iend,
    KOKKOS_LAMBDA (int k, int j, int i) {
      solution(k, j, i) = solution(k, j, i) - dt * res(k,j,i);
    });

  // Update residual
  this->SetRes();

  // Test convergence
  this->TestErrorL2();

  idfx::popRegion();
}

template <class C>
void Jacobi<C>::ShowConfig() {
  idfx::cout << "Jacobi: TargetError: " << this->targetError << std::endl;
  idfx::cout << "Jacobi: Maximum iterations: " << this->maxiter << std::endl;
  idfx::cout << "Jacobi: step: " << this->dt << std::endl;
  return;
}

#endif // UTILS_ITERATIVESOLVER_JACOBI_HPP_
