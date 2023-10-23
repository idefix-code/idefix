// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef UTILS_ITERATIVESOLVER_CG_HPP_
#define UTILS_ITERATIVESOLVER_CG_HPP_
#include <vector>
#include "idefix.hpp"
#include "vector.hpp"
#include "iterativesolver.hpp"

// The conjugate gradient derives from the iterativesolver class
template <class T>
class Cg : public IterativeSolver<T> {
 public:
  Cg(T &op, real error, int maxIter,
           std::array<int,3> ntot, std::array<int,3> beg, std::array<int,3> end);

  int Solve(IdefixArray3D<real> &guess, IdefixArray3D<real> &rhs);

  void PerformIter();
  void InitSolver();
  void ShowConfig();

 private:
  IdefixArray3D<real> p1; // Search direction for gradient descent
  IdefixArray3D<real> s1; // Search direction for gradient descent
};

template <class T>
Cg<T>::Cg(T &op, real error, int maxiter,
            std::array<int,3> ntot, std::array<int,3> beg, std::array<int,3> end) :
            IterativeSolver<T>(op, error, maxiter, ntot, beg, end) {
  // CG scalars initialisation

  this->p1 = IdefixArray3D<real> ("p1", this->ntot[KDIR],
                                                this->ntot[JDIR],
                                                this->ntot[IDIR]);


  this->s1 = IdefixArray3D<real> ("s1", this->ntot[KDIR],
                                                this->ntot[JDIR],
                                                this->ntot[IDIR]);
}

template <class T>
int Cg<T>::Solve(IdefixArray3D<real> &guess, IdefixArray3D<real> &rhs) {
  idfx::pushRegion("Cg::Solve");
  this->solution = guess;
  this->rhs = rhs;

  // Re-initialise convStatus
  this->convStatus = false;
  this->InitSolver();
  int n = 0;

  while(this->convStatus != true && n < this->maxiter) {
    this->PerformIter();


    n++;
  }

  if(n == this->maxiter) {
    idfx::cout << "Cg:: Reached max iter." << std::endl;
    IDEFIX_ERROR("Cg:: Failed to converge before reaching max iter."
                    "You should consider to use a preconditionner.");
  }

  idfx::popRegion();
  return(n);
}

template <class T>
void Cg<T>::InitSolver() {
  idfx::pushRegion("Cg::InitSolver");
  // Residual initialisation
  this->SetRes();

  Kokkos::deep_copy(this->p1, this->res); // (Re)setting reference residual

  idfx::popRegion();
}

template <class T>
void Cg<T>::PerformIter() {
  idfx::pushRegion("Cg::PerformIter");

  // Loading needed attributes
  auto x = this->solution;
  auto r = this->res;
  auto p1 = this->p1;
  auto s1 = this->s1;

  int ibeg, iend, jbeg, jend, kbeg, kend;
  ibeg = this->beg[IDIR];
  iend = this->end[IDIR];
  jbeg = this->beg[JDIR];
  jend = this->end[JDIR];
  kbeg = this->beg[KDIR];
  kend = this->end[KDIR];

  // ***** Step 1.
  this->linearOperator(p1, s1);

  real rr = this->ComputeDotProduct(r,r);
  //idfx::cout << "rr=" << rr << std::endl;
  real alpha = rr / (this->ComputeDotProduct(p1,s1));

  // Checking for Nans
  if(std::isnan(alpha)) {
    idfx::cout << "Cg:: alpha is nan in step 1." << std::endl;
    this->restart = true;
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

  real beta = this->ComputeDotProduct(r,r) / rr;

  // Checking for Nans
  if(std::isnan(beta)) {
    idfx::cout << "Cg:: beta is nan in step 2." << std::endl;
    idfx::popRegion();
    return;
    // IDEFIX_ERROR("rho is nan in step 1");
  }

  idefix_for("UpdateDir", kbeg, kend, jbeg, jend, ibeg, iend,
    KOKKOS_LAMBDA (int k, int j, int i) {
      p1(k,j,i) = r(k,j,i) + beta * p1(k,j,i);
    });

  idfx::popRegion();
}



template <class T>
void Cg<T>::ShowConfig() {
  idfx::pushRegion("Cg::ShowConfig");
  idfx::cout << "Cg: TargetError: " << this->targetError << std::endl;
  idfx::cout << "Cg: Maximum iterations: " << this->maxiter << std::endl;
  idfx::popRegion();
  return;
}


#endif // UTILS_ITERATIVESOLVER_CG_HPP_
