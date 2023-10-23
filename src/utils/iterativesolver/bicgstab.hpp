// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef UTILS_ITERATIVESOLVER_BICGSTAB_HPP_
#define UTILS_ITERATIVESOLVER_BICGSTAB_HPP_
#include <vector>
#include "idefix.hpp"
#include "vector.hpp"
#include "iterativesolver.hpp"

// The bicgstab derives from the iterativesolver class
template <class T>
class Bicgstab : public IterativeSolver<T> {
 public:
  Bicgstab(T &op, real error, int maxIter,
           std::array<int,3> ntot, std::array<int,3> beg, std::array<int,3> end);

  int Solve(IdefixArray3D<real> &guess, IdefixArray3D<real> &rhs);

  void PerformIter();
  void InitSolver();
  void ShowConfig();

 private:
  real rho;           // BICSTAB parameter
  real alpha;         // BICGSTAB parameter
  real omega;         // BICGSTAB parameter


  IdefixArray3D<real> res0; // Reference (initial) residual
  IdefixArray3D<real> dir; // Search direction for gradient descent
  IdefixArray3D<real> work1; // work array
  IdefixArray3D<real> work2; // work array
  IdefixArray3D<real> work3; // work array
};

template <class T>
Bicgstab<T>::Bicgstab(T &op, real error, int maxiter,
            std::array<int,3> ntot, std::array<int,3> beg, std::array<int,3> end) :
            IterativeSolver<T>(op, error, maxiter, ntot, beg, end) {
  // BICGSTAB scalars initialisation
  this->rho = 1.0;
  this->alpha = 1.0;
  this->omega = 1.0;



  this->dir = IdefixArray3D<real> ("Direction", this->ntot[KDIR],
                                                this->ntot[JDIR],
                                                this->ntot[IDIR]);

  this->res0 = IdefixArray3D<real> ("InitialResidual", this->ntot[KDIR],
                                                        this->ntot[JDIR],
                                                        this->ntot[IDIR]);

  this->work1 = IdefixArray3D<real> ("WorkingArray1", this->ntot[KDIR],
                                                      this->ntot[JDIR],
                                                      this->ntot[IDIR]);

  this->work2 = IdefixArray3D<real> ("WorkingArray2", this->ntot[KDIR],
                                                      this->ntot[JDIR],
                                                      this->ntot[IDIR]);

  this->work3 = IdefixArray3D<real> ("WorkingArray3", this->ntot[KDIR],
                                                      this->ntot[JDIR],
                                                      this->ntot[IDIR]);
}

template <class T>
int Bicgstab<T>::Solve(IdefixArray3D<real> &guess, IdefixArray3D<real> &rhs) {
  idfx::pushRegion("Bicgstab::Solve");
  this->solution = guess;
  this->rhs = rhs;

  // Re-initialise convStatus
  this->convStatus = false;

  this->InitSolver();

  int n = 0;
  while(this->convStatus != true && n < this->maxiter) {
    this->PerformIter();
    if(this->restart) {
      this->restart=false;
      // Resetting parameters
      this->rho = 1.0;
      this->alpha = 1.0;
      this->omega = 1.0;
      n = -1;
      idfx::popRegion();
      return(n);
    }
    n++;
  }

  if(n == this->maxiter) {
    idfx::cout << "Bicgstab:: Reached max iter." << std::endl;
    IDEFIX_WARNING("Bicgstab:: Failed to converge before reaching max iter."
                    "You should consider to use a preconditionner.");
  }

  idfx::popRegion();
  return(n);
}

template <class T>
void Bicgstab<T>::InitSolver() {
  idfx::pushRegion("Bicgstab::InitSolver");
  // Residual initialisation
  this->SetRes();

  Kokkos::deep_copy(this->res0, this->res); // (Re)setting reference residual
  Kokkos::deep_copy(this->dir, this->res); // (Re)setting initial searching direction
  this->linearOperator(this->dir, this->work1); // (Re)setting associated laplacian

  // // Resetting parameters
  // this->rho = 1.0;
  // this->alpha = 1.0;
  // this->omega = 1.0;

  idfx::popRegion();
}

template <class T>
void Bicgstab<T>::PerformIter() {
  idfx::pushRegion("Bicgstab::PerformIter");

  // Loading needed attributes
  IdefixArray3D<real> solution = this->solution;
  IdefixArray3D<real> res = this->res;
  IdefixArray3D<real> dir = this->dir;
  IdefixArray3D<real> res0 = this->res0; // Reference residual, do not evolve through the loop

  // The following variables are named following wikipedia's page nomenclature of BICGSTAB algorithm
  IdefixArray3D<real> v = this->work1; // Working array, for laplacian dir calculation
  IdefixArray3D<real> s = this->work2; // Working array, for intermediate dir calculation
  IdefixArray3D<real> t = this->work3; // Working array, for laplacian intermediate dir calculation
  real omega;
  real &alpha = this->alpha;
  real &rhoOld = this->rho;
  real &omegaOld = this->omega;

  int ibeg, iend, jbeg, jend, kbeg, kend;
  ibeg = this->beg[IDIR];
  iend = this->end[IDIR];
  jbeg = this->beg[JDIR];
  jend = this->end[JDIR];
  kbeg = this->beg[KDIR];
  kend = this->end[KDIR];

  // ***** Step 1.
  real rho = this->ComputeDotProduct(res0, res);
  #ifdef DEBUG_BICGSTAB
    idfx::cout << "rho= " << rho;
  #endif

  // Checking for Nans
  if(std::isnan(rho)) {
    idfx::cout << "Bicgstab:: rho is nan in step 1." << std::endl;
    this->restart = true;
    idfx::popRegion();
    return;
    // IDEFIX_ERROR("rho is nan in step 1");
  }

  // ***** Step 2.
  real beta = rho / rhoOld * alpha / omegaOld;
  #ifdef DEBUG_BICGSTAB
    idfx::cout << " ; beta= " << beta;
  #endif

  // Checking for Nans
  if(std::isnan(beta)) {
    idfx::cout << "Bicgstab:: beta is nan in step 2." << std::endl;
    // IDEFIX_ERROR("beta is nan in step 2");
    this->restart = true;
    idfx::popRegion();
    return;
  }

  // ***** Step 3.
  idefix_for("UpdateDir", kbeg, kend, jbeg, jend, ibeg, iend,
    KOKKOS_LAMBDA (int k, int j, int i) {
      dir(k,j,i) = res(k,j,i) + beta * (dir(k,j,i) - omegaOld * v(k,j,i));
    });

  // From now dir is updated

  // ***** Step 4.
  this->linearOperator(dir, v);

  // from now v is updated (laplacian of dir)

  // ******* Step 5.

  alpha = rho / this->ComputeDotProduct(res0, v);
  #ifdef DEBUG_BICGSTAB
    idfx::cout << " ; alpha= " << alpha;
  #endif

  // Checking Nans
  if(std::isnan(alpha)) {
    idfx::cout << "Bicgstab:: alpha is nan in step 5." << std::endl;
    // IDEFIX_ERROR("alpha is nan in step 5");
    this->restart = true;
    idfx::popRegion();
    return;
  }

  // *********** Step 6.
  // Assumes solution = x_i-1
  idefix_for("FirstUpdatePot", kbeg, kend, jbeg, jend, ibeg, iend,
    KOKKOS_LAMBDA (int k, int j, int i) {
      solution(k,j,i) = solution(k,j,i) + alpha * dir(k,j,i);
    });

  // From here solution = h_i

  // ********** Step.7.
  // Store current residual
  Kokkos::deep_copy(s, res); // s is momentarily oldRes to recycle arrays

  // Update residual
  this->SetRes();

  // Test intermediate guess h_i
  this->TestErrorL2();

  // The loop continues if no convergence
  if(this->convStatus == false) {
    // ***************** Step. 8.
    idefix_for("FillIntermediateDir", kbeg, kend, jbeg, jend, ibeg, iend,
      KOKKOS_LAMBDA (int k, int j, int i) {
        s(k,j,i) = s(k,j,i) - alpha * v(k,j,i); // s in RHS is oldRes before update
      });

    // From here s is updated

    // ************** Step 9.
    this->linearOperator(s, t);

    // From here t is updated

    // ************* Step 10.
    omega = this->ComputeDotProduct(t, s) / this->ComputeDotProduct(t, t);

    #ifdef DEBUG_BICGSTAB
      idfx::cout << " ; omega= " << omega;
    #endif

    // Checking Nans
    if(std::isnan(omega)) {
      idfx::cout << "Bicgstab:: omega is nan in step 10." << std::endl;
      // IDEFIX_ERROR("omega is nan in step 10");
      this->restart = true;
      idfx::popRegion();
      return;
    }

    // ************ Step 11.
    // solution is h_i from step 6.
    idefix_for("SecondUpdatePot", kbeg, kend, jbeg, jend, ibeg, iend,
      KOKKOS_LAMBDA (int k, int j, int i) {
        solution(k,j,i) = solution(k,j,i) + omega * s(k,j,i);
      });

    // From here, solution = x_i

    // *********** Step 12.
    // Update residual
    this->SetRes();

    // Test final guess x_i
    this->TestErrorL2();

    // Last task if no convergence : update res
    if(this->convStatus == false) {
      idefix_for("UpdateRes", kbeg, kend, jbeg, jend, ibeg, iend,
        KOKKOS_LAMBDA (int k, int j, int i) {
          res(k,j,i) = s(k,j,i) - omega * t(k,j,i);
        });

      // Saving rho and omega for next iteration
      rhoOld = rho;
      omegaOld = omega;
    } else {
        #ifdef DEBUG_BICGSTAB
        idfx::cout << "Bicgstab:: BICGSTAB converged in step 12." << std::endl;
        #endif
    }
  }  else {
      #ifdef DEBUG_BICGSTAB
      idfx::cout << "Bicgstab:: BICGSTAB converged in step 7." << std::endl;
      #endif
  }
  #ifdef DEBUG_BICGSTAB
  idfx::cout << std::endl;
  #endif

  idfx::popRegion();
}



template <class T>
void Bicgstab<T>::ShowConfig() {
  idfx::pushRegion("Bicgstab::ShowConfig");
  idfx::cout << "Bicgstab: TargetError: " << this->targetError << std::endl;
  idfx::cout << "Bicgstab: Maximum iterations: " << this->maxiter << std::endl;
  idfx::popRegion();
  return;
}


#endif // UTILS_ITERATIVESOLVER_BICGSTAB_HPP_
