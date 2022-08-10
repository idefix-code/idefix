/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Sod shock tube problem.

  The Sod shock tube problem is one of the most used benchmark for
  shock-capturing schemes.
  It is a one-dimensional problem with initial condition given
  by a discontinuity separating two constant states:
  \f[
     \begin{array}{lcll}
   \left(\rho,\, v_x,\, p\right)_L &=&
    \left(1,  0,  1\right)   & \qquad\mathrm{for}\quad x < 0.5
     \\ \noalign{\medskip}
   \left(\rho,\, v_x,\, p\right)_R &=&
    \left(\frac{1}{8},  0,  \frac{1}{10}\right)
      & \qquad\mathrm{for}\quad x > 0.5
   \end{array}
  \f]
  The evolved structured at \c t=0.2 is shown in the panels below and consists of a
  left-going rarefaction wave, a right-going contact discontinutity and a
  right-going shock wave.
  The results shown here were carried with \c PARABOLIC interpolation,
  \c CHARACTERISIC_TRACING time stepping and the \c two_shock Riemann Solver
  on 400 zones (configuration #04).

  \image html hd_sod.04.jpg "Flow profiles for the Sod shock tube at t = 0.2  using configuration #04".

  \author A. Mignone (mignone@ph.unito.it)
  \date   June 08, 2014

  \b References
     - Sod, G. A. (1978).
       "A Survey of Several Finite Difference Methods for Systems of
        Nonlinear Hyperbolic Conservation Laws".
        JCP (1978) 27:1-31

*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
/*
 *
 *********************************************************************** */
{
  #if EOS == IDEAL
   g_gamma = 1.4;
  #endif

  if (fabs(x1) < 0.5) {
    v[RHO] = 1.0;
    v[PRS] = 1.0;
  }else{
    v[RHO] = 0.125;
    v[PRS] = 0.1;
  }
  v[VX1] = 0.0;

}

/* ********************************************************************* */
void InitDomain (Data *d, Grid *grid)
/*!
 * Assign initial condition by looping over the computational domain.
 * Called after the usual Init() function to assign initial conditions
 * on primitive variables.
 * Value assigned here will overwrite those prescribed during Init().
 *
 *
 *********************************************************************** */
{
}

/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/*
 *
 *
 *********************************************************************** */
{
}
/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid)
/*
 *
 *********************************************************************** */
{
}
