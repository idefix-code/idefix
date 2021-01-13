/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Orszag-Tang MHD vortex.

  The Orszag Tang vortex system describes a doubly periodic fluid
  configuration leading to two-dimensional supersonic MHD turbulence.
  Although an analytical solution is not known, its simple and reproducible
  set of initial conditions has made it a widespread benchmark for
   inter-scheme comparison.

  The computational domain is the periodic box \f$[0,2\pi]^D\f$ where
  \c D is the number of spatial dimensions.
  In 2D, the initial condition is given by
  \f[
     \vec{v} = \left(-\sin y,\, \sin x, 0\right) \,,\qquad
     \vec{B} = \left(-\sin y,\, \sin 2x, 0\right) \,,\qquad
     \rho = 25/9\,,\qquad
     p    = 5/3
  \f]

  This test problem does not have any input parameter.

  A snapshot of the solution on a \c 512x512 grid is shown below.

  \image html mhd_ot.02.jpg "Density at t=3.1 (configuration #02)."

  \author A. Mignone (mignone@ph.unito.it)
  \date   April 13, 2014

  \b References
     - "Comparison of some Flux Corrected Transport and TVD numerical
        schemes for hydrodynamic and magnetohydrodynamic problems",
        Toth & Odstrcil, JCP (1996) 128, 82
     - "High-order conservative finite difference GLM-MHD schemes for
        cell-centered MHD", Mignone, Tzeferacos & Bodo, JCP (2010) 229, 5896.
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *
 *********************************************************************** */
{
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
  int i,j,k;
  double *x1 = grid->x[IDIR];
  double *x2 = grid->x[JDIR];
  double *x3 = grid->x[KDIR];

  double x,y,z;
  double B0=1.0/sqrt(4.0*M_PI);

  TOT_LOOP(k,j,i) {
    x=x1[i];
    y=x2[j];
    z=x3[k];

    d->Vc[RHO][k][j][i] = 25.0/(36.0*M_PI);
    d->Vc[PRS][k][j][i] = 5.0/(12.0*M_PI);
    d->Vc[VX1][k][j][i] = - sin(2.0*M_PI*y);
    d->Vc[VX2][k][j][i] =   sin(2.0*M_PI*x);
    d->Vs[BX1s][k][j][i] = -B0*sin(2.0*M_PI*y);
    d->Vs[BX2s][k][j][i] = B0*sin(4.0*M_PI*x);

  }
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
 *
 *********************************************************************** */
{
}
