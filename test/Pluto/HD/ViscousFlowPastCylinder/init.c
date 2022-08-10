/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Viscous compressible flow past a cylinder.

  Set initial and boundary conditions for a flow past a cylinder in 2D
  cylindrical polar coordinates \f$(r,\phi)\f$.
  The cylinder has radius 1 and the domain is initially filled with
  constant-density and pressure gas with value \f$\rho = 1,\, p = 1/\Gamma\f$.

  The velocity field is initialized using the potential flow solution
  for an inviscid incompressible flow around a cylinder,
  \f[
      V_r    =  U\left(1 - \frac{1}{r^2}\right)\cos\theta \,,\qquad
      V_\phi = -U\left(1 + \frac{1}{r^2}\right)\sin\theta
  \f]
  where \c U, the far-field velocity, is given by the Mach number.
  The boundary conditions in \c phi are periodic while the outer radial
  boundary is set to inflow for negative values of x while outflow
  for positive values.
  A no-slip boundary condition is used at the fluid-solid interface.

  The flow past the cylinder, no matter how small the viscosity, will
  acquire vorticity in a thin boundary layer adjacent to the cylinder.
  Boundary layer separation may occur leading to the formation of a trailing
  wake behind the cylinder.

  The input parameters are:
  - \c g_inputParam[MACH]: sets the upstream Mach number
  - \c g_inputParam[NU_VISC]: sets the viscosity

  \image html hd_flow_past_cylinder.01.jpg "Density map for configuration #01"
  \image html hd_flow_past_cylinder.04.jpg "Entropy distribution for configuration #04 using 3 levels of refinement."

  \author A. Mignone (mignone@ph.unito.it)
  \date   Aug 18, 2014
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
/*
 *
 *
 *********************************************************************** */
{
  double x, y;
  double U = g_inputParam[MACH];

  v[RHO] = 1.0;
  g_isoSoundSpeed = 1/U;

  //v[PRS] = 1/g_gamma/U/U;
  v[TRC] = 0.0;

  v[VX1] =  cos(x2);
  v[VX2] =  -sin(x2);
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
 *********************************************************************** */
{}

/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid)
/*
 *
 *
 *********************************************************************** */
{
  int   i, j, k, nv;
  double  *x1, *x2, *x3;
  double x, y, r, rnd, Mach;
  double c, s;

  x1 = grid->x[IDIR];
  x2 = grid->x[JDIR];
  x3 = grid->x[KDIR];

  if (side == 0) {    /* -- check solution inside domain -- */
    TOT_LOOP(k,j,i){
    }
  }

  if (side == X1_BEG){  /* -- X1_BEG boundary -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){
        d->Vc[RHO][k][j][i] =  1.0;
        //d->Vc[PRS][k][j][i] =  d->Vc[PRS][k][j][2*IBEG - i - 1];
        d->Vc[VX1][k][j][i] = -d->Vc[VX1][k][j][2*IBEG - i - 1];
        d->Vc[VX2][k][j][i] = -d->Vc[VX2][k][j][2*IBEG - i - 1];
      }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }

  if (side == X1_END){  /* -- X1_END boundary -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){
        c = cos(x2[j]);
        s = sin(x2[j]);
        x = x1[i]*c;
        y = x1[i]*s;
        if (x < 0.0){
/*          rnd = RandomNumber(-1.0, 1.0);  */
          rnd = 0.0;
          Mach = g_inputParam[MACH]*(1.0 + 0.1*rnd);
          d->Vc[RHO][k][j][i] = 1.0;
          //d->Vc[PRS][k][j][i] = 1.0/g_gamma/Mach/Mach;
          d->Vc[VX1][k][j][i] =  c;
          d->Vc[VX2][k][j][i] = -s;
        }else{
          d->Vc[RHO][k][j][i] = d->Vc[RHO][k][j][IEND];
          //d->Vc[PRS][k][j][i] = d->Vc[PRS][k][j][IEND];
          d->Vc[VX1][k][j][i] = d->Vc[VX1][k][j][IEND];
          d->Vc[VX2][k][j][i] = d->Vc[VX2][k][j][IEND];
        }

      }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }
}

#if BODY_FORCE != NO
/* ********************************************************************* */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
/*
 *
 *
 *********************************************************************** */
{
  g[IDIR] = 0.0;
  g[JDIR] = 0.0;
  g[KDIR] = 0.0;
}
/* ********************************************************************* */
double BodyForcePotential(double x1, double x2, double x3)
/*
 *
 *
 *********************************************************************** */
{
  return 0.0;
}
#endif
