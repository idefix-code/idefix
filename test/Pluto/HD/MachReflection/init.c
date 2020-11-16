/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Double Mach reflection test problem.
 
  Sets the initial condition for a planar shock front making an angle of
  \f$ \pi/3 \f$ with a reflecting wall:
  \f[
     \left(\rho,v_x,v_y,p\right) 
     =
    \left\{\begin{array}{ll}
    (1.4,0,0,1) & \qquad {\rm for} \quad x >x_s(0)  \\ \noalign{\medskip}
    (8,8,25,-8.25,116.5) & \qquad {\rm otherwise}
    \end{array}\right.
  \f]
  where \f$x_s(t) = 10t/\sin\alpha + 1/6 + y/\tan\alpha\f$ is the shock 
  position.
  The ideal equation of state with \f$ \Gamma = 1.4\f$ is used.
  The wedge is represented by a reflecting boundary starting at x=1/6 
  along the lower y-boundary.
  As the shock reflects off the lower wall, a complex flow structure 
  develops with two curved reflected shocks propagating at directions
  almost orthogonal to each other and a tangential discontinuity
  separating them. 
  At the wall, a pressure gradient sets up a denser fluid jet 
  propagating along the wall. Kelvin-Helmholtz instability
  patterns may be identified with the rolls developing at
  the slip line. 
  This feature is very sensitive to numerical diffusion and it 
  becomes visible at high resolution and/or low dissipative schemes.

  In the frames below we show configuration \# 02 using a high-order
  finite difference scheme with 5-th order WENO-Z reconstruction and
  \c RK3 time stepping.
  The resolution is 960 x 240.

  \image html dmr_rho.gif "Density contours for the double Mach reflection test" 
  \image html dmr_prs.gif "Pressure contours for the double Mach reflection test" 

  \author A. Mignone (mignone@ph.unito.it)
  \date   May 23, 2014

  \b References
    - "The Numerical Simulation of Two-Dimensional Fluid Flow with Strong Shocks"\n
       Woodward, P.R., & Colella, P., JCP (1984) 54, 115.
    - "The PLUTO Code for computational astrophysics"\n
       Mignone et al., ApJS(2007) 170,228
    - "Comparison of some FLux Corrected Transport and Total variatiion
       diminishing numerical schemes for hydrodynamics and magnetohydrodynamic
       problems" \n
       Toth, G., Odstrcil, D., JCP (1996) 128, 82.

*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *
 *********************************************************************** */
{
  double alpha, xs;

  g_gamma = 1.4;
  alpha   = 1.0/3.0*CONST_PI;
  xs      = 1.0/6.0 + x2/tan(alpha);

  if (x1 > xs){
    us[RHO] = 1.4;
    us[VX1] = 0.0;
    us[VX2] = 0.0;
    us[PRS] = 1.0;
  }else{
    us[RHO] = 8.0;
    us[VX1] =  8.25*sin(alpha);
    us[VX2] = -8.25*cos(alpha);
    us[PRS] = 116.5;
  }      

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
{

}
/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/*! 
 *  Assign user-defined boundary conditions:
 *  - left side (x beg):   constant shocked values
 *  - bottom side (y beg): constant shocked values for 
 *                          <tt> x < 1/6 </tt> and reflective boundary
 *                         otherwise.
 *  - top side (y end):   time-dependent boundary:
 *    for \f$ x < x_s(t) = 10 t/\sin\alpha + 1/6 + 1.0/\tan\alpha \f$ we use 
 *    fixed (post-shock) values. Unperturbed values otherwise.
 *
 *********************************************************************** */
{
  int     i, j, k;
  double  x1, x2, x3;
  double  alpha,xs;
  
  alpha = 60./180.*CONST_PI;

  if (side == X1_BEG){

    X1_BEG_LOOP(k,j,i) {  
      d->Vc[RHO][k][j][i] = 8.0;
      d->Vc[VX1][k][j][i] =   8.25*sin(alpha);
      d->Vc[VX2][k][j][i] = - 8.25*cos(alpha);
      d->Vc[PRS][k][j][i] = 116.5;
    }

  } else if (side == X2_BEG){

    X2_BEG_LOOP(k,j,i) {  
      x1 = grid->x[IDIR][i];
      if (x1 < 1.0/6.0){
        d->Vc[RHO][k][j][i] = 8.0;
        d->Vc[VX1][k][j][i] =   8.25*sin(alpha);
        d->Vc[VX2][k][j][i] = - 8.25*cos(alpha);
        d->Vc[PRS][k][j][i] = 116.5;
      }else{                                   /* reflective boundary */
        d->Vc[RHO][k][j][i] =   d->Vc[RHO][k][2*JBEG - j - 1][i];
        d->Vc[VX1][k][j][i] =   d->Vc[VX1][k][2*JBEG - j - 1][i];
        d->Vc[VX2][k][j][i] = - d->Vc[VX2][k][2*JBEG - j - 1][i];
        d->Vc[PRS][k][j][i] =   d->Vc[PRS][k][2*JBEG - j - 1][i];
      }                    
    }

  } else if (side == X2_END) {

    X2_END_LOOP(k,j,i){

      x1 = grid->x[IDIR][i];
      xs = 10.0*g_time/sin(alpha) + 1.0/6.0 + 1.0/tan(alpha);
      if (x1 < xs){
        d->Vc[RHO][k][j][i] = 8.0;
        d->Vc[VX1][k][j][i] =   8.25*sin(alpha);
        d->Vc[VX2][k][j][i] = - 8.25*cos(alpha);
        d->Vc[PRS][k][j][i] = 116.5;
      }else{                                   /* reflective boundary */
        d->Vc[RHO][k][j][i] = 1.4;
        d->Vc[VX1][k][j][i] = 0.;
        d->Vc[VX2][k][j][i] = 0.;
        d->Vc[PRS][k][j][i] = 1.;
      }                    
    }
  }
}

