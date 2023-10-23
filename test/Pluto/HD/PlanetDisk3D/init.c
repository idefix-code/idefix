/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Disk-Planet interaction problem.

  Simulate the interaction of a planet embedded in a disk as described
  in section 3.4 of Mignone et al., A&A (2012) 545, A152.
  This test is a nice benchmark for the FARGO module and the
  \c ROTATING_FRAME switch.
  For testing-purposes  no viscosity is used here.
  The initial condition consists of a locally isothermal configuration
  with temperature profile \f$\propto T^{-1}\f$  yielding a disk vertical
  height to radius of \c 0.05.
  The gravitational potential due to the presence of the star and the
  planet is defined in BodyForcePotential() function.

  The conventions used throught the implementation are the following:

  - \c r  = spherical radius
  - \c R  = cylindrical radius
  - \c z  = cylindrical height
  - \c th = meridional angle

  The test can be carried out in polar (2D or 3D) or spherical (3D)
  coordinates and the following parameters determine the initial configuration:

  -# <tt>g_inputParam[Mstar]</tt>: controls the star mass (in solar masses)
  -# <tt>g_inputParam[Mdisk]</tt>: controls the disk mass (in solar masses)
  -# <tt>g_inputParam[Mplanet]</tt>: sets the planet mass (in earth masses)
  -# <tt>g_inputParam[Viscosity]</tt>: sets the amount of viscosity


  Computation can be carried in the rotating or in the observer's frame
  of reference (\c ROTATING_FRAME to \c YES or \c NO, respectively).
  In particular:

  - Configurations #01 and #02 are in 2D polar coordinates without and with
    FARGO, in the rotating frame.
  - Configurations #03, #04 and #05 are in spherical 3D coordinates with
    and without FARGO in the rotating frame
  - Configurations #06 and #07 are in 2D polar coordinates but in the
    observer's frame
  - Configuration #08 employs static AMR  (grid levels are spatial
    dependent but not dependent on time) in the rotating frame.


  \image html hd_disk_planet.08.png "Density map for configuration #08 using AMR" width=1cm

  \author A. Mignone (mignone@ph.unito.it)
  \date   Aug 16, 2012

  \b References:
     - "A Conservative orbital advection scheme for simulations
        of magnetized shear flows with the PLUTO Code"
        Mignone et al, A&A (2012)
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

#define MIN_DENSITY 1e-8

static void NormalizeDensity (const Data *d, Grid *g);
#if ROTATING_FRAME == NO
 #define g_OmegaZ  0.0
#endif

/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *
 *
 *
 *********************************************************************** */
{
  double r, th, R, z, H, OmegaK, cs;
  double scrh;
  double sigma0=g_inputParam[SIGMA0];
  double sigmaSlope=g_inputParam[SIGMA_SLOPE];
  double h0=g_inputParam[H0];

  #if EOS == IDEAL
   g_gamma = 1.01;
  #endif



  #if GEOMETRY == POLAR
   R  = x1;
   #if DIMENSIONS == 2
    z  = 0.0;
    r  = R;
    th = 0.5*CONST_PI;
   #else
    z  = x3;
    r  = sqrt(R*R + z*z);
    th = atan2(R,z);
   #endif
  #endif

  H      = h0*R;
  OmegaK = 1.0/(R*sqrt(R));
  cs     = H*OmegaK;

  scrh   = (0.5*CONST_PI - th)*r/H;
  // us[RHO] = sigma0/sqrt(2.0*M_PI)/H*pow(R,-sigmaSlope) ;
  us[RHO] = sigma0/sqrt(2.0*M_PI)/(h0*R)*pow(R,-sigmaSlope) * exp(1.0/ (cs*cs) * (1.0/sqrt(R*R+z*z)-1.0/R)) ;

  us[VX1] = us[VX2] = us[VX3] = 0.0;

  us[iVPHI] = R*OmegaK*sqrt( R / sqrt(R*R + z*z) -(2.0+sigmaSlope)*h0*h0);

  if(us[RHO]<1e-6) us[RHO] = 1e-6;

  #if EOS == IDEAL
   us[PRS] = us[RHO]*cs*cs;
  #elif EOS == ISOTHERMAL
    // g_isoSoundSpeed = cs;
    g_isoSoundSpeed = h0;
  #endif


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
  int   i, j, k, nv;
  double *x1, *x2, *x3, R, z, OmegaK, v[256];
  static int do_once = 1;

  double h0 = g_inputParam[H0];
  double sigmaSlope = g_inputParam[SIGMA_SLOPE];
  double sigma0 = g_inputParam[SIGMA0];

  x1 = grid->x[IDIR];
  x2 = grid->x[JDIR];
  x3 = grid->x[KDIR];

  if (side == X1_BEG){
    X1_BEG_LOOP(k,j,i){
      R=x1[i];
      z=x3[k];
      double Vk = 1/sqrt(R);
      double cs2=h0*h0*Vk*Vk;

      d->Vc[RHO][k][j][i] = sigma0/sqrt(2.0*M_PI)/(h0*R)*pow(R,-sigmaSlope) * exp(1.0/(cs2) * (1.0/sqrt(R*R+z*z)-1.0/R)) ;
      d->Vc[VX1][k][j][i] = -d->Vc[VX1][k][j][IBEG];
      d->Vc[VX2][k][j][i] = Vk*sqrt(R/sqrt(R*R + z*z)-(2.0+sigmaSlope)*h0*h0);
      d->Vc[VX3][k][j][i] = 0.0;
    }
  }

  if (side == X1_END){
    X1_END_LOOP(k,j,i){
      NVAR_LOOP(nv)  d->Vc[nv][k][j][i] = d->Vc[nv][k][j][IEND];
      #if GEOMETRY == POLAR
       R = x1[i];
//       d->Vc[iVR][k][j][i] = 0.0;
      #elif GEOMETRY == SPHERICAL
       R = x1[i]*sin(x2[j]);
       d->Vc[iVR][k][j][i]  = 0.0;
       d->Vc[iVTH][k][j][i] = 0.0;
      #endif
      OmegaK = 1.0/(R*sqrt(R));
      d->Vc[iVPHI][k][j][i] = R*(OmegaK - g_OmegaZ);

    }
  }

  if (side == X3_BEG){
    X3_BEG_LOOP(k,j,i){
      R=x1[i];
      z=x3[k];
      double Vk = 1/sqrt(R);
      double cs2=h0*h0*Vk*Vk;
      d->Vc[RHO][k][j][i] = d->Vc[RHO][KBEG][j][i]*exp(1.0/(cs2) * (1.0/sqrt(R*R+z*z) - 1.0/sqrt(R*R+x3[KBEG]*x3[KBEG])));
      d->Vc[VX1][k][j][i] = d->Vc[VX1][KBEG][j][i];
      d->Vc[VX2][k][j][i] = Vk*sqrt(R/sqrt(R*R + z*z)-(1+sigmaSlope)*h0*h0);
      d->Vc[VX3][k][j][i] = -d->Vc[VX3][KBEG][j][i];
    }
  }
  if (side == X3_END){
    X3_END_LOOP(k,j,i){
      R=x1[i];
      z=x3[k];
      double Vk = 1/sqrt(R);
      double cs2=h0*h0*Vk*Vk;
      d->Vc[RHO][k][j][i] = d->Vc[RHO][KEND][j][i]*exp(1.0/(cs2) * (1.0/sqrt(R*R+z*z) - 1.0/sqrt(R*R+x3[KEND]*x3[KEND])));
      d->Vc[VX1][k][j][i] = d->Vc[VX1][KEND][j][i];
      d->Vc[VX2][k][j][i] = Vk*sqrt(R/sqrt(R*R + z*z)-(1+sigmaSlope)*h0*h0);
      d->Vc[VX3][k][j][i] = -d->Vc[VX3][KEND][j][i];
    }
  }


}



#if (BODY_FORCE & VECTOR)
/* ************************************************************************ */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
/*
 *
 *
 *
 *************************************************************************** */
{
  g[IDIR] = 0.0;
  g[JDIR] = 0.0;
  g[KDIR] = 0.0;
}
#endif

#if (BODY_FORCE & POTENTIAL)
/* ************************************************************************ */
double BodyForcePotential(double x1, double x2, double x3)
/*
 *
 *
 *
 *************************************************************************** */
{
  double d, R, r, z, th, x, y, phiplanet, rsm;
  double xp, yp, t, phi;
  double h0=g_inputParam[H0];
  double thicknessSmoothing = g_inputParam[THICKNESS_SMOOTHING];
  double smoothing = h0*thicknessSmoothing;

  #if GEOMETRY == POLAR
   R  = x1;
   #if DIMENSIONS == 2
    z  = 0.0;
    r  = R;
    th = 0.5*CONST_PI;
   #else
    z  = x3;
    r  = sqrt(R*R + z*z);
    th = atan2(R,z);
   #endif
   x  = R*cos(x2);
   y  = R*sin(x2);
  #elif (GEOMETRY == SPHERICAL)
   r  = x1;
   th = x2;
   R = r*sin(th);
   z = r*cos(th);
   x = r*sin(th)*cos(x3);
   y = r*sin(th)*sin(x3);
  #endif

/* ---------------------------------------------
             planet position
   --------------------------------------------- */


 double OmegaZ;
 t = g_time;
 // This is not present in Idefix yet
 //if (g_stepNumber == 2) t += g_dt;
 OmegaZ  = sqrt(1.0 + g_inputParam[Mplanet]/g_inputParam[Mstar]);

 xp = cos(OmegaZ*t);
 yp = sin(OmegaZ*t);


 d = (x-xp)*(x-xp) + (y-yp)*(y-yp) + z*z;
  phiplanet = g_inputParam[Mplanet]/g_inputParam[Mstar]/sqrt(d+smoothing*smoothing);

  phi  = - 1/r;
  phi += -phiplanet;
  phi += g_inputParam[Mplanet]/g_inputParam[Mstar]*(x*xp+y*yp); // /Rp^3 with Rp=1

  return phi;
}
#endif
