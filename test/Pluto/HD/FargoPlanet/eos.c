#include "pluto.h"

/* **************************************************************** */
void SoundSpeed2 (const State *p, int beg, int end,
                  int pos, Grid *grid)
/*
 *
 *    Define the square of the sound speed for different EOS
 *
 ****************************************************************** */
{
  int  i;

  #if EOS == IDEAL
   for (i = beg; i <= end; i++) p->a2[i] = g_gamma*p->v[i][PRS]/p->v[i][RHO];
  #elif EOS == ISOTHERMAL
  {
    int    j,k;  /* -- used as multidimensional indices -- */
    double *x1, *x2, *x3;

    x1 = grid->x[IDIR];
    x2 = grid->x[JDIR];
    x3 = grid->x[KDIR];

    i = g_i;
    j = g_j;
    k = g_k;

    if (g_dir == IDIR) {
      double R;

      x1 = (pos == FACE_CENTER ? grid->xr[IDIR] : grid->x[IDIR]);
      for (i = beg; i <= end; i++){
        #if GEOMETRY == POLAR
         R = x1[i];
        #elif GEOMETRY == SPHERICAL
         R = x1[i]*sin(x2[j]);
        #endif
        p->a2[i] = g_isoSoundSpeed*g_isoSoundSpeed/R;
      }

    }else if (g_dir == JDIR){
      double R;

      x2 = (pos == FACE_CENTER ? grid->xr[JDIR] : grid->x[JDIR]);
      for (j = beg; j <= end; j++) {
        #if GEOMETRY == POLAR
         R = x1[i];
        #elif GEOMETRY == SPHERICAL
         R = x1[i]*sin(x2[j]);
        #endif
        p->a2[j] = g_isoSoundSpeed*g_isoSoundSpeed/R;
      }
    }else if (g_dir == KDIR){
      double R;

      x3 = (pos == FACE_CENTER ? grid->xr[KDIR] : grid->x[KDIR]);
      for (k = beg; k <= end; k++){
        #if GEOMETRY == POLAR
         R = x1[i];
        #elif GEOMETRY == SPHERICAL
         R = x1[i]*sin(x2[j]);
        #endif
        p->a2[k] = g_isoSoundSpeed*g_isoSoundSpeed/R;
      }
    }
  }
  #else
   print ("! SoundSpeed2: not defined for this EoS\n");
   QUIT_PLUTO(1);
  #endif
}


/* *************************************************************** */
void Enthalpy (real **uprim, real *h, int beg, int end)
/*
 *
 *
 *
 ***************************************************************** */
{
  int i;
  double g_gammar;

  #if EOS == IDEAL
   g_gammar = g_gamma/(g_gamma - 1.0);
   for (i = beg; i <= end; i++){
     h[i] = g_gammar*uprim[i][PRS]/uprim[i][RHO];
   }
  #elif EOS == ISOTHERMAL
   print (" Enthalpy not defined for isothermal EoS\n");
   QUIT_PLUTO(1);
  #endif
}
/* *************************************************************** */
void ENTROPY (real **v, real *s, int is, int ie)
/*
 *
 *
 *
 ***************************************************************** */
{
  int i;
  double rho;

  #if EOS == IDEAL
   for (i = is; i <= ie; i++){
     rho  = v[i][RHO];
     s[i] = v[i][PRS]/pow(rho,g_gamma);
   }
  #elif EOS == ISOTHERMAL || EOS == BAROTROPIC
   print (" Entropy not defined in isothermal or barotropic MHD\n");
   QUIT_PLUTO(1);
  #endif
}
