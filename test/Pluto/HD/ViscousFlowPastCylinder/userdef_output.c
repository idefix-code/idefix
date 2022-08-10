#include "pluto.h"

/* *************************************************************** */
void ComputeUserVar (const Data *d, Grid *grid)
/*
 *
 *
 ***************************************************************** */
{
  int i, j, k;
  double x1, x2, x3, v[3];
  double ***vx, ***vy, ***vz;

  vx = GetUserVar("vx");
  vy = GetUserVar("vy");

  DOM_LOOP(k,j,i){
    x1 = grid->x[IDIR][i];
    x2 = grid->x[JDIR][j];
    x3 = grid->x[KDIR][k];

    v[0] = d->Vc[VX1][k][j][i];
    v[1] = d->Vc[VX2][k][j][i];

    VectorCartesianComponents(v, x1, x2, x3);
    EXPAND(vx[k][j][i] = v[0];  ,
           vy[k][j][i] = v[1];  ,
           vz[k][j][i] = v[2];)
  }
}
/* ************************************************************* */
void ChangeOutputVar ()
/*
 *
 *
 *************************************************************** */
{
  Image *image;

//  SetDumpVar("vx1",VTK_OUTPUT,NO);
//  SetDumpVar("vx2",VTK_OUTPUT,NO);
}
