// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef REAL_TYPES_HPP_
#define REAL_TYPES_HPP_

#include <math.h>



#ifdef SINGLE_PRECISION
  using real = float;
  #ifdef WITH_MPI
    #define realMPI      MPI_FLOAT
  #endif
#else
  using real = double;
  #ifdef WITH_MPI
    #define realMPI     MPI_DOUBLE
  #endif
#endif // SINGLE_PRECISION

// math function
#ifdef SINGLE_PRECISION

#define FMAX(x,y) fmaxf(x,y)
#define FMIN(x,y) fminf(x,y)
#define FABS(x) fabsf(x)
#define TAN(x) tanf(x)
#define SIN(x) sinf(x)
#define COS(x) cosf(x)
#define COPYSIGN(x,y) copysignf(x,y)
#define ISNAN(x) isnanf(x)
#define FMOD(x,y) fmodf(x,y)
#define ZERO_F (0.0f)
#define HALF_F (0.5f)
#define ONE_FOURTH_F (0.25f)
#define ONE_F  (1.0f)
#define TWO_F  (2.0f)
#define THREE_F (3.0f)
#define FOUR_F  (4.0f)

#else

#define FMAX(x,y) fmax(x,y)
#define FMIN(x,y) fmin(x,y)
#define FABS(x) fabs(x)
#define TAN(x) tan(x)
#define SIN(x) sin(x)
#define COS(x) cos(x)
#define COPYSIGN(x,y) copysign(x,y)
#define ISNAN(x) isnan(x)
#define FMOD(x,y) fmod(x,y)
#define ZERO_F (0.0)
#define HALF_F (0.5)
#define ONE_FOURTH_F (0.25)
#define ONE_F   (1.0)
#define TWO_F   (2.0)
#define THREE_F (3.0)
#define FOUR_F  (4.0)


#endif // SINGLE_PRECISION

#endif // REAL_TYPES_HPP_
