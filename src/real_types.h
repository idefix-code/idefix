#ifndef REAL_TYPE_H
#define REAL_TYPE_H

#include <math.h>

/**
 * \typedef real_t (alias to float or double)
 */
#ifdef USE_DOUBLE
using real = double;
#else
using real = float;
#endif // USE_DOUBLE

// math function
#if defined(USE_DOUBLE) ||  defined(USE_MIXED_PRECISION)
#define FMAX(x,y) fmax(x,y)
#define FMIN(x,y) fmin(x,y)
#define SQRT(x) sqrt(x)
#define FABS(x) fabs(x)
#define COPYSIGN(x,y) copysign(x,y)
#define ISNAN(x) isnan(x)
#define FMOD(x,y) fmod(x,y)
#define ZERO_F (0.0)
#define HALF_F (0.5)
#define ONE_FOURTH_F (0.25)
#define ONE_F  (1.0)
#define TWO_F  (2.0)
#define FOUR_F  (4.0)
#else
#define FMAX(x,y) fmaxf(x,y)
#define FMIN(x,y) fminf(x,y)
#define SQRT(x) sqrtf(x)
#define FABS(x) fabsf(x)
#define COPYSIGN(x,y) copysignf(x,y)
#define ISNAN(x) isnanf(x)
#define FMOD(x,y) fmodf(x,y)
#define ZERO_F (0.0f)
#define HALF_F (0.5f)
#define ONE_FOURTH_F (0.25f)
#define ONE_F  (1.0f)
#define TWO_F  (2.0f)
#define FOUR_F  (4.0f)
#endif // USE_DOUBLE


#endif // REAL_TYPE_H
