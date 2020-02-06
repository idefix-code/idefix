
/* *********************************************************
    Set flow variable indices.
    Extra vector components, when not needed, point to the
    last element (255) of the array stored by startup.c.  
   ********************************************************* */

#define  RHO 0
#define  MX1 1
#define  MX2 (COMPONENTS >= 2 ? 2: 255)
#define  MX3 (COMPONENTS == 3 ? 3: 255)
#if HAVE_ENERGY
  #define ENG  (COMPONENTS + 1)
  #define PRS  ENG
#endif

#define VX1   MX1
#define VX2   MX2
#define VX3   MX3

#define NFLX (1 + COMPONENTS + HAVE_ENERGY)