#define  PHYSICS                        HD
#define  DIMENSIONS                     2
#define  COMPONENTS                     2
#define  GEOMETRY                       POLAR
#define  BODY_FORCE                     NO
#define  FORCED_TURB                    NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  RK2
#define  DIMENSIONAL_SPLITTING          NO
#define  NTRACER                        0
#define  USER_DEF_PARAMETERS            2

/* -- physics dependent declarations -- */

#define  EOS                            ISOTHERMAL
#define  ENTROPY_SWITCH                 NO
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      RK_LEGENDRE
#define  ROTATING_FRAME                 NO
#define  DUST                           NO

/* -- user-defined parameters (labels) -- */

#define  MACH                           0
#define  NU_VISC                        1

/* [Beg] user-defined constants (do not change this line) */

#define  LIMITER                        MC_LIM
#define  PRINT_TO_SCREEN                YES

/* [End] user-defined constants (do not change this line) */
