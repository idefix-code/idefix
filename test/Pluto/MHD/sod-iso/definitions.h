#define  PHYSICS                        MHD
#define  DIMENSIONS                     1
#define  COMPONENTS                     3
#define  GEOMETRY                       CARTESIAN
#define  BODY_FORCE                     NO
#define  FORCED_TURB                    NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  RK2
#define  DIMENSIONAL_SPLITTING          NO
#define  NTRACER                        0
#define  USER_DEF_PARAMETERS            17

/* -- physics dependent declarations -- */

#define  EOS                            ISOTHERMAL
#define  ENTROPY_SWITCH                 NO
#define  DIVB_CONTROL                   NO
#define  BACKGROUND_FIELD               NO
#define  AMBIPOLAR_DIFFUSION            NO
#define  RESISTIVITY                    NO
#define  HALL_MHD                       NO
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      NO
#define  ROTATING_FRAME                 NO

/* -- user-defined parameters (labels) -- */

#define  RHO_LEFT                       0
#define  VX_LEFT                        1
#define  VY_LEFT                        2
#define  VZ_LEFT                        3
#define  BY_LEFT                        4
#define  BZ_LEFT                        5
#define  PR_LEFT                        6
#define  RHO_RIGHT                      7
#define  VX_RIGHT                       8
#define  VY_RIGHT                       9
#define  VZ_RIGHT                       10
#define  BY_RIGHT                       11
#define  BZ_RIGHT                       12
#define  PR_RIGHT                       13
#define  BX_CONST                       14
#define  GAMMA_EOS                      15
#define  DIVIDE_BY_4PI                  16

/* [Beg] user-defined constants (do not change this line) */

#define  LIMITER                        VANLEER_LIM

/* [End] user-defined constants (do not change this line) */
