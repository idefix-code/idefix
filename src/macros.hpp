/* *******************************************************
    Expand dimension- or component-dependent expressions
   *******************************************************  */

/*! \def EXPAND(a,b,c)
    Allows to write component-independent expressions involving vectors 
    by evaluating as many arguments as the value of COMPONENTS. 
    The result is that only the first argument will be compiled in 1D, 
    the first two arguments in 2D and all of them in 3D. 
    As an example: 
    \code
    EXPAND( v[VX1] =  1.0;  ,
            v[VX2] =  5.0;  , 
            v[VX3] = -4.0; )
    \endcode
    becomes
    \code
     v[VX1] = 1.0;   
    \endcode
    when \c COMPONENTS is equal to 1 or
    \code
     v[VX1] = 1.0;  
     v[VX2] = 5.0;
    \endcode
    when \c COMPONENTS is equal to 2 or
    \code
     v[VX1] =  1.0;  
     v[VX2] =  5.0;
     v[VX3] = -4.0;
    \endcode
    when \c COMPONENTS is equal to 3.
*/    
    
/*! \def D_EXPAND(a,b,c)
    Similar to the EXPAND() macro but the expansion depends on DIMENSIONS.  */

/*! \def SELECT(a,b,c)
    Expand only the 1st, 2nd or 3rd argument based on the value of
    COMPONENTS.                                                       */

/*! \def D_SELECT(a,b,c)
    Expand only the 1st, 2nd or 3rd argument based on the value of
    DIMENSIONS.                                                       */

#if COMPONENTS == 1
 #define EXPAND(a,b,c) a
 #define SELECT(a,b,c) a
#endif

#if COMPONENTS == 2
 #define EXPAND(a,b,c) a b
 #define SELECT(a,b,c) b
#endif

#if COMPONENTS == 3
 #define EXPAND(a,b,c) a b c
 #define SELECT(a,b,c) c
#endif

#if DIMENSIONS == 1
 #define D_EXPAND(a,b,c)  a
 #define D_SELECT(a,b,c)  a
 #define IOFFSET 1   /* ** so we know when to loop in boundary zones ** */
 #define JOFFSET 0
 #define KOFFSET 0
#endif

#if DIMENSIONS == 2
 #define D_EXPAND(a,b,c) a b
 #define D_SELECT(a,b,c) b
 #define IOFFSET 1   /* ** so we know when to loop in boundary zones ** */
 #define JOFFSET 1
 #define KOFFSET 0
#endif

#if DIMENSIONS == 3
 #define D_EXPAND(a,b,c) a b c
 #define D_SELECT(a,b,c) c
 #define IOFFSET 1   /* ** so we know when to loop in boundary zones ** */
 #define JOFFSET 1
 #define KOFFSET 1
#endif

/*! \name Spatial averages macros.
    The following set of macros provide a compact way to perform multi-D
    averages from cell centered values to interfaces.
    For instance, \C AVERAGE_X(q,k,j,i) will simply take the
    arithmetic average betwen q(i) and q(i+1) at the i+1/2 interface.
    Likewise, AVERAGE_YZ(q,k,j,i) will produce an average at the
    j+1/2 and k+1/2 edge.
*/
/**@{ */

/* macros for 3D arrays (can't do VA_ARGS because the extra variabl appears first in the parameter list) */
#define AVERAGE_3D_X(q,k,j,i)   (0.5*(q(k,j,i) + q(k,j,i-1)))

#if DIMENSIONS == 1

  #define AVERAGE_3D_Y(q,k,j,i)    (q(0,0,i))
  #define AVERAGE_3D_Z(q,k,j,i)    (q(0,0,i))

  #define AVERAGE_3D_XY(q,k,j,i)   AVERAGE_3D_X(q,0,0,i)
  #define AVERAGE_3D_XZ(q,k,j,i)   AVERAGE_3D_X(q,0,0,i)
  #define AVERAGE_3D_YZ(q,k,j,i)   (q(0,0,i))
 
  #define AVERAGE_3D_XYZ(q,k,j,i)  0.5*(q(0,0,i) + q(0,0,i-1))

#elif DIMENSIONS == 2


  #define AVERAGE_3D_Y(q,k,j,i)   (0.5*(q(k,j,i) + q(k,j-1,i)))
  #define AVERAGE_3D_Z(q,k,j,i)    (q(0,j,i))

  #define AVERAGE_3D_XY(q,k,j,i)   ( 0.25*(  q(k,j,i)   + q(k,j,i-1) \
                                        + q(k,j-1,i) + q(k,j-1,i-1)) )
  #define AVERAGE_3D_XZ(q,k,j,i)   (0.5*(q(0,j,i) + q(0,j,i-1)))
  #define AVERAGE_3D_YZ(q,k,j,i)   (0.5*(q(0,j,i) + q(0,j-1,i)))

  #define AVERAGE_3D_XYZ(q,k,j,i)  (0.25*(  q(0,j,i)   + q(0,j,i-1)        \
                                       + q(0,j-1,i) + q(0,j-1,i-1)))

#elif DIMENSIONS == 3

  #define AVERAGE_3D_Y(q,k,j,i)   (0.5*(q(k,j,i) + q(k,j-1,i)))
  #define AVERAGE_3D_Z(q,k,j,i)   (0.5*(q(k,j,i) + q(k-1,j,i)))

  #define AVERAGE_3D_XY(q,k,j,i)  (0.25*(  q(k,j,i)   + q(k,j,i-1) \
                                       + q(k,j-1,i) + q(k,j-1,i-1)) )
  #define AVERAGE_3D_XZ(q,k,j,i)  (0.25*(  q(k,j,i)   + q(k,j,i-1) \
                                       + q(k-1,j,i) + q(k-1,j,i-1)) )
  #define AVERAGE_3D_YZ(q,k,j,i)  (0.25*(  q(k,j,i)   + q(k,j-1,i) \
                                       + q(k-1,j,i) + q(k-1,j-1,i)) )
  #define AVERAGE_3D_XYZ(q,k,j,i) (0.125*(  q(k,j,i)   + q(k,j,i-1)        \
                                       + q(k,j-1,i) + q(k,j-1,i-1)      \
                                       + q(k-1,j,i)   + q(k-1,j,i-1)    \
                                       + q(k-1,j-1,i) + q(k-1,j-1,i-1)))
#endif

// Same but for 4D views
#define AVERAGE_4D_X(q,n,k,j,i)   (0.5*(q(n,k,j,i) + q(n,k,j,i-1)))

#if DIMENSIONS == 1

  #define AVERAGE_4D_Y(q,n,k,j,i)    (q(n,0,0,i))
  #define AVERAGE_4D_Z(q,n,k,j,i)    (q(n,0,0,i))

  #define AVERAGE_4D_XY(q,n,k,j,i)   AVERAGE_4D_X(q,n,0,0,i)
  #define AVERAGE_4D_XZ(q,n,k,j,i)   AVERAGE_4D_X(q,n,0,0,i)
  #define AVERAGE_4D_YZ(q,n,k,j,i)   (q(n,0,0,i))
 
  #define AVERAGE_4D_XYZ(q,n,k,j,i)  0.5*(q(n,0,0,i) + q(n,0,0,i-1))

#elif DIMENSIONS == 2


  #define AVERAGE_4D_Y(q,n,k,j,i)   (0.5*(q(n,k,j,i) + q(n,k,j-1,i)))
  #define AVERAGE_4D_Z(q,n,k,j,i)    (q(n,0,j,i))

  #define AVERAGE_4D_XY(q,n,k,j,i)   ( 0.25*(  q(n,k,j,i)   + q(n,k,j,i-1) \
                                        + q(n,k,j-1,i) + q(n,k,j-1,i-1)) )
  #define AVERAGE_4D_XZ(q,n,k,j,i)   (0.5*(q(n,0,j,i) + q(n,0,j,i-1)))
  #define AVERAGE_4D_YZ(q,n,k,j,i)   (0.5*(q(n,0,j,i) + q(n,0,j-1,i)))

  #define AVERAGE_4D_XYZ(q,n,k,j,i)  (0.25*(  q(n,0,j,i)   + q(n,0,j,i-1)        \
                                       + q(n,0,j-1,i) + q(n,0,j-1,i-1)))

#elif DIMENSIONS == 3

  #define AVERAGE_4D_Y(q,n,k,j,i)   (0.5*(q(n,k,j,i) + q(n,k,j-1,i)))
  #define AVERAGE_4D_Z(q,n,k,j,i)   (0.5*(q(n,k,j,i) + q(n,k-1,j,i)))

  #define AVERAGE_4D_XY(q,n,k,j,i)  (0.25*(  q(n,k,j,i)   + q(n,k,j,i-1) \
                                       + q(n,k,j-1,i) + q(n,k,j-1,i-1)) )
  #define AVERAGE_4D_XZ(q,n,k,j,i)  (0.25*(  q(n,k,j,i)   + q(n,k,j,i-1) \
                                       + q(n,k-1,j,i) + q(n,k-1,j,i-1)) )
  #define AVERAGE_4D_YZ(q,n,k,j,i)  (0.25*(  q(n,k,j,i)   + q(n,k,j-1,i) \
                                       + q(n,k-1,j,i) + q(n,k-1,j-1,i)) )
  #define AVERAGE_4D_XYZ(q,n,k,j,i) (0.125*(  q(n,k,j,i)   + q(n,k,j,i-1)        \
                                       + q(n,k,j-1,i) + q(n,k,j-1,i-1)      \
                                       + q(n,k-1,j,i)   + q(n,k-1,j,i-1)    \
                                       + q(n,k-1,j-1,i) + q(n,k-1,j-1,i-1)))
#endif

/**@} */
#if MHD == YES
  #if COMPONENTS == 1
    #if HAVE_ENERGY
      #define ARG_EXPAND(a,b,c,f,g,h,e) a,f,e
    #else
      #define ARG_EXPAND(a,b,c,f,g,h,e) a,f
    #endif
  #endif

  #if COMPONENTS == 2
    #if HAVE_ENERGY
      #define ARG_EXPAND(a,b,c,f,g,h,e) a,b,f,g,e
    #else
      #define ARG_EXPAND(a,b,c,f,g,h,e) a,b,f,g
    #endif
  #endif

  #if COMPONENTS == 3
    #if HAVE_ENERGY
      #define ARG_EXPAND(a,b,c,f,g,h,e) a,b,c,f,g,h,e
    #else
      #define ARG_EXPAND(a,b,c,f,g,h,e) a,b,c,f,g,h
    #endif
  #endif
#else
  #if COMPONENTS == 1
    #if HAVE_ENERGY
      #define ARG_EXPAND(a,b,c,e) a,e
    #else
      #define ARG_EXPAND(a,b,c,e) a
    #endif
  #endif

  #if COMPONENTS == 2
    #if HAVE_ENERGY
      #define ARG_EXPAND(a,b,c,e) a,b,e
    #else
      #define ARG_EXPAND(a,b,c,e) a,b
    #endif
  #endif

  #if COMPONENTS == 3
    #if HAVE_ENERGY
      #define ARG_EXPAND(a,b,c,e) a,b,c,e
    #else
      #define ARG_EXPAND(a,b,c,e) a,b,c
    #endif
  #endif
#endif

#define MPI_SAFE_CALL(cmd)                                                    \
   {                                                                          \
        int mpiErrNo = (cmd);                                                 \
        if (MPI_SUCCESS != mpiErrNo) {                                        \
            char msg[MPI_MAX_ERROR_STRING];                                   \
            int len;                                                          \
            std::stringstream stream;                                         \
            MPI_Error_string(mpiErrNo, msg, &len);                            \
            stream << "MPI failed with error code :" << mpiErrNo              \
                                << " " << msg << std::endl;                   \
            IDEFIX_ERROR(stream);                                             \
        }                                                                     \
    }

  
