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
