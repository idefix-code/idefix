/*
//@HEADER
// ************************************************************************
//
//                        IDEFIX v 0.0-alpha
//
// ************************************************************************
//@HEADER
*/

#include <limits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sys/time.h>

#define  INDEX_LOOP

#include <Kokkos_Core.hpp>

#include "idefix.hpp"

#define  I_STRIDE
//#define  J_STRIDE
//#define  K_STRIDE
//#define  NO_STRIDE



//#define USE_LEFT_ITERATE





int main( int argc, char* argv[] )
{
  const int NX=64;
  const int NY=64;
  const int NZ=64;
  const int NVAR=5;

#ifdef I_STRIDE
  const int is=1;
  const int ie=NX-2;
#else
  const int is=0;
  const int ie=NX-1;
#endif

#ifdef J_STRIDE
  const int js=1;
  const int je=NY-2;
#else
  const int js=1;
  const int je=NY-2;
#endif

#ifdef K_STRIDE
  const int ks=1;
  const int ke=NZ-2;
#else
  const int ks=0;
  const int ke=NZ-1;
#endif


  int nrepeat=10000;

  Kokkos::initialize( argc, argv );
  {

    Input input = Input();
    Grid grid = Grid(input);
    
  // Allocate the arrays
  IdefixArray4D<real> Q("Q",NVAR,NZ,NY,NX);
  IdefixArray4D<real> V("V",NVAR,NZ,NY,NX);
  IdefixArray4D<real> S("S",NVAR,NZ,NY,NX);
  IdefixArray4D<real>::HostMirror Q_H = Kokkos::create_mirror_view(Q);
  IdefixArray4D<real>::HostMirror V_H = Kokkos::create_mirror_view(V);
  IdefixArray4D<real>::HostMirror S_H = Kokkos::create_mirror_view(S);

  Globals g;

  printf("1\n");
  // Initialize vector on host
  for( int nv = 0; nv < NVAR ; nv++) {
    for( int k = 0; k < NZ ; k++) {
      for(int j = 0; j < NY ; j++) {
        for(int i = 0; i < NX ; i++) {
          Q_H(nv,k,j,i) = i+j+k+nv;
          V_H(nv,k,j,i) = i-j+k-nv;
          S_H(nv,k,j,i) = 0.0;
        }
      }
    }
  }
  printf("2\n");

  // Deep copy host views to device views.
  Kokkos::deep_copy( Q, Q_H );
  Kokkos::deep_copy( V, V_H );
  Kokkos::deep_copy( S, S_H );

  // Timer products.
  Kokkos::Timer timer;

  printf("3\n");


  for ( int repeat = 0; repeat < nrepeat; repeat++ ) {

    idefix_for("product",0,NVAR-1,ks,ke,js,je,is,ie,
                    KOKKOS_LAMBDA (int nv, int k, int j, int i) {
                           #ifdef NO_STRIDE
                           S(nv,k,j,i) = S(nv,k,j,i) + Q(nv,k,j,i) - g.inputParam(0)*(V(nv,k,j,i));
                           #endif

                           #ifdef I_STRIDE
                           S(nv,k,j,i) = S(nv,k,j,i) + Q(nv,k,j,i+1)-Q(nv,k,j,i-1) - g.inputParam(0)*(V(nv,k,j,i+1)+V(nv,k,j,i-1));
                           #endif

                           #ifdef J_STRIDE
                           S(nv,k,j,i) = S(nv,k,j,i) + Q(nv,k,j+1,i)-Q(nv,k,j-1,i) - g.inputParam(0)*(V(nv,k,j+1,i)+V(nv,k,j-1,i));
                           #endif

                           #ifdef K_STRIDE
                           S(nv,k,j,i) = S(nv,k,j,i) + Q(nv,k+1,j,i)-Q(nv,k-1,j,i) - g.inputParam(0)*(V(nv,k+1,j,i)+V(nv,k-1,j,i));
                           #endif
                           //num(0)++;
    });




  }


  Kokkos::fence();

  printf("4\n");

  // Calculate time.
  double time = timer.seconds();
  double th_value=0;

  printf("5\n");
  // Check solution
  Kokkos::deep_copy( S_H, S );

  printf("6\n");
  for( int nv = 0; nv < NVAR ; nv++) {
    for( int k = ks; k <= ke ; k++) {
      for(int j = js; j <= je ; j++) {
        for(int i = is; i <= ie ; i++) {
          #ifdef NO_STRIDE
          th_value=nrepeat*((i+j+k+nv)-0.5*(i-j+k-nv));
          #else
          th_value=nrepeat*(2-(i-j+k-nv));
          #endif

          if( S_H(nv,k,j,i) != th_value)
            printf("Solution mismatch at (nv,i,j,k)=(%d, %d, %d, %d): %g instead of %g \n",nv,i,j,k,S_H(nv,k,j,i),th_value);
        }
      }
    }
  }
  printf("7\n");



  // Calculate bandwidth.
  // Each matrix Q is read once, V is read twise, S is written once
  #ifdef NO_STRIDE
  double Gbytes = 1.0e-9 * double( sizeof(double) * ( 4*NVAR*NZ*NY*NX ) );
  #else
  double Gbytes = 1.0e-9 * double( sizeof(double) * ( 6*NVAR*NZ*NY*NX ) );
  #endif
  // Print results (problem size, time and bandwidth in GB/s).
  //if(num(0)-ref!=0) printf("!!!!!!!!!!! diff num=%ld %% \n",100*num(0)/ref);
  printf( "  NVAR( %d ) Resolution( %d x %d x %d ) nrepeat ( %d ) time( %g s ) bandwidth( %g GB/s )\n",
          NVAR, NX,NY,NZ, nrepeat, time, Gbytes * nrepeat / time );

  }
  Kokkos::finalize();

  return 0;
}
