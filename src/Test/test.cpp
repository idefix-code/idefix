#include <cstdio>
#include "idefix.hpp"
#include "test.hpp"

Test::Test(DataBlock &data) {
    // Allocate the arrays
    this->dataHost = DataBlockHost(data);
    this->data=data;
}

Test::Test() {
    // Do nothing
}

void Test::MakeTest(int stride, int nrepeat) {
    // Do nothing
    int NX=data.np_int[0];
    int NY=data.np_int[1];
    int NZ=data.np_int[2];

    int is,ie,js,je,ks,ke;

    

    if(stride==IDIR) {
        is=1;
        ie=NX-1;
    }
    else {
        is=0;
        ie=NX;
    }

    if(stride==JDIR) {
        js=1;
        je=NY-1;
    }
    else {
        js=0;
        je=NY;
    }

    if(stride==KDIR) {
        ks=1;
        ke=NZ-1;
    }
    else {
        ks=0;
        ke=NZ;
    }

    // Initialize vector on host
    for( int k = 0; k < NZ ; k++) {
        for(int j = 0; j < NY ; j++) {
            for(int i = 0; i < NX ; i++) {
            dataHost.Vc(0,k,j,i) = i+j+k;
            dataHost.Vc(1,k,j,i) = i-j+k;
            dataHost.Vc(2,k,j,i) = 0.0;
            }
        }
    }

    // Deep copy host views to device views.
    dataHost.SyncToDevice();

    // Get a local pointer to the needed views since we can't access classes from Kokkos kernels
    IdefixArray4D<real> Vc = data.Vc;
    
    // Timer products.
    Kokkos::Timer timer;


    if(stride<IDIR) {
        for ( int repeat = 0; repeat < nrepeat; repeat++ ) {
            
        idefix_for("product",ks,ke,js,je,is,ie,
                        KOKKOS_LAMBDA (int k, int j, int i) {
                            Vc(2,k,j,i) = Vc(2,k,j,i) + Vc(0,k,j,i) - HALF_F*(Vc(1,k,j,i));
                            
                        });
        }
    }

    if(stride==IDIR) {
        for ( int repeat = 0; repeat < nrepeat; repeat++ ) {

        idefix_for("product",ks,ke,js,je,is,ie,
                        KOKKOS_LAMBDA (int k, int j, int i) {
                            Vc(2,k,j,i) = Vc(2,k,j,i) + Vc(0,k,j,i+1)-Vc(0,k,j,i-1) - HALF_F*(Vc(1,k,j,i+1)+Vc(1,k,j,i-1));
                            
                        });
        }
    }
    if(stride==JDIR) {
        for ( int repeat = 0; repeat < nrepeat; repeat++ ) {

        idefix_for("product",ks,ke,js,je,is,ie,
                        KOKKOS_LAMBDA (int k, int j, int i) {
                            Vc(2,k,j,i) = Vc(2,k,j,i) + Vc(0,k,j+1,i)-Vc(0,k,j-1,i) - HALF_F*(Vc(1,k,j+1,i)+Vc(1,k,j-1,i));                     

                        });
        }
    }
    if(stride==KDIR) {
        for ( int repeat = 0; repeat < nrepeat; repeat++ ) {

        idefix_for("product",ks,ke,js,je,is,ie,
                        KOKKOS_LAMBDA (int k, int j, int i) {
                            Vc(2,k,j,i) = Vc(2,k,j,i) + Vc(0,k+1,j,i)-Vc(0,k-1,j,i) - HALF_F*(Vc(1,k+1,j,i)+Vc(1,k-1,j,i));
                        });
        }
    }
    


    Kokkos::fence();


    // Calculate time.
    double time = timer.seconds();
    double th_value=0;


    // Check solution
    dataHost.SyncFromDevice();


    for( int k = ks; k < ke ; k++) {
        for(int j = js; j < je ; j++) {
            for(int i = is; i < ie ; i++) {
            if(stride<IDIR)
                th_value=nrepeat*((i+j+k)-0.5*(i-j+k));
            else
                th_value=nrepeat*(2-(i-j+k));

            if( dataHost.Vc(2,k,j,i) != th_value)
                printf("Solution mismatch at (i,j,k)=(%d, %d, %d, %d): %g instead of %g \n",i,j,k,dataHost.Vc(2,k,j,i),th_value);
            }
        }
    }



    // Calculate bandwidth.
    // Each matrix Q is read once, V is read twise, S is written once
    double Gbytes;
    if(stride<IDIR) Gbytes = 1.0e-9 * double( sizeof(real) * ( 4*NZ*NY*NX ) );
    else Gbytes = 1.0e-9 * double( sizeof(real) * ( 6*NZ*NY*NX ) );

    // Print results (problem size, time and bandwidth in GB/s).
    //if(num(0)-ref!=0) printf("!!!!!!!!!!! diff num=%ld %% \n",100*num(0)/ref);
    printf( " Stride %d Resolution( %d x %d x %d ) nrepeat ( %d ) time( %g s ) bandwidth( %g GB/s )\n",
            stride, NX,NY,NZ, nrepeat, time, Gbytes * nrepeat / time );

}
