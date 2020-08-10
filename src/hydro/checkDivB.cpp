
#include "../idefix.hpp"
#include "hydro.hpp"

real Hydro::CheckDivB(DataBlock &data) {
    real divB;
    IdefixArray4D<real> Vs = data.Vs;
    IdefixArray1D<real> dx1 = data.dx[IDIR];
    IdefixArray1D<real> dx2 = data.dx[JDIR];
    IdefixArray1D<real> dx3 = data.dx[KDIR];
    IdefixArray3D<real> Ax1 = data.A[IDIR];
    IdefixArray3D<real> Ax2 = data.A[JDIR];
    IdefixArray3D<real> Ax3 = data.A[KDIR];
    IdefixArray3D<real> dV = data.dV;
    

    Kokkos::parallel_reduce("CheckDivB",
                                Kokkos::MDRangePolicy<Kokkos::Rank<3, Kokkos::Iterate::Right, Kokkos::Iterate::Right>>
                                ({data.beg[KDIR],data.beg[JDIR],data.beg[IDIR]},{data.end[KDIR], data.end[JDIR], data.end[IDIR]}),
                                KOKKOS_LAMBDA (int k, int j, int i, real &divBmax) {
                real dB1,dB2,dB3;

                dB1=dB2=dB3=ZERO_F;
                
                D_EXPAND( dB1=(Ax1(k,j,i+1)*Vs(BX1s,k,j,i+1)-Ax1(k,j,i)*Vs(BX1s,k,j,i)); ,
                          dB2=(Ax2(k,j+1,i)*Vs(BX2s,k,j+1,i)-Ax2(k,j,i)*Vs(BX2s,k,j,i)); ,
                          dB3=(Ax3(k+1,j,i)*Vs(BX3s,k+1,j,i)-Ax3(k,j,i)*Vs(BX3s,k,j,i));  )
                
                /*D_EXPAND( dB1=(Vs(BX1s,k,j,i+1)-Vs(BX1s,k,j,i))/(dx1(i)); ,
                          dB2=(Vs(BX2s,k,j+1,i)-Vs(BX2s,k,j,i))/(dx2(j)); ,
                          dB3=(Vs(BX3s,k+1,j,i)-Vs(BX3s,k,j,i))/(dx3(k));  )*/

                divBmax=FMAX(FABS(D_EXPAND(dB1, +dB2, +dB3))/dV(k,j,i),divBmax);

            }, Kokkos::Max<real>(divB) );

    #ifdef WITH_MPI
    if(idfx::psize>1) {
        MPI_Allreduce(MPI_IN_PLACE, &divB, 1, realMPI, MPI_MAX, MPI_COMM_WORLD);
    }
    #endif
    return(divB);
}


/*
real Hydro::CheckDivB(DataBlock &data) {

    real divB=0;
    IdefixArray4D<real> Vs = data.Vs;
    IdefixArray1D<real> dx1 = data.dx[IDIR];
    IdefixArray1D<real> dx2 = data.dx[JDIR];
    IdefixArray1D<real> dx3 = data.dx[KDIR];

    int iref,jref,kref;
    
    for(int k = data.beg[KDIR] ; k < data.end[KDIR] ; k++) {
        for(int j = data.beg[JDIR] ; j < data.end[JDIR] ; j++) {
            for(int i = data.beg[IDIR] ; i < data.end[IDIR] ; i++) {
                real dB1,dB2,dB3;

                dB1=dB2=dB3=ZERO_F;

                D_EXPAND( dB1=(Vs(BX1s,k,j,i+1)-Vs(BX1s,k,j,i))/(dx1(i)); ,
                          dB2=(Vs(BX2s,k,j+1,i)-Vs(BX2s,k,j,i))/(dx2(j)); ,
                          dB3=(Vs(BX3s,k+1,j,i)-Vs(BX3s,k,j,i))/(dx3(k));  )
                
                if(FABS(D_EXPAND(dB1, +dB2, +dB3)) > divB) {
                    iref=i;
                    jref=j;
                    kref=k;
                    divB=FABS(D_EXPAND(dB1, +dB2, +dB3));
                }
            }
        }
    }
    //idfx::cout << "divB=" << divB << "(i,j,k)=(" << iref << "," << jref << "," << kref << ")" << std::endl;
    return(divB);

}

*/
