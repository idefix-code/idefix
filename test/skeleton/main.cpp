#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sys/time.h>
#include <Kokkos_Core.hpp>

#include "idefix.hpp"

// minimal skeleton to use idfx basic functions
void testReduction();


// main function
int main( int argc, char* argv[] )
{
  Kokkos::initialize( argc, argv );
  {
    idfx::initialize();
    idfx::cout << "Beginning of skeleton" << std::endl;

    // This example demonstrate a reduction with kokkos
    testReduction();

    idfx::cout << "End" << std::endl;

  }
Kokkos::finalize();

return 0;

}

void testReduction() {
  const int nx1 = 64;
  const int nx2 = 64;
  const int nx3 = 64;

  // Init an array on device, and image on host
  IdefixArray3D<real> rho = IdefixArray3D<real>("rho",nx3,nx2,nx1);
  IdefixArray3D<real>::HostMirror rhoHost = Kokkos::create_mirror_view(rho);

  // Fill the host array, and compute the theoretical results
  real theoreticalResult = 0;

  for(int k = 0 ; k < nx3 ; k++) {
    for(int j = 0 ; j < nx2 ; j++) {
      for(int i = 0 ; i < nx1 ; i++) {
        rhoHost(k,j,i) = 1.0+k-j;
        theoreticalResult += rhoHost(k,j,i);
      }
    }
  }

  // Send it on the device
  Kokkos::deep_copy(rho, rhoHost);

  real kokkosResult = 0;
  // Try a sum reduction with Kokkos
  Kokkos::parallel_reduce("standard_reduction",
    Kokkos::MDRangePolicy<Kokkos::Rank<3, Kokkos::Iterate::Right, Kokkos::Iterate::Right>>
    ({0,0,0},{nx3,nx2,nx1}),
    KOKKOS_LAMBDA (int k, int j, int i, real &sum) {

        sum += rho(k,j,i);


    }, Kokkos::Sum<real>(kokkosResult) );

  idfx::cout << "Host=" << theoreticalResult << " ; kokkos=" << kokkosResult << std::endl;
}
