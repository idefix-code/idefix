#include "idefix.hpp"
#include "grid.hpp"
#include "gridHost.hpp"

GridHost::GridHost(Grid &grid, Input &input) {
    std::cout << "building gridHost\n";
    for(int dir = 0 ; dir < 3 ; dir++) {
        std::cout << "nghost[" << dir << "]=" << grid.nghost[dir] << "\n";
        nghost[dir] = grid.nghost[dir];
        np_tot[dir] = grid.np_tot[dir];
        np_int[dir] = grid.np_int[dir];

        lbound[dir] = grid.lbound[dir]; 
        rbound[dir] = grid.rbound[dir];

        lbeg[dir] = grid.lbeg[dir];
        lend[dir] = grid.lend[dir];

    }

    // Create mirrors on host
    for(int dir = 0 ; dir < 3 ; dir++) {
        x[dir] = Kokkos::create_mirror_view(grid.x[dir]);
        xr[dir] = Kokkos::create_mirror_view(grid.xr[dir]);
        xl[dir] = Kokkos::create_mirror_view(grid.xl[dir]);
        dx[dir] = Kokkos::create_mirror_view(grid.dx[dir]);
        A[dir] = Kokkos::create_mirror_view(grid.A[dir]);
    }
    dV = Kokkos::create_mirror_view(grid.dV);


}

void GridHost::MakeGrid(Grid &grid, Input &input) {
    std::cout << "GridHost::Makegrid\n";
    // Create the grid
    // TBD

    // Sync with the device
    for(int dir = 0 ; dir < 3 ; dir++) {
        Kokkos::deep_copy(grid.x[dir],x[dir]);
        Kokkos::deep_copy(grid.xr[dir],xr[dir]);
        Kokkos::deep_copy(grid.xl[dir],xl[dir]);
        Kokkos::deep_copy(grid.dx[dir],dx[dir]);

        Kokkos::deep_copy(grid.A[dir],A[dir]);
    }

    Kokkos::deep_copy(grid.dV,dV);
}