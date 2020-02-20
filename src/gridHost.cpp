#include "idefix.hpp"
#include "grid.hpp"
#include "gridHost.hpp"

GridHost::GridHost(Grid &grid) {

    for(int dir = 0 ; dir < 3 ; dir++) {

        nghost[dir] = grid.nghost[dir];
        np_tot[dir] = grid.np_tot[dir];
        np_int[dir] = grid.np_int[dir];

        lbound[dir] = grid.lbound[dir]; 
        rbound[dir] = grid.rbound[dir];


    }

    // Create mirrors on host
    for(int dir = 0 ; dir < 3 ; dir++) {
        x[dir] = Kokkos::create_mirror_view(grid.x[dir]);
        xr[dir] = Kokkos::create_mirror_view(grid.xr[dir]);
        xl[dir] = Kokkos::create_mirror_view(grid.xl[dir]);
        dx[dir] = Kokkos::create_mirror_view(grid.dx[dir]);
    }


}

void GridHost::MakeGrid(Grid &grid, Input &input) {

    // Create the grid

    for(int dir = 0 ; dir < 3 ; dir++) {
        for(int i = 0 ; i < np_tot[dir] ; i++) {
            dx[dir](i) = (input.xend[dir]-input.xstart[dir])/(np_int[dir]+1);
            x[dir](i)=input.xstart[dir] + (i-nghost[dir]+HALF_F)*dx[dir](i);
            xl[dir](i)=input.xstart[dir] + (i-nghost[dir])*dx[dir](i);
            xr[dir](i)=input.xstart[dir] + (i-nghost[dir]+1)*dx[dir](i);
        }
    }

    
}

void GridHost::SyncFromDevice() {
    for(int dir = 0 ; dir < 3 ; dir++) {
        Kokkos::deep_copy(x[dir],grid.x[dir]);
        Kokkos::deep_copy(xr[dir],grid.xr[dir]);
        Kokkos::deep_copy(xl[dir],grid.xl[dir]);
        Kokkos::deep_copy(dx[dir],grid.dx[dir]);
    }
}

void GridHost::SyncToDevice() {
    // Sync with the device
    for(int dir = 0 ; dir < 3 ; dir++) {
        Kokkos::deep_copy(grid.x[dir],x[dir]);
        Kokkos::deep_copy(grid.xr[dir],xr[dir]);
        Kokkos::deep_copy(grid.xl[dir],xl[dir]);
        Kokkos::deep_copy(grid.dx[dir],dx[dir]);
    }
}