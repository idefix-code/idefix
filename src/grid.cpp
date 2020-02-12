#include "idefix.hpp"
#include "gridHost.hpp"
#include "grid.hpp"

Grid::Grid() {
    // Do nothing
}

Grid::Grid(Input &input) {
    std::cout << "Building grid \n";
    // Compute the number of grid points in each direction

    for(int dir = 0 ; dir < 3 ; dir++) {
        nghost[dir] = input.nghost[dir];
        np_tot[dir] = input.npoints[dir] + 2*nghost[dir];
        np_int[dir] = input.npoints[dir];

        // Boundary conditions are yet to be defined
        lbound[dir] = 0; 
        rbound[dir] = 0;

    }

    // Allocate the grid structure on device. Initialisation will come from GridHost
    
    for(int dir = 0 ; dir < 3 ; dir++) {
        x[dir] = IdefixArray1D<real>("Grid_x",np_tot[dir]);
        xr[dir] = IdefixArray1D<real>("Grid_xr",np_tot[dir]);
        xl[dir] = IdefixArray1D<real>("Grid_xl",np_tot[dir]);
        dx[dir] = IdefixArray1D<real>("Grid_dx",np_tot[dir]);
        
    }

}

