#include "idefix.hpp"
#include "gridHost.hpp"
#include "grid.hpp"

Grid::Grid() {
    // Do nothing
}

Grid::Grid(Input &input) {
    Kokkos::Profiling::pushRegion("Grid::Grid(Input)");
    int npoints[3];


    std::cout << "Building Grid... " << std::endl;
    
    // Get grid size from input file, block [Grid]
    for(int dir = 0 ; dir < 3 ; dir++) {
        npoints[dir] = 1;
        nghost[dir] = 0;
        std::string label = std::string("X")+std::to_string(dir+1)+std::string("-grid");
        int numPatch = input.GetInt("Grid",label,0);
        if(numPatch>1) {
            std::stringstream msg;
            msg << "While creating Grid: this version of idefix doesn't handle more than one grid patch in each direction" << std::endl;
            msg << "We will build a grid based only on the first patch.";
            IDEFIX_WARNING(msg); 
        }

        if(dir<DIMENSIONS) {
            npoints[dir] = input.GetInt("Grid",label,2);
            nghost[dir] = 2;
        }

    }


    for(int dir = 0 ; dir < 3 ; dir++) {
        np_tot[dir] = npoints[dir] + 2*nghost[dir];
        np_int[dir] = npoints[dir];

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

    Kokkos::Profiling::popRegion();

}

