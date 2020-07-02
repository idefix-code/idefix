#include "idefix.hpp"
#include "gridHost.hpp"
#include "grid.hpp"

Grid::Grid() {
    // Do nothing
}

Grid::Grid(Input &input) {
    idfx::pushRegion("Grid::Grid(Input)");
    int npoints[3];


    idfx::cout << "Building Grid... " << std::endl;
    
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
        
        std::string label = std::string("X")+std::to_string(dir+1)+std::string("-beg");
        std::string boundary = input.GetString("Boundary",label,0);
        if(boundary.compare("outflow") == 0) lbound[dir] = outflow;
        else if(boundary.compare("periodic") == 0) lbound[dir] = periodic;
        else if(boundary.compare("internal") == 0) lbound[dir] = internal;
        else if(boundary.compare("shearingbox") == 0) lbound[dir] = shearingbox;
        else if(boundary.compare("userdef") == 0) lbound[dir] = userdef;
        else {
            std::stringstream msg;
            msg << "Unknown boundary type " << boundary;
            IDEFIX_ERROR(msg);
        }

        label = std::string("X")+std::to_string(dir+1)+std::string("-end");
        boundary = input.GetString("Boundary",label,0);
        if(boundary.compare("outflow") == 0) rbound[dir] = outflow;
        else if(boundary.compare("periodic") == 0) rbound[dir] = periodic;
        else if(boundary.compare("internal") == 0) rbound[dir] = internal;
        else if(boundary.compare("shearingbox") == 0) rbound[dir] = shearingbox;
        else if(boundary.compare("userdef") == 0) rbound[dir] = userdef;
        else {
            std::stringstream msg;
            msg << "Unknown boundary type " << boundary;
            IDEFIX_ERROR(msg);
        }

    }

    // Allocate the grid structure on device. Initialisation will come from GridHost
    
    for(int dir = 0 ; dir < 3 ; dir++) {
        x[dir] = IdefixArray1D<real>("Grid_x",np_tot[dir]);
        xr[dir] = IdefixArray1D<real>("Grid_xr",np_tot[dir]);
        xl[dir] = IdefixArray1D<real>("Grid_xl",np_tot[dir]);
        dx[dir] = IdefixArray1D<real>("Grid_dx",np_tot[dir]);
        
    }

    // Allocate proc structure (by default: one proc in each direction, size one)
    for(int i=0 ; i < 3; i++) {
        nproc[i] = 1;
        xproc[i] = 0;
    }

#ifdef WITH_MPI 
    // Domain decomposition required for the grid

    // Init periodicity array
    int period[3];
    for(int i=0 ; i < 3 ; i++) period[i] = 0;

    // Check if the dec option has been passed when number of procs > 1
    if(idfx::psize>1) {

        if(input.CheckEntry("CommandLine","dec")  != DIMENSIONS) {
            IDEFIX_ERROR("When using MPI with np > 1, -dec option is mandatory.");
        }
        int ntot=1;
        for(int dir=0 ; dir < DIMENSIONS; dir++) {
            nproc[dir] = input.GetInt("CommandLine","dec",dir);
            ntot = ntot * nproc[dir];
        }
        if(ntot != idfx::psize) {
            std::stringstream msg;
            msg << "The number of MPI process (" << idfx::psize << ") does not match your domain decomposition (";
            for(int dir=0 ; dir < DIMENSIONS; dir++) {
                msg << nproc[dir];
                if(dir<DIMENSIONS-1) msg << ", ";
            }
            msg << ").";
            IDEFIX_ERROR(msg);
        }
    }

    // Add periodicity indications
    for(int dir=0 ; dir < DIMENSIONS; dir++) {
        if(rbound[dir] == periodic || rbound[dir] == shearingbox) period[dir] = 1;
    }

    // Create cartesian communicator along with cartesian coordinates.
    MPI_Cart_create(MPI_COMM_WORLD, 3, nproc, period, 0, &CartComm);
    MPI_Cart_coords(CartComm, idfx::prank, 3, xproc);
    
    MPI_Barrier(MPI_COMM_WORLD);
    idfx::cout << "Proc dimensions (";

    for(int dir = 0; dir < 3; dir++) {
        idfx::cout << xproc[dir];
        if(dir < 2) idfx::cout << ", ";
    }
    idfx::cout << ")" << std::endl;

#endif


    idfx::popRegion();

}

/*
Grid& Grid::operator=(const Grid& grid) {
    for(int dir = 0 ; dir < 3 ; dir++) {
        x[dir] = grid.x[dir];
        xr[dir] = grid.xr[dir];
        xl[dir] = grid.xl[dir];
        dx[dir] = grid.dx[dir];
        xbeg[dir] = grid.xbeg[dir];
        xend[dir] = grid.xend[dir];
        np_tot[dir] = grid.np_tot[dir];
        np_int[dir] = grid.np_int[dir];
        nghost[dir] = grid.nghost[dir];
        lbound[dir] = grid.lbound[dir];
        rbound[dir] = grid.rbound[dir];
    }

    return *this;
}
*/

