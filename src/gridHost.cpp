#include "idefix.hpp"
#include "grid.hpp"
#include "gridHost.hpp"

GridHost::GridHost() {
    // Void
}
GridHost::GridHost(Grid &grid) {

    Kokkos::Profiling::pushRegion("GridHost::GridHost(Grid)");
    this->grid=grid;
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

    Kokkos::Profiling::popRegion();
}

void GridHost::MakeGrid(Input &input) {

    Kokkos::Profiling::pushRegion("GridHost::MakeGrid");
    real xstart[3];
    real xend[3];
    // Create the grid

    // Get grid parameters from input file, block [Grid]
    std::cout << "GridHost::MakeGrid" << std::endl;
    for(int dir = 0 ; dir < 3 ; dir++) {
        std::string label = std::string("X")+std::to_string(dir+1)+std::string("-grid");
        int numPatch = input.GetInt("Grid",label,0);

        xstart[dir] = input.GetReal("Grid",label,1);
        xend[dir] = input.GetReal("Grid",label,4);

        if(dir<DIMENSIONS) {

            std::string gridType = input.GetString("Grid",label,3);
            if(gridType.compare("u")) {
                std::stringstream msg;
                msg << "While creating Grid: this version of idefix doesn't handle non-uniform grid." << std::endl;
                msg << "Will assume uniform grid.";
                IDEFIX_WARNING(msg);
            }
            std::string lboundString, rboundString;
            switch(lbound[dir]) {
                case outflow:
                    lboundString="outlow";
                    break;
                case periodic:
                    lboundString="periodic";
                    break;
                case internal:
                    lboundString="internal";
                    break;
                case userdef:
                    lboundString="userdef";
                    break;
                default:
                    lboundString="unknown";
            }
            switch(rbound[dir]) {
                case outflow:
                    rboundString="outlow";
                    break;
                case periodic:
                    rboundString="periodic";
                    break;
                case internal:
                    rboundString="internal";
                    break;
                case userdef:
                    rboundString="userdef";
                    break;
                default:
                    rboundString="unknown";
            }

            std::cout << "\t Direction X" << (dir+1) << ": " << lboundString << "\t" << xstart[dir] << "...." << np_int[dir] << "...." << xend[dir] << "\t" << rboundString << std::endl;
        }

    }


    for(int dir = 0 ; dir < 3 ; dir++) {
        for(int i = 0 ; i < np_tot[dir] ; i++) {
            dx[dir](i) = (xend[dir]-xstart[dir])/(np_int[dir]+1);
            x[dir](i)=xstart[dir] + (i-nghost[dir]+HALF_F)*dx[dir](i);
            xl[dir](i)=xstart[dir] + (i-nghost[dir])*dx[dir](i);
            xr[dir](i)=xstart[dir] + (i-nghost[dir]+1)*dx[dir](i);
        }
    }

    Kokkos::Profiling::popRegion();
}

void GridHost::SyncFromDevice() {
    Kokkos::Profiling::pushRegion("GridHost::SyncFromDevice");
    for(int dir = 0 ; dir < 3 ; dir++) {
        Kokkos::deep_copy(x[dir],grid.x[dir]);
        Kokkos::deep_copy(xr[dir],grid.xr[dir]);
        Kokkos::deep_copy(xl[dir],grid.xl[dir]);
        Kokkos::deep_copy(dx[dir],grid.dx[dir]);
    }
    Kokkos::Profiling::popRegion();
}

void GridHost::SyncToDevice() {
    Kokkos::Profiling::pushRegion("GridHost::SyncToDevice");
    // Sync with the device
    for(int dir = 0 ; dir < 3 ; dir++) {
        Kokkos::deep_copy(grid.x[dir],x[dir]);
        Kokkos::deep_copy(grid.xr[dir],xr[dir]);
        Kokkos::deep_copy(grid.xl[dir],xl[dir]);
        Kokkos::deep_copy(grid.dx[dir],dx[dir]);
    }
    Kokkos::Profiling::popRegion();
}