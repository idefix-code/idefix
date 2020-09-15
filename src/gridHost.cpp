// ********************************************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ********************************************************************************************************

#include "idefix.hpp"
#include "grid.hpp"
#include "gridHost.hpp"

GridHost::GridHost() {
    // Void
}
GridHost::GridHost(Grid &grid) {

    idfx::pushRegion("GridHost::GridHost(Grid)");
    this->grid=&grid;
    for(int dir = 0 ; dir < 3 ; dir++) {

        nghost[dir] = grid.nghost[dir];
        np_tot[dir] = grid.np_tot[dir];
        np_int[dir] = grid.np_int[dir];

        lbound[dir] = grid.lbound[dir]; 
        rbound[dir] = grid.rbound[dir];

        xbeg[dir] = grid.xbeg[dir];
        xend[dir] = grid.xend[dir];


    }

    // Create mirrors on host
    for(int dir = 0 ; dir < 3 ; dir++) {
        x[dir] = Kokkos::create_mirror_view(grid.x[dir]);
        xr[dir] = Kokkos::create_mirror_view(grid.xr[dir]);
        xl[dir] = Kokkos::create_mirror_view(grid.xl[dir]);
        dx[dir] = Kokkos::create_mirror_view(grid.dx[dir]);
    }

    idfx::popRegion();
}

void GridHost::MakeGrid(Input &input) {

    idfx::pushRegion("GridHost::MakeGrid");
    real xstart[3];
    real xend[3];
    // Create the grid

    // Get grid parameters from input file, block [Grid]
    idfx::cout << "GridHost::MakeGrid" << std::endl;
    for(int dir = 0 ; dir < 3 ; dir++) {
        std::string label = std::string("X")+std::to_string(dir+1)+std::string("-grid");
        int numPatch = input.GetInt("Grid",label,0);

        xstart[dir] = input.GetReal("Grid",label,1);
        xend[dir] = input.GetReal("Grid",label,4);

        this->xbeg[dir] = xstart[dir];
        this->xend[dir] = xend[dir];


        if(dir<DIMENSIONS) {

            std::string gridType = input.GetString("Grid",label,3);
            if(gridType.compare("u")) {
                std::stringstream msg;
                msg << "While creating Grid: this version of idefix doesn't handle non-uniform grid->" << std::endl;
                msg << "Will assume uniform grid->";
                IDEFIX_WARNING(msg);
            }
            std::string lboundString, rboundString;
            switch(lbound[dir]) {
                case outflow:
                    lboundString="outflow";
                    break;
                case periodic:
                    lboundString="periodic";
                    break;
                case internal:
                    lboundString="internal";
                    break;
                case shearingbox:
                    lboundString="shearingbox";
                    break;
                case userdef:
                    lboundString="userdef";
                    break;
                default:
                    lboundString="unknown";
            }
            switch(rbound[dir]) {
                case outflow:
                    rboundString="outflow";
                    break;
                case periodic:
                    rboundString="periodic";
                    break;
                case internal:
                    rboundString="internal";
                    break;
                case shearingbox:
                    rboundString="shearingbox";
                    break;
                case userdef:
                    rboundString="userdef";
                    break;
                default:
                    rboundString="unknown";
            }

            idfx::cout << "\t Direction X" << (dir+1) << ": " << lboundString << "\t" << xstart[dir] << "...." << np_int[dir] << "...." << xend[dir] << "\t" << rboundString << std::endl;
        }

    }


    for(int dir = 0 ; dir < 3 ; dir++) {
        for(int i = 0 ; i < np_tot[dir] ; i++) {
            dx[dir](i) = (xend[dir]-xstart[dir])/(np_int[dir]);
            x[dir](i)=xstart[dir] + (i-nghost[dir]+HALF_F)*dx[dir](i);
            xl[dir](i)=xstart[dir] + (i-nghost[dir])*dx[dir](i);
            xr[dir](i)=xstart[dir] + (i-nghost[dir]+1)*dx[dir](i);
        }
    }

    idfx::popRegion();
}

void GridHost::SyncFromDevice() {
    idfx::pushRegion("GridHost::SyncFromDevice");
    for(int dir = 0 ; dir < 3 ; dir++) {
        Kokkos::deep_copy(x[dir],grid->x[dir]);
        Kokkos::deep_copy(xr[dir],grid->xr[dir]);
        Kokkos::deep_copy(xl[dir],grid->xl[dir]);
        Kokkos::deep_copy(dx[dir],grid->dx[dir]);

        xbeg[dir] = grid->xbeg[dir];
        xend[dir] = grid->xend[dir];
    }
    idfx::popRegion();
}

void GridHost::SyncToDevice() {
    idfx::pushRegion("GridHost::SyncToDevice");
    // Sync with the device
    for(int dir = 0 ; dir < 3 ; dir++) {
        Kokkos::deep_copy(grid->x[dir],x[dir]);
        Kokkos::deep_copy(grid->xr[dir],xr[dir]);
        Kokkos::deep_copy(grid->xl[dir],xl[dir]);
        Kokkos::deep_copy(grid->dx[dir],dx[dir]);

        grid->xbeg[dir] = xbeg[dir];
        grid->xend[dir] = xend[dir];

    }
    idfx::popRegion();
}