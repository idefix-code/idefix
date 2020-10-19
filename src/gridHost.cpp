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
        xend[dir] = input.GetReal("Grid",label,4+(numPatch-1)*3);

        this->xbeg[dir] = xstart[dir];
        this->xend[dir] = xend[dir];


        if(dir<DIMENSIONS) {
            // Loop on all the patches
            int idxstart = nghost[dir];
            for(int patch = 0 ; patch < numPatch ; patch++) {
                std::string patchType = input.GetString("Grid",label,3+patch*3);
                real patchStart = input.GetReal("Grid",label,1+patch*3);
                real patchEnd = input.GetReal("Grid",label,4+patch*3);
                int patchSize = input.GetInt("Grid",label,2+patch*3);

                // If this is the first or last patch, also define ghost cells
                int ghostStart = 0;
                int ghostEnd = 0;
                if(patch == 0) ghostStart = nghost[dir];
                if(patch == numPatch-1) ghostEnd = nghost[dir];
                // Define the grid depending on patch type
                if(patchType.compare("u")==0) {
                    // Uniform patch
                    for(int i = 0 - ghostStart ; i < patchSize + ghostEnd ; i++) {
                        dx[dir](i+idxstart) = (patchEnd-patchStart)/(patchSize);
                        x[dir](i+idxstart)=patchStart + (i+HALF_F)*dx[dir](i+idxstart);
                        xl[dir](i+idxstart)=patchStart  + i*dx[dir](i+idxstart);
                        xr[dir](i+idxstart)=patchStart  + (i+1)*dx[dir](i+idxstart);
                    }
                }
                else if(patchType.compare("l")==0) {
                    // log-increasing patch
                    double alpha = (patchEnd + fabs(patchStart) - patchStart)/fabs(patchStart);

                    for(int i = 0 - ghostStart ; i < patchSize + ghostEnd ; i++) {
                        xl[dir](i+idxstart) = patchStart * pow(alpha, (double) i / (double(patchSize)));
                        xr[dir](i+idxstart) = patchStart * pow(alpha, (double) (i+1) / (double(patchSize)));
                        dx[dir](i+idxstart) = xr[dir](i+idxstart) - xl[dir](i+idxstart);
                        x[dir](i+idxstart)= 0.5*(xr[dir](i+idxstart) + xl[dir](i+idxstart));
                    }
                }
                else {
                    std::stringstream msg;
                    msg << "GridHost::MakeGrid: Unknown grid type :" << patchType << std::endl;
                    IDEFIX_ERROR(msg); 
                }

                // Increment offset
                idxstart += patchSize;

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

        else { // dir >= DIMENSIONS/ Init simple uniform grid            
            for(int i = 0 ; i < np_tot[dir] ; i++) {
                dx[dir](i) = (xend[dir]-xstart[dir])/(np_int[dir]);
                x[dir](i)=xstart[dir] + (i-nghost[dir]+HALF_F)*dx[dir](i);
                xl[dir](i)=xstart[dir] + (i-nghost[dir])*dx[dir](i);
                xr[dir](i)=xstart[dir] + (i-nghost[dir]+1)*dx[dir](i);
            }
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