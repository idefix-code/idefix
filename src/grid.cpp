#include "idefix.hpp"
#include "grid.hpp"

Grid::Grid(Input &i) {
    for(int dir = 0 ; dir < 3 ; dir++) {
        this->x[dir] = IdefixArray1D<real>("Grid_x",i.npoints[dir]);
        this->xH[dir] = create_mirror(this->x[dir]);

    }
}