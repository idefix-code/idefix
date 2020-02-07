#include "idefix.hpp"
#include "data.hpp"

Data::Data() {
    // Do nothing
}
Data::Data(Grid &grid) {

    Vc = IdefixArray4D<real>("Vc", NVAR, grid.np_tot[KDIR], grid.np_tot[JDIR], grid.np_tot[IDIR]);
    Uc = IdefixArray4D<real>("Uc", NVAR, grid.np_tot[KDIR], grid.np_tot[JDIR], grid.np_tot[IDIR]);

}

Data::Data(const Data &data) {

    this->Vc=data.Vc;
    this->Uc=data.Uc;

}

Data& Data::operator=(const Data& data) {

    if (this != & data) {
        Vc=data.Vc;
        Uc=data.Uc;
    }

    return *this;
}