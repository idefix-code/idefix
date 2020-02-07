#include "idefix.hpp"
#include "data.hpp"

Data::Data() {
    // Do nothing
}
Data::Data(Grid &grid) {
    std::cout << "Building data Array with size " << NVAR << "x" << grid.np_tot[KDIR] << "x" << grid.np_tot[JDIR] << "x" << grid.np_tot[IDIR] << std::endl;

    Vc = IdefixArray4D<real>("Vc", NVAR, grid.np_tot[KDIR], grid.np_tot[JDIR], grid.np_tot[IDIR]);
    Uc = IdefixArray4D<real>("Uc", NVAR, grid.np_tot[KDIR], grid.np_tot[JDIR], grid.np_tot[IDIR]);
}

Data::Data(const Data &data) {
    std::cout << "Data copy constructor";

    Vc=data.Vc;
    Uc=data.Uc;
}

Data& Data::operator=(const Data& data) {
    std::cout << "Assignement operator on Data" << std::endl;
    if (this != & data) {
        Vc=data.Vc;
        Uc=data.Uc;
    }

    return *this;
}