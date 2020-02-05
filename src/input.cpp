#include "idefix.hpp"
#include "input.hpp"

Input::Input() {
    std::cout << "Creating Input\n";
    // Default constructor
    npoints[IDIR] = 64;
    npoints[JDIR] = 64;
    npoints[KDIR] = 64;

    xstart[IDIR] = 0.0;
    xstart[JDIR] = 0.0;
    xstart[KDIR] = 0.0;

    xend[IDIR] = 1.0;
    xend[JDIR] = 1.0;
    xend[KDIR] = 1.0;

    nghost[IDIR] = 2;
    nghost[JDIR] = 2;
    nghost[KDIR] = 2;
}