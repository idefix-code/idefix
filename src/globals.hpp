#ifndef GLOBALS_HPP
#define GLOBALS_HPP
#include "idefix.hpp"


class Globals {
public:
    IdefixArray1D<real> inputParam;
    IdefixArray1D<real>::HostMirror inputParamH;
    Globals();
};

#endif