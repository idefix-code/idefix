#ifndef INPUT_HPP
#define INPUT_HPP
#include "idefix.hpp"


class Input {
public:
    int npoints[3];
    real xstart[3];
    real xend[3];
    int nghost[3];
    int nstages;
    real firstDt;
    real tperiodVTK;
    real tfinal;
    Input();
private:
    void PrintLogo();
};

#endif