#ifndef SETUP_HPP
#define SETUP_HPP
#include "idefix.hpp"

class Setup {
public:
    Setup();
    Setup(Input &, Grid &, DataBlock &, Hydro &);

    void InitFlow(DataBlock &);
    void MakeAnalysis(DataBlock&, real);

};



#endif