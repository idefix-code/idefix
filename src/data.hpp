#ifndef DATA_HPP
#define DATA_HPP
#include "idefix.hpp"

class Data {
public:
    IdefixArray4D<real> Vc;     // Main cell-centered primitive variables index

    IdefixArray4D<real> Uc;     // Main cell-centered conservative variables

    // Constructor
    Data(Grid &);

    // Copy constructor
    Data(const Data &);

    // Assignement operator
    Data& operator=(const Data&);


    Data();
};



#endif