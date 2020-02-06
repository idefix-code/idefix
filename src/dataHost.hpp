#ifndef DATAHOST_HPP
#define DATAHOST_HPP
#include "idefix.hpp"

class DataHost {
public:
    IdefixArray4D<real>::HostMirror Vc;     // Main cell-centered primitive variables index

    IdefixArray4D<real>::HostMirror Uc;     // Main cell-centered conservative variables

    // Constructor
    DataHost(Data &);

    // Default constructor
    DataHost();

    // Synchronisation routines
    void SyncToDevice();
    void SyncFromDevice();


private:
    // Data object to which we are the mirror
    Data data;
};



#endif