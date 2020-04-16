#ifndef ElectromotiveForce_HPP
#define ElectromotiveForce_HPP
#include "idefix.hpp"

class DataBlock;

class ElectroMotiveForce {
    public:
        /* Face centered emf components */
        IdefixArray3D<real>     exj;
        IdefixArray3D<real>     exk;
        IdefixArray3D<real>     eyi;
        IdefixArray3D<real>     eyk;
        IdefixArray3D<real>     ezi;
        IdefixArray3D<real>     ezj;

        /* Edge centered emf components */
        IdefixArray3D<real>     ex;
        IdefixArray3D<real>     ey;
        IdefixArray3D<real>     ez;

        /* Range of existence */
        int  ibeg, jbeg, kbeg;
        int  iend, jend, kend;

        // Constructor from datablock structure
        ElectroMotiveForce(DataBlock *);

        // Default constructor
        ElectroMotiveForce();

};

#endif



    