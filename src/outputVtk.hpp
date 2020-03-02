#ifndef OUTPUTVTK_HPP
#define OUTPUTVTK_HPP
#include "idefix.hpp"




class OutputVTK {
public:

    OutputVTK(Grid &);                     // Create Output Object
    int Write(DataBlock &);         // Create a VTK from the current DataBlock

private:
    GridHost grid;
    int vtkFileNumber;

    // dimensions
    long int nx1,nx2,nx3;

    // Coordinates needed by VTK outputs
    float *node_coord, *xnode, *ynode, *znode, *Vwrite;

    // Array designed to store the temporary vector array
    IdefixHostArray3D<float> vect3D;

    // Endianness swaping function and variables
    
    int doneEndianTest, shouldSwapEndian;

    void WriteHeader(FILE *fvtk);
    void WriteScalar(FILE *, IdefixHostArray3D<float> &,  std::string &);
    float BigEndian(float);

};

#endif