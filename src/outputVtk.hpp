#ifndef OUTPUTVTK_HPP
#define OUTPUTVTK_HPP
#include "idefix.hpp"

class OutputVTK {
public:

    OutputVTK(Grid &);                     // Create Output Object
    int Write(DataBlock datain);         // Create a VTK from the current DataBlock

private:
    GridHost grid;
    int vtkFileNumber;

    // dimensions
    int nx1,nx2,nx3;

    // Coordinates needed by VTK outputs
    float **node_coord, *xnode, *ynode, *znode;

    // Endianness swaping function and variables
    float BigEndian(float);
    int doneEndianTest, shouldSwapEndian;
};

#endif