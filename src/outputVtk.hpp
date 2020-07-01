#ifndef OUTPUTVTK_HPP
#define OUTPUTVTK_HPP
#include "idefix.hpp"


// File handler depends on the type of I/O we use
#ifdef WITH_MPI
using IdfxFileHandler = MPI_File;
#else
using IdfxFileHandler = FILE*;
#endif

class OutputVTK {
public:

    OutputVTK(Input &, DataBlock &, real);                     // Create Output Object
    int Write(DataBlock &, real);         // Create a VTK from the current DataBlock

private:
    int vtkFileNumber;
    real tperiod, tnext;

    // dimensions
    long int nx1,nx2,nx3;
    long int nx1loc,nx2loc,nx3loc;

    // number of ghost zones
    long int ngx1,ngx2,ngx3;

    // Coordinates needed by VTK outputs
    float *node_coord, *xnode, *ynode, *znode, *Vwrite;

    // Array designed to store the temporary vector array
    float *vect3D;

    // Endianness swaping function and variable
    int doneEndianTest, shouldSwapEndian;

    // Timer
    Kokkos::Timer timer;

    // File offset
#ifdef WITH_MPI
    MPI_Offset offset;
    MPI_Datatype view;
#endif

    void WriteHeader(IdfxFileHandler);
    void WriteScalar(IdfxFileHandler, float*,  std::string &);
    float BigEndian(float);
    void WriteHeaderString(char* , IdfxFileHandler );
    void WriteHeaderFloat(float* , long int, IdfxFileHandler);

};

#endif