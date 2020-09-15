// ********************************************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ********************************************************************************************************

#ifndef OUTPUTVTK_HPP
#define OUTPUTVTK_HPP
#include "../idefix.hpp"


class OutputVTK {
    friend class OutputDump;
    
public:

    OutputVTK(Input &, DataBlock &, real);                     // Create Output Object
    int Write(DataBlock &, real);         // Create a VTK from the current DataBlock
    int CheckForWrite(DataBlock &, real);

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