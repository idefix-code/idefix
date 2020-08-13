#ifndef OUTPUTDMP_HPP
#define OUTPUTDMP_HPP
#include "../idefix.hpp"


enum DataType {DoubleType, SingleType, IntegerType};

class OutputDump {
public:

    OutputDump(Input &, DataBlock &, real);                             // Create Output Object
    int Write(Grid&, DataBlock &, TimeIntegrator&, OutputVTK&);         // Create a Dump file from the current state of the code
    int Read(Grid&, DataBlock &, TimeIntegrator&, OutputVTK&, int);     // Read and load a dump file as current state of the code

private:
    int dumpFileNumber;
    real tperiod, tnext;

    real *scrch;                                                      // Scratch array in host space


    // Timer
    Kokkos::Timer timer;

    // File offset
#ifdef WITH_MPI
    MPI_Offset offset;
    MPI_Datatype view;
#endif

    // len of field names
    const int nameSize  =  16;


    void WriteString(IdfxFileHandler, char *);
    void WriteSerial(IdfxFileHandler, int, int *, DataType, char*, void*);
    void WriteDistributed(IdfxFileHandler, int, int*, char*, real*);
    void ReadNextFieldProperties(IdfxFileHandler, int&, int*, DataType&, std::string&);
    void ReadSerial(IdfxFileHandler, int, int*, DataType, void*);
    void ReadDistributed(IdfxFileHandler, int, int*, void*);
};

#endif