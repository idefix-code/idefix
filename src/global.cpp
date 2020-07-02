#include "idefix.hpp"

namespace idfx {

    int prank;
    int psize;

    IdefixOstream cout;

#ifdef WITH_MPI
    extern MPI_Comm CartComm;
#endif

    int initialize() {
        #ifdef WITH_MPI
        MPI_Comm_size(MPI_COMM_WORLD,&psize);
        MPI_Comm_rank(MPI_COMM_WORLD,&prank);
        #else
        psize=1;
        prank=0;
        #endif
        cout.init(prank);
        return(0);

    }   // Initialisation routine for idefix

    // Init the iostream with defined rank
    void IdefixOstream::init(int rank) {
        char logFileName[20];
        sprintf(logFileName,"idefix.%d.log",rank);
        this->my_fstream = std::ofstream(logFileName);
        if(rank==0) this->toscreen=true;
        else this->toscreen=false;
    }

}
