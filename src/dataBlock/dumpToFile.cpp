#include "../idefix.hpp"
#include "dataBlock.hpp"

// dump the current dataBlock to a file (mainly used for debug purposes)

void DataBlock::DumpToFile(std::string filebase)  {

    FILE *fileHdl;

    real nfield;
    real nx1;
    real nx2;
    real nx3;
    real nv;

    IdefixArray4D<real>::HostMirror locVc = Kokkos::create_mirror_view(this->Vc);
    Kokkos::deep_copy(locVc,this->Vc);
#if MHD == YES
    IdefixArray4D<real>::HostMirror locVs = Kokkos::create_mirror_view(this->Vs);
    Kokkos::deep_copy(locVs,this->Vs);
#endif

    // Data format:
    // First groupe is nfield (number of 4D arrays written)
    // Then for each field:
    // nx1, nx2, nx3, nvar then the data, then off you go 
    std::string dot = std::string(".");
    std::string ext = std::string("idfx");
    std::string filename = filebase + dot + std::to_string(idfx::prank) + dot + ext;

    fileHdl = fopen(filename.c_str(),"wb");
#if MHD== YES
    nfield = 2;
#else
    nfield = 1;
#endif
    fwrite(&nfield, sizeof(real), 1, fileHdl);

    // Write Vc
    nx1=this->np_tot[IDIR];
    nx2=this->np_tot[JDIR];
    nx3=this->np_tot[KDIR];
    nv=NVAR;

    fwrite(&nx1, sizeof(real),1,fileHdl);
    fwrite(&nx2, sizeof(real),1,fileHdl);
    fwrite(&nx3, sizeof(real),1,fileHdl);
    fwrite(&nv, sizeof(real),1,fileHdl);

    fwrite(locVc.data(), sizeof(real), nx1*nx2*nx3*nv, fileHdl);

    // Write Vs
#if MHD == YES
    nx1=nx1+IOFFSET;
    nx2=nx2+JOFFSET;
    nx3=nx3+KOFFSET;
    nv=DIMENSIONS;

    fwrite(&nx1, sizeof(real),1,fileHdl);
    fwrite(&nx2, sizeof(real),1,fileHdl);
    fwrite(&nx3, sizeof(real),1,fileHdl);
    fwrite(&nv, sizeof(real),1,fileHdl);

    fwrite(locVs.data(), sizeof(real), nx1*nx2*nx3*nv, fileHdl);
#endif

    fclose(fileHdl);


}