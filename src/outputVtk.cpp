#include "outputVtk.hpp"
#include "gitversion.h"

#define VTK_RECTILINEAR_GRID    14
#define VTK_STRUCTURED_GRID     35

#ifndef VTK_FORMAT
  #if GEOMETRY == CARTESIAN || GEOMETRY == CYLINDRICAL
    #define VTK_FORMAT  VTK_RECTILINEAR_GRID
  #else
    #define VTK_FORMAT  VTK_STRUCTURED_GRID
  #endif
#endif

#define WRITE_STAGGERED_FIELD
#define WRITE_EMF

/* ---------------------------------------------------------
    The following macros are specific to this file only 
    and are used to ease up serial/parallel implementation
    for writing strings and real arrays 
   --------------------------------------------------------- */   
    
#ifdef PARALLEL
 #define VTK_HEADER_WRITE_STRING(header) \
         TOBEDEFINED(header, strlen(header), MPI_CHAR, SZ_Float_Vect);
 #define VTK_HEADER_WRITE_FLTARR(arr, nelem) \
         TOBEDEFINED (arr, nelem, MPI_FLOAT, SZ_Float_Vect);
 #define VTK_HEADER_WRITE_DBLARR(arr, nelem) \
         TOBEDEFINED (arr, nelem, MPI_DOUBLE, SZ_Float_Vect);
#else
 #define VTK_HEADER_WRITE_STRING(header) \
         fprintf (fvtk,header);
 #define VTK_HEADER_WRITE_FLTARR(arr,nelem) \
         fwrite(arr, sizeof(float), nelem, fvtk);
 #define VTK_HEADER_WRITE_DBLARR(arr,nelem) \
         fwrite(arr, sizeof(double), nelem, fvtk);
#endif


/* Main constructor */
OutputVTK::OutputVTK(Input &input, DataBlock &datain, real t)
{
    // Init the output period
    this->tperiod=input.GetReal("Output","vtk",0);
    this->tnext = t;

    // Initialize the output structure
    // Create a local datablock as an image of gridin
    DataBlockHost data(datain);
    data.SyncFromDevice();


    // Create the coordinate array required in VTK files
    this->nx1 = data.np_int[IDIR];
    this->nx2 = data.np_int[JDIR];
    this->nx3 = data.np_int[KDIR];

    // Vector array where we store the pencil before write
    this->Vwrite = new float[nx1+IOFFSET];

    // Temporary storage on host for 3D arrays
    this->vect3D = IdefixHostArray3D<real>("vect3D",data.np_tot[KDIR],data.np_tot[JDIR],data.np_tot[IDIR]);

    // Essentially does nothing
    this->vtkFileNumber = 0;

    // Test endianness
    int tmp1 = 1;
    this->shouldSwapEndian = 0;
    unsigned char *tmp2 = (unsigned char *) &tmp1;
    if (*tmp2 != 0)
        this->shouldSwapEndian = 1;
    
    
#if VTK_FORMAT == VTK_RECTILINEAR_GRID
    this->xnode = new float[nx1+IOFFSET];
    this->ynode = new float[nx2+JOFFSET];
    this->znode = new float[nx3+KOFFSET];

    for (long int i = 0; i < nx1 + IOFFSET; i++) {
        xnode[i] = BigEndian(data.xl[IDIR](i + data.nghost[IDIR]));
    }
    for (long int j = 0; j < nx2 + JOFFSET; j++)    {
        ynode[j] = BigEndian(data.xl[JDIR](j + data.nghost[JDIR]));
    }
    for (long int k = 0; k < nx3 + KOFFSET; k++)
    {
        if(DIMENSIONS==2) znode[k] = BigEndian(0.0);
        else znode[k] = BigEndian(data.xl[KDIR](k + data.nghost[KDIR]));
    }
#else   // VTK_FORMAT
        /* -- Allocate memory for node_coord which is later used -- */
    node_coord = new float[(nx1+IOFFSET)*3]; 
    
#endif

}

int OutputVTK::Write(DataBlock &datain, real t)
{
    FILE *fileHdl;
    char filename[256];

    // Do we need an output?
    if(t<this->tnext) return(0);

    this->tnext+= this->tperiod;
    Kokkos::Profiling::pushRegion("OutputVTK::Write");

    idfx::cout << "OutputVTK::Write file n " << vtkFileNumber << "..." << std::flush;

    timer.reset();

    // Create a copy of the dataBlock on Host, and sync it.
    DataBlockHost data(datain);
    data.SyncFromDevice();
    if(idfx::psize>1) sprintf (filename, "data.%04d.%04d.vtk", idfx::prank, vtkFileNumber);
    else sprintf (filename, "data.%04d.vtk", vtkFileNumber);
    fileHdl = fopen(filename,"wb");
    WriteHeader(fileHdl);
    for(int nv = 0 ; nv < NVAR ; nv++) {
        for(int k = 0; k < data.np_tot[KDIR] ; k++ ) {
            for(int j = 0; j < data.np_tot[JDIR] ; j++ ) {
                for(int i = 0; i < data.np_tot[IDIR] ; i++ ) {
                    vect3D(k,j,i) = data.Vc(nv,k,j,i);
                }
            }
        }
        WriteScalar(fileHdl, vect3D, datain.VcName[nv],datain);
    }

#if MHD == YES
#ifdef WRITE_STAGGERED_FIELD
    for(int nv = 0 ; nv < DIMENSIONS ; nv++) {
        for(int k = 0; k < data.np_tot[KDIR] ; k++ ) {
            for(int j = 0; j < data.np_tot[JDIR] ; j++ ) {
                for(int i = 0; i < data.np_tot[IDIR] ; i++ ) {
                    vect3D(k,j,i) = data.Vs(nv,k,j,i);
                }
            }
        }
        WriteScalar(fileHdl, vect3D, datain.VsName[nv],datain);
    }
#endif // WRITE_STAGGERED_FIELD

#ifdef WRITE_EMF
    std::string varname;
#if DIMENSIONS == 3
    Kokkos::deep_copy(vect3D,datain.emf.ex);
    varname="Ex";
    WriteScalar(fileHdl, vect3D, varname,datain);

    Kokkos::deep_copy(vect3D,datain.emf.ey);
    varname="Ey";
    WriteScalar(fileHdl, vect3D, varname,datain);
#endif // DIMENSIONS

    Kokkos::deep_copy(vect3D,datain.emf.ez);
    varname="Ez";
    WriteScalar(fileHdl, vect3D, varname,datain);
#endif // WRITE_EMF
#endif// MHD


    fclose(fileHdl);

    vtkFileNumber++;
    // Make file number

    idfx::cout << "done in " << timer.seconds() << " s." << std::endl;
    Kokkos::Profiling::popRegion();
    // One day, we will have a return code.
    return(0);
}



/* ********************************************************************* */
void OutputVTK::WriteHeader(FILE *fvtk)
/*!
 * Write VTK header in parallel or serial mode.
 *
 * \param [in]  fvtk  pointer to file
 * \param [in]  grid  pointer to an array of Grid structures
 *
 * \todo  Write the grid using several processors.
 *********************************************************************** */
{
    char header[1024];
    float x1, x2, x3;


    /* -------------------------------------------
   1. File version and identifier
   ------------------------------------------- */

    sprintf(header, "# vtk DataFile Version 2.0\n");

    /* -------------------------------------------
   2. Header
   ------------------------------------------- */

    sprintf(header + strlen(header), "Idefix %s VTK Data\n", GITVERSION);

    /* ------------------------------------------
   3. File format
   ------------------------------------------ */

    sprintf(header + strlen(header), "BINARY\n");

    /* ------------------------------------------
   4. Dataset structure
   ------------------------------------------ */

#if VTK_FORMAT == VTK_RECTILINEAR_GRID
    sprintf(header + strlen(header), "DATASET %s\n", "RECTILINEAR_GRID");
#elif VTK_FORMAT == VTK_STRUCTURED_GRID
    sprintf(header + strlen(header), "DATASET %s\n", "STRUCTURED_GRID");
#endif

    
    VTK_HEADER_WRITE_STRING(header);

    /* -- Generate time info (VisIt reader only) -- */

    /*
#if VTK_TIME_INFO == YES
  sprintf (header,"FIELD FieldData 1\n");
  sprintf (header+strlen(header),"TIME 1 1 double\n");
  double tt=g_time;
  if (IsLittleEndian()) SWAP_VAR(tt);
  VTK_HEADER_WRITE_STRING(header);
  VTK_HEADER_WRITE_DBLARR(&tt, 1);
  VTK_HEADER_WRITE_STRING("\n");
#endif /* VTK_TIME_INFO */

    sprintf(header, "DIMENSIONS %d %d %d\n",
            nx1 + IOFFSET, nx2 + JOFFSET, nx3 + KOFFSET);
    VTK_HEADER_WRITE_STRING(header);

#if VTK_FORMAT == VTK_RECTILINEAR_GRID

    /* -- Write rectilinear grid information -- */

    sprintf(header, "X_COORDINATES %d float\n", nx1 + IOFFSET);
    VTK_HEADER_WRITE_STRING(header);
    VTK_HEADER_WRITE_FLTARR(xnode, nx1 + IOFFSET);

    sprintf(header, "\nY_COORDINATES %d float\n", nx2 + JOFFSET);
    VTK_HEADER_WRITE_STRING(header);
    VTK_HEADER_WRITE_FLTARR(ynode, nx2 + JOFFSET);

    sprintf(header, "\nZ_COORDINATES %d float\n", nx3 + KOFFSET);
    VTK_HEADER_WRITE_STRING(header);
    VTK_HEADER_WRITE_FLTARR(znode, nx3 + KOFFSET);

#elif VTK_FORMAT == VTK_STRUCTURED_GRID

    /* -- define node_coord -- */

    sprintf(header, "POINTS %d float\n", (nx1 + IOFFSET) * (nx2 + JOFFSET) * (nx3 + KOFFSET));
    VTK_HEADER_WRITE_STRING(header);

    /* -- Write structured grid information -- */

    x1 = x2 = x3 = 0.0;
    for (long int k = 0; k < nx3 + KOFFSET; k++)
    {
        for (long int j = 0; j < nx2 + JOFFSET; j++)
        {
            for (long int i = 0; i < nx1 + IOFFSET; i++)
            {
                D_EXPAND(x1 = data.xl[IDIR](i + data.nghost[IDIR]);,
                         x2 = data.xl[JDIR](j + data.nghost[JDIR]);,
                         x3 = data.xl[KDIR](k + data.nghost[KDIR]);)

#if (GEOMETRY == CARTESIAN) || (GEOMETRY == CYLINDRICAL)
                node_coord[3*i+IDIR] = BigEndian(x1);
                node_coord[3*i+JDIR] = BigEndian(x2);
                node_coord[3*i+KDIR] = BigEndian(x3);
#elif GEOMETRY == POLAR
                node_coord[3*i+IDIR] = BigEndian(x1 * cos(x2));
                node_coord[3*i+JDIR] = BigEndian(x1 * sin(x2));
                node_coord[3*i+KDIR] = BigEndian(x3);
#elif GEOMETRY == SPHERICAL
#if DIMENSIONS == 2
                node_coord[3*i+IDIR] = BigEndian(x1 * sin(x2));
                node_coord[3*i+JDIR] = BigEndian(x1 * cos(x2));
                node_coord[3*i+KDIR] = BigEndian(0.0);
#elif DIMENSIONS == 3
                node_coord[3*i+IDIR] = BigEndian(x1 * sin(x2) * cos(x3));
                node_coord[3*i+JDIR] = BigEndian(x1 * sin(x2) * sin(x3));
                node_coord[3*i+KDIR] = BigEndian(x1 * cos(x2));
#endif
#endif
            }
            VTK_HEADER_WRITE_FLTARR(node_coord, 3 * (nx1 + IOFFSET));
        }
    }

#endif

    /* -----------------------------------------------------
   5. Dataset attributes [will continue by later calls
      to WriteVTK_Vector() or WriteVTK_Scalar()...]
   ----------------------------------------------------- */

    sprintf(header, "\nCELL_DATA %d\n", nx1 * nx2 * nx3);
    VTK_HEADER_WRITE_STRING(header);
}
#undef VTK_STRUCTERED_GRID
#undef VTK_RECTILINEAR_GRID



/* ********************************************************************* */
void OutputVTK::WriteScalar(FILE *fvtk, IdefixHostArray3D<real> &Vin,  std::string &var_name, DataBlock& data)
/*!
 * Write VTK scalar field.
 *
 * \param [in]   fvtk       pointer to file (handle)
 * \param [in]   V          pointer to 3D data array
 * \param [in] unit     the corresponding cgs unit (if specified, 1 otherwise)
 * \param [in]   var_name   the variable name appearing in the VTK file
 * \param [in]   grid       pointer to an array of Grid structures
 *********************************************************************** */
{
    int i, j, k;
    char header[512];


    sprintf(header, "\nSCALARS %s float\n", var_name.c_str());
    sprintf(header + strlen(header), "LOOKUP_TABLE default\n");


    fprintf(fvtk, "%s", header);


    for(long int k = 0 ; k < nx3 ; k++ ) {
        for(long int j = 0 ; j < nx2 ; j++ ) {
            for(long int i = 0 ; i < nx1 ; i++ ) {
                Vwrite[i] = BigEndian(float(Vin(k + data.nghost[KDIR],j + data.nghost[JDIR],i + data.nghost[IDIR])));
            }
            fwrite(Vwrite, sizeof(float), nx1, fvtk);
        }
    }
}

/* ****************************************************************************/
/** Determines if the machine is little-endian.  If so, 
    it will force the data to be big-endian. 
	@param in_number floating point number to be converted in big endian */
/* *************************************************************************** */

float OutputVTK::BigEndian(float in_number)
{
    if (shouldSwapEndian)
    {
		unsigned char *bytes = (unsigned char*) &in_number;
        unsigned char tmp = bytes[0];
        bytes[0] = bytes[3];
        bytes[3] = tmp;
        tmp = bytes[1];
        bytes[1] = bytes[2];
        bytes[2] = tmp;
    }
	return(in_number);
}