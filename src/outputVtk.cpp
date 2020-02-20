#include "outputVtk.hpp"

OutputVTK::OutputVTK(Grid &gridin)
{
    // Initialize the output structure
    // Create a local gridhost as an image of gridin
    grid = GridHost(gridin);
    grid.SyncFromDevice();

    // Essentially does nothing
    vtkFileNumber = 0;

    // Test endianness
    int tmp1 = 1;
    shouldSwapEndian = 0;
    unsigned char *tmp2 = (unsigned char *) &tmp1;
    if (*tmp2 != 0)
        shouldSwapEndian = 1;
    
    // Create the coordinate array required in VTK files
    nx1 = grid.np_tot[IDIR] - 2 * grid.nghost[IDIR];
    nx2 = grid.np_tot[JDIR] - 2 * grid.nghost[JDIR];
    nx3 = grid.np_tot[KDIR] - 2 * grid.nghost[KDIR];

#if VTK_FORMAT == VTK_RECTILINEAR_GRID
    xnode = new float[nx1+IOFFSET];
    ynode = new float[nx2+JOFFSET];
    znode = new float[nx3+KOFFSET];

    for (int i = 0; i < nx1 + IOFFSET; i++) {
        xnode[i] = BigEndian(grid.xl[IDIR](i + grid.nghost[IDIR]));
    }
    for (int j = 0; j < nx2 + JOFFSET; j++)    {
        ynode[j] = BigEndian(grid.xl[JDIR](j + grid.nghost[JDIR]));
    }
    for (int k = 0; k < nx3 + KOFFSET; k++)
    {
        if(DIMENSIONS==2) znode[k] = BigEndian(0.0);
        else znode[k] = BigEndian(grid.xl[KDIR](k + grid.nghost[KDIR]));
    }
#else
        /* -- Allocate memory and define node_coord -- */

    if (node_coord == NULL)
        node_coord = ARRAY_2D(nx1 + IOFFSET, 3, float);

    sprintf(header, "POINTS %d float\n", (nx1 + IOFFSET) * (nx2 + JOFFSET) * (nx3 + KOFFSET));
    VTK_HEADER_WRITE_STRING(header);

    /* -- Write structured grid information -- */

    x1 = x2 = x3 = 0.0;
    for (k = 0; k < nx3 + KOFFSET; k++)
    {
        for (j = 0; j < nx2 + JOFFSET; j++)
        {
            for (i = 0; i < nx1 + IOFFSET; i++)
            {
                D_EXPAND(x1 = grid->xl_glob[IDIR][IBEG + i];,
                                                            x2 = grid->xl_glob[JDIR][JBEG + j];
                         ,
                         x3 = grid->xl_glob[KDIR][KBEG + k];)

#if (GEOMETRY == CARTESIAN) || (GEOMETRY == CYLINDRICAL)
                node_coord[i][0] = x1;
                node_coord[i][1] = x2;
                node_coord[i][2] = x3;
#elif GEOMETRY == POLAR
                node_coord[i][0] = x1 * cos(x2);
                node_coord[i][1] = x1 * sin(x2);
                node_coord[i][2] = x3;
#elif GEOMETRY == SPHERICAL
#if DIMENSIONS == 2
                node_coord[i][0] = x1 * sin(x2);
                node_coord[i][1] = x1 * cos(x2);
                node_coord[i][2] = 0.0;
#elif DIMENSIONS == 3
                node_coord[i][0] = x1 * sin(x2) * cos(x3);
                node_coord[i][1] = x1 * sin(x2) * sin(x3);
                node_coord[i][2] = x1 * cos(x2);
#endif
#endif

                if (IsLittleEndian())
                {
                    SWAP_VAR(node_coord[i][0]);
                    SWAP_VAR(node_coord[i][1]);
                    SWAP_VAR(node_coord[i][2]);
                }
            }
            VTK_HEADER_WRITE_FLTARR(node_coord[0], 3 * (nx1 + IOFFSET));
        }
    }
#endif

}

int OutputVTK::Write(Grid &grid, DataBlock &data)
{
}



/* ********************************************************************* */
void OutputVTK::WriteHeader(FILE *fvtk, GridHost &grid)
/*!
 * Write VTK header in parallel or serial mode.
 *
 * \param [in]  fvtk  pointer to file
 * \param [in]  grid  pointer to an array of Grid structures
 *
 * \todo  Write the grid using several processors.
 *********************************************************************** */
{
    long int i, j, k;
    long int nx1, nx2, nx3;
    char header[1024];
    float x1, x2, x3;

    /* -- Get global domain sizes -- */

    

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

    /* -- Allocate memory and define node_coord -- */

    if (node_coord == NULL)
        node_coord = ARRAY_2D(nx1 + IOFFSET, 3, float);

    sprintf(header, "POINTS %d float\n", (nx1 + IOFFSET) * (nx2 + JOFFSET) * (nx3 + KOFFSET));
    VTK_HEADER_WRITE_STRING(header);

    /* -- Write structured grid information -- */

    x1 = x2 = x3 = 0.0;
    for (k = 0; k < nx3 + KOFFSET; k++)
    {
        for (j = 0; j < nx2 + JOFFSET; j++)
        {
            for (i = 0; i < nx1 + IOFFSET; i++)
            {
                D_EXPAND(x1 = grid->xl_glob[IDIR][IBEG + i];,
                                                            x2 = grid->xl_glob[JDIR][JBEG + j];
                         ,
                         x3 = grid->xl_glob[KDIR][KBEG + k];)

#if (GEOMETRY == CARTESIAN) || (GEOMETRY == CYLINDRICAL)
                node_coord[i][0] = x1;
                node_coord[i][1] = x2;
                node_coord[i][2] = x3;
#elif GEOMETRY == POLAR
                node_coord[i][0] = x1 * cos(x2);
                node_coord[i][1] = x1 * sin(x2);
                node_coord[i][2] = x3;
#elif GEOMETRY == SPHERICAL
#if DIMENSIONS == 2
                node_coord[i][0] = x1 * sin(x2);
                node_coord[i][1] = x1 * cos(x2);
                node_coord[i][2] = 0.0;
#elif DIMENSIONS == 3
                node_coord[i][0] = x1 * sin(x2) * cos(x3);
                node_coord[i][1] = x1 * sin(x2) * sin(x3);
                node_coord[i][2] = x1 * cos(x2);
#endif
#endif

                if (IsLittleEndian())
                {
                    SWAP_VAR(node_coord[i][0]);
                    SWAP_VAR(node_coord[i][1]);
                    SWAP_VAR(node_coord[i][2]);
                }
            }
            VTK_HEADER_WRITE_FLTARR(node_coord[0], 3 * (nx1 + IOFFSET));
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
void WriteVTK_Vector(FILE *fvtk, Data_Arr V, double unit,
                     char *var_name, Grid *grid)
/*!
 * Write VTK vector field data.
 * This is enabled only when VTK_VECTOR_DUMP is set to \c YES.
 * For generality purposes, vectors are written always with 3
 * components, even when there're only 2 being used.
 *
 * The following Maple script has been used to find vector
 * components from cyl/sph to cartesian:
 *
 * \code
   restart;
   with(linalg);
   Acyl := matrix (3,3,[ cos(phi), sin(phi),  0,
                      -sin(phi), cos(phi),  0,
                        0     ,    0    , 1]);
   Asph := matrix (3,3,[ sin(theta)*cos(phi), sin(theta)*sin(phi),  cos(theta),
                   cos(theta)*cos(phi), cos(theta)*sin(phi), -sin(theta),
                   -sin(phi)          , cos(phi)           , 0]);
   Bcyl := simplify(inverse(Acyl));
   Bsph := simplify(inverse(Asph));
 * \endcode
 *
 * \param [in]  fvtk    pointer to file
 * \param [in]  V       a 4D array [nv][k][j][i] containing the vector
 *                      components (nv) defined at cell centers (k,j,i).
 *                      The index nv = 0 marks the vector first component.
 * \param [in] unit     the corresponding cgs unit (if specified, 1 otherwise)
 * \param [in] var_name the variable name appearing in the VTK file
 * \param [in]    grid  pointer to an array of Grid structures
 *********************************************************************** */
{
    int i, j, k, ndust;
    int vel_field, mag_field;
    int dust_field, dust_num;
    char header[512];
    char dustname[512];
    static Float_Vect ***vect3D;
    double v[3], x1, x2, x3;

    if (vect3D == NULL)
    {
        vect3D = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, Float_Vect);
    }

    /* --------------------------------------------------------
               Write VTK vector fields
   -------------------------------------------------------- */

    v[0] = v[1] = v[2] = 0.0;
    x1 = x2 = x3 = 0.0;
    vel_field = (strcmp(var_name, "vx1") == 0);
    mag_field = (strcmp(var_name, "Bx1") == 0);
    dust_field = 0;
#if DUST == FLUID
    for (ndust = 0; ndust < NDUST; ndust++)
    {
        sprintf(dustname, "vx1_d%d", ndust);

        if (strcmp(var_name, dustname) == 0)
        {
            dust_field = 1;
            dust_num = ndust;
        }
    }
#endif
    if (vel_field || mag_field || dust_field)
    {
        DOM_LOOP(k, j, i)
        {
            D_EXPAND(v[0] = V[0][k][j][i]; x1 = grid->x[IDIR][i];,
                                                                 v[1] = V[1][k][j][i];
                     x2 = grid->x[JDIR][j];,
                                           v[2] = V[2][k][j][i];
                     x3 = grid->x[KDIR][k];)

            VectorCartesianComponents(v, x1, x2, x3);
            vect3D[k][j][i].v1 = (float)v[0] * unit;
            vect3D[k][j][i].v2 = (float)v[1] * unit;
            vect3D[k][j][i].v3 = (float)v[2] * unit;

            if (IsLittleEndian())
            {
                SWAP_VAR(vect3D[k][j][i].v1);
                SWAP_VAR(vect3D[k][j][i].v2);
                SWAP_VAR(vect3D[k][j][i].v3);
            }

        } /* endfor DOM_LOOP(k,j,i) */

        if (vel_field)
            sprintf(header, "\nVECTORS Velocity_Field float\n");
#if DUST == FLUID
        else if (dust_field)
            sprintf(header, "\nVECTORS Velocity_Field_D%d float\n", dust_num);
#endif
        else
            sprintf(header, "\nVECTORS Magnetic_Field float\n");

        VTK_HEADER_WRITE_STRING(header);
        FileWriteData(vect3D[0][0], sizeof(Float_Vect), SZ_Float_Vect, fvtk, -1);
    }
}

/* ********************************************************************* */
void WriteVTK_Scalar(FILE *fvtk, double ***V, double unit,
                     char *var_name, Grid *grid)
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
    float ***Vflt;

    sprintf(header, "\nSCALARS %s float\n", var_name);
    sprintf(header + strlen(header), "LOOKUP_TABLE default\n");

#ifdef PARALLEL
    MPI_Barrier(MPI_COMM_WORLD);
    AL_Write_header(header, strlen(header), MPI_CHAR, SZ_float);
    MPI_Barrier(MPI_COMM_WORLD);
#else
    fprintf(fvtk, "%s", header);
#endif

    Vflt = Convert_dbl2flt(V, unit, IsLittleEndian());
    FileWriteData(Vflt[0][0], sizeof(float), SZ_float, fvtk, -1);
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