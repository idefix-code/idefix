// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2021 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef HYDRO_ELECTROMOTIVEFORCE_ELECTROMOTIVEFORCE_HPP_
#define  HYDRO_ELECTROMOTIVEFORCE_ELECTROMOTIVEFORCE_HPP_

#include "idefix.hpp"

// Forward declarations
class Hydro;
class DataBlock;

class ElectroMotiveForce {
 public:
  // Face centered emf components
  IdefixArray3D<real>     exj;
  IdefixArray3D<real>     exk;
  IdefixArray3D<real>     eyi;
  IdefixArray3D<real>     eyk;
  IdefixArray3D<real>     ezi;
  IdefixArray3D<real>     ezj;

  // Edge centered emf components
  IdefixArray3D<real>     ex;
  IdefixArray3D<real>     ey;
  IdefixArray3D<real>     ez;

#if EMF_AVERAGE == UCT_CONTACT
  IdefixArray3D<int>      svx;
  IdefixArray3D<int>      svy;
  IdefixArray3D<int>      svz;

#elif EMF_AVERAGE == UCT_HLL
  // Signal velocities
  IdefixArray3D<real> SxL;
  IdefixArray3D<real> SxR;
  IdefixArray3D<real> SyL;
  IdefixArray3D<real> SyR;
  IdefixArray3D<real> SzL;
  IdefixArray3D<real> SzR;

  // Staggered magnetic field and velocity slopes
  IdefixArray3D<real> dbx_dy, dby_dx;
  #if DIMENSIONS == 3
  IdefixArray3D<real> dbz_dx, dbz_dy;
  IdefixArray3D<real> dbx_dz, dby_dz;
  #endif

  IdefixArray3D<real> dvx_dx, dvx_dy;
  IdefixArray3D<real> dvy_dx, dvy_dy;
  #if DIMENSIONS == 3
  IdefixArray3D<real> dvx_dz, dvy_dz;
  IdefixArray3D<real> dvz_dx, dvz_dy, dvz_dz;
  #endif
#endif

  IdefixArray3D<real>     Ex1;
  IdefixArray3D<real>     Ex2;
  IdefixArray3D<real>     Ex3;

  // Helper arrays for shearing box
  IdefixArray2D<real>     sbEyL;
  IdefixArray2D<real>     sbEyR;
  IdefixArray2D<real>     sbEyRL;

  // Range of existence
  int  ibeg, jbeg, kbeg;
  int  iend, jend, kend;

  // Init from Hydro class
  void Init(Hydro *);

  void EvolveMagField(real, real, IdefixArray4D<real>&);
  void CalcCornerEMF(real );
  void calcRiemannEmf();
  // Enforce boundary conditions on the EMFs.
  void EnforceEMFBoundary();
  void CalcNonidealEMF(real );

  // Specific routines to symmetrize EMFs with shearing box BCs
  void SymmetrizeEMFShearingBox();
  void ExtrapolateEMFShearingBox(BoundarySide,
                            IdefixArray2D<real>,
                            IdefixArray2D<real>);



#ifdef WITH_MPI
  // Exchange surface EMFs to remove interprocess round off errors
  void ExchangeAll();
  void ExchangeX1();
  void ExchangeX2();
  void ExchangeX3();
#endif

  // Default constructor
  ElectroMotiveForce();

  // Destructor
  ~ElectroMotiveForce();

 private:
  DataBlock *data;
  Hydro *hydro;

#ifdef WITH_MPI
  enum {faceRight, faceLeft};

  // Buffers for MPI calls
  IdefixArray1D<real> BufferSendX1[2];
  IdefixArray1D<real> BufferSendX2[2];
  IdefixArray1D<real> BufferSendX3[2];
  IdefixArray1D<real> BufferRecvX1[2];
  IdefixArray1D<real> BufferRecvX2[2];
  IdefixArray1D<real> BufferRecvX3[2];

  IdefixArray1D<int>  mapVars;

  int bufferSizeX1;
  int bufferSizeX2;
  int bufferSizeX3;

  // Requests for MPI persistent communications
  MPI_Request sendRequestX1[2];
  MPI_Request sendRequestX2[2];
  MPI_Request sendRequestX3[2];
  MPI_Request recvRequestX1[2];
  MPI_Request recvRequestX2[2];
  MPI_Request recvRequestX3[2];

  Kokkos::Timer timer;    // Internal MPI timer
#endif
};

#endif // HYDRO_ELECTROMOTIVEFORCE_ELECTROMOTIVEFORCE_HPP_
