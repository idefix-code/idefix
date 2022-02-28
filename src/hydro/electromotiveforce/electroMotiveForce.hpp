// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2022 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef HYDRO_ELECTROMOTIVEFORCE_ELECTROMOTIVEFORCE_HPP_
#define  HYDRO_ELECTROMOTIVEFORCE_ELECTROMOTIVEFORCE_HPP_

#include "idefix.hpp"
#include "input.hpp"

// Forward declarations
class Hydro;
class DataBlock;

class ElectroMotiveForce {
 public:
  enum AveragingType {none, arithmetic, uct0, uct_contact, uct_hll, uct_hlld};

  // Type of averaging
  AveragingType averaging{none};

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

// Required by uct_contact averaging
  IdefixArray3D<int>      svx;
  IdefixArray3D<int>      svy;
  IdefixArray3D<int>      svz;

// required by uct_hll averaging
  IdefixArray3D<real> axL;
  IdefixArray3D<real> axR;
  IdefixArray3D<real> ayL;
  IdefixArray3D<real> ayR;
  #if DIMENSIONS == 3
  IdefixArray3D<real> azL;
  IdefixArray3D<real> azR;

  IdefixArray3D<real> dzL;
  IdefixArray3D<real> dzR;
  #endif
  IdefixArray3D<real> dxL;
  IdefixArray3D<real> dxR;
  IdefixArray3D<real> dyL;
  IdefixArray3D<real> dyR;


  IdefixArray3D<real>     Ex1;
  IdefixArray3D<real>     Ex2;
  IdefixArray3D<real>     Ex3;

  // Helper arrays for shearing box
  IdefixArray2D<real>     sbEyL;
  IdefixArray2D<real>     sbEyR;
  IdefixArray2D<real>     sbEyRL;

  // Range of existence

  // Init from Hydro class
  void Init(Input &, Hydro *);

  void EvolveMagField(real, real, IdefixArray4D<real>&);
  void CalcCornerEMF(real );
  void ShowConfig();

  // Different flavors of EMF average schemes
  void CalcRiemannAverage();
  void CalcArithmeticAverage();
  void CalcCellCenteredEMF();
  void CalcUCT0Average();
  void CalcContactAverage();

  // Enforce boundary conditions on the EMFs.
  void EnforceEMFBoundary();
  void CalcNonidealEMF(real );

  // Specific routines to symmetrize EMFs with shearing box BCs
  void SymmetrizeEMFShearingBox();
  void ExtrapolateEMFShearingBox(BoundarySide,
                            IdefixArray2D<real>,
                            IdefixArray2D<real>);

  // Routines for evolving the magnetic potential (only available when EVOLVE_VECTOR_POTENTIAL)
  void EvolveVectorPotential(real, IdefixArray4D<real> &);
  void ComputeMagFieldFromA(IdefixArray4D<real> &Vein, IdefixArray4D<real> &Vsout);

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
