// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef FLUID_BOUNDARY_AXIS_HPP_
#define FLUID_BOUNDARY_AXIS_HPP_

#include <vector>
#include "idefix.hpp"
#include "grid.hpp"

// Forward class hydro declaration
#include "physics.hpp"
template <typename Phys> class Fluid;
template <typename Phys> class ConstrainedTransport;
template <typename Phys> class Boundary;
class DataBlock;

// Whether we use athena++ procedure to regularise BX2s
#define AXIS_BX2S_USE_ATHENA_REGULARISATION

class Axis {
 public:
  template <typename Phys>
  explicit Axis(Boundary<Phys> *);  // Initialisation
  void RegularizeEMFs();                 // Regularize the EMF sitting on the axis
  void RegularizeCurrent();             // Regularize the currents along the axis
  void EnforceAxisBoundary(int side);   // Enforce the boundary conditions (along X2)
  void ReconstructBx2s();               // Reconstruct BX2s in the ghost zone using divB=0
  void ShowConfig();


  void SymmetrizeEx1Side(int);         // Symmetrize on a specific side (internal method)
  void RegularizeEx3side(int);         // Regularize Ex3 along the axis (internal method)
  void RegularizeCurrentSide(int);      // Regularize J along the axis (internal method)
  void FixBx2sAxis(int side);           // Fix BX2s on the axis using the field around it (internal)
  void FixBx2sAxisGhostAverage(int side); //Fix BX2s on the axis using the average of neighbouring
                                          // cell in theta direction (like Athena)
  void ExchangeMPI(int side);           // Function has to be public for GPU, but its technically
                                        // a private function


 private:
  bool isTwoPi = false;
  bool axisRight = false;
  bool axisLeft = false;
  bool needMPIExchange = false;
  bool haveMHD = false;
  int nVar;

  enum {faceTop, faceBot};
#ifdef WITH_MPI
  MPI_Request sendRequest;
  MPI_Request recvRequest;

  IdefixArray1D<real> bufferSend;
  IdefixArray1D<real> bufferRecv;

  int bufferSize;

  IdefixArray1D<int>  mapVars;
  int mapNVars{0};

#endif
  void InitMPI();

  IdefixArray1D<real> Ex1Avg;
  IdefixArray2D<real> BAvg;
  bool haveCurrent;
  IdefixArray2D<real> JAvg;
  IdefixArray1D<int> symmetryVc;
  IdefixArray1D<int> symmetryVs;

  IdefixArray3D<real> ex;
  IdefixArray3D<real> ey;
  IdefixArray3D<real> ez;
  IdefixArray4D<real> J;

  IdefixArray4D<real> Vc;
  IdefixArray4D<real> Vs;

  DataBlock *data;
};

#include "constrainedTransport.hpp"

template<typename Phys>
Axis::Axis(Boundary<Phys> *boundary) {
  Vc = boundary->Vc;
  Vs = boundary->Vs;
  J = boundary->fluid->J;
  if constexpr(Phys::mhd) {
    ex = boundary->fluid->emf->ex;
    ey = boundary->fluid->emf->ey;
    ez = boundary->fluid->emf->ez;
  }

  data = boundary->data;
  haveMHD = Phys::mhd;
  nVar = boundary->nVar;
  haveCurrent = boundary->fluid->haveCurrent;


  #if GEOMETRY != SPHERICAL
    IDEFIX_ERROR("Axis boundary conditions are only designed to handle spherical geometry");
  #endif


  if(fabs((data->mygrid->xend[KDIR] - data->mygrid->xbeg[KDIR] -2.0*M_PI)) < 1e-10) {
    this->isTwoPi = true;
    #ifdef WITH_MPI
      // Check that there is a domain decomposition in phi
      if(data->mygrid->nproc[KDIR]>1) {
        if(data->mygrid->nproc[KDIR]%2==1) {
          IDEFIX_ERROR("The numbre of processes in the phi direction should"
                        " be even for axis decomposition");
        }
        needMPIExchange = true;
      }
    #endif
  } else {
    this->isTwoPi = false;
  }

  // Check where the axis is lying.
  if(data->lbound[JDIR] == axis) axisLeft = true;
  if(data->rbound[JDIR] == axis) axisRight = true;

  // Init the symmetry array (used to flip the signs of arrays accross the axis)
  symmetryVc = IdefixArray1D<int>("Axis:SymmetryVc",nVar);
  IdefixArray1D<int>::HostMirror symmetryVcHost = Kokkos::create_mirror_view(symmetryVc);
  // Init the array
  for (int nv = 0; nv < nVar; nv++) {
    symmetryVcHost(nv) = 1;
    if (nv == VX2)
      symmetryVcHost(nv) = -1;
    if (nv == VX3)
      symmetryVcHost(nv) = -1;
    if (nv == BX2)
      symmetryVcHost(nv) = -1;
    if (nv == BX3)
      symmetryVcHost(nv) = -1;
  }
  Kokkos::deep_copy(symmetryVc, symmetryVcHost);

  if constexpr(Phys::mhd) {
    idfx::cout << "Phys MHD" << std::endl;
    symmetryVs = IdefixArray1D<int>("Axis:SymmetryVs",DIMENSIONS);
    IdefixArray1D<int>::HostMirror symmetryVsHost = Kokkos::create_mirror_view(symmetryVs);
    // Init the array
    for(int nv = 0; nv < DIMENSIONS; nv++) {
      symmetryVsHost(nv) = 1;
      if (nv == BX2s)
        symmetryVsHost(nv) = -1;
      if (nv == BX3s)
        symmetryVsHost(nv) = -1;
    }
    Kokkos::deep_copy(symmetryVs, symmetryVsHost);

    this->Ex1Avg = IdefixArray1D<real>("Axis:Ex1Avg",data->np_tot[IDIR]);
    this->BAvg = IdefixArray2D<real>("Axis:BxAvg",data->np_tot[IDIR],2);
    if(haveCurrent) {
      this->JAvg = IdefixArray2D<real>("Axis:JAvg",data->np_tot[IDIR],3);
    }
  }

  #ifdef WITH_MPI
    if(needMPIExchange) {
      // Make MPI exchange datatypes
      InitMPI();
    }
  #endif
}



#endif // FLUID_BOUNDARY_AXIS_HPP_
