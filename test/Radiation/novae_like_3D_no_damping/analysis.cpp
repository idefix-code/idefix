// ***********************************************************************************
// Idefix MHD astrophysical code
//
// Source file test/Radiation/novae_like_3D_no_damping/analysis.cpp
//
// Last modified : 05/2025
//
// Copyright(C) by :
// - Nicolas Scepi <nicolas.scepi@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2026)
// and other code contributors
//
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include "analysis.hpp"
#include "idefix.hpp"
#include "fluid.hpp"
#include <iostream>
#include <fstream>

Analysis::Analysis(Grid &grid, DataBlock &data, std::string filename, real epsilon) : grid(grid), d(data), filename(filename) {
  maskDisk = MakeDiskMask(epsilon);
  maskTop = MakeSurfaceMask(epsilon,1);
  maskBottom = MakeSurfaceMask(epsilon,-1);
  maskMid = MakeMidMask();

  // This hack ensures that d.Vc is an array distinct from data.hydro->Vc, even on CPUs
  this->d.Vc = Kokkos::create_mirror(data.hydro->Vc);
  this->Ex3Ideal = IdefixHostArray3D<real>("Ex3Ideal", data.np_tot[KDIR],data.np_tot[JDIR],data.np_tot[IDIR]);
}

std::vector<real> Analysis::MakeDiskMask(real epsilon) {
    std::vector<real> mask(d.np_tot[JDIR]);

    for(int j = 0; j < d.np_tot[JDIR] ; j++) {
      real th = d.x[JDIR](j);
      if( FABS(cos(th)) <= epsilon) {
        mask[j] = sin(th)*d.dx[JDIR](j);
      } else {
        mask[j] = 0.0;
      }
    }
    return(mask);

}

std::vector<real> Analysis::MakeSurfaceMask(real epsilon, int sign) {
    GridHost gh(grid);
    gh.SyncFromDevice();
    std::vector<real> mask(d.np_tot[JDIR], 0.0);
    int j;
    for(int index = gh.nghost[JDIR]; index < gh.np_tot[JDIR] - gh.nghost[JDIR] ; index++) {
      if(sign==1) {
        // We do it from the north
        j = index;
      } else {
        // We do it from the south
        j = gh.np_tot[JDIR] - 1 - index;
      }
      real th = gh.x[JDIR](j);
      if( FABS(cos(th)) <= epsilon) {
        idfx::cout << "Analysis: Surface at th=" << th << std::endl;
        break;
      }
    }
    int jloc =  j - d.gbeg[JDIR] + d.beg[JDIR];
    if( jloc >= d.beg[JDIR] && jloc < d.end[JDIR] ) {
      mask[jloc] = 1.0;
    }

    return(mask);
}

std::vector<real> Analysis::MakeMidMask() {
  GridHost gh(grid);
  gh.SyncFromDevice();
  std::vector<real> mask(d.np_tot[JDIR], 0.0);
  int jmid = 0;
  for(int j = gh.nghost[JDIR]; j < gh.np_tot[JDIR] - gh.nghost[JDIR] ; j++) {
      real th = gh.x[JDIR](j);
      if(cos(th) < 0.0) {
        real thm1 = gh.x[JDIR](j-1);
        if(FABS(cos(th)) > FABS(cos(thm1))) {
          jmid = j-1;
        } else {
          jmid = j;
        }
        idfx::cout << "Analysis: midplane at th=" << th << std::endl;
        break;
      }
  }
  int jloc =  jmid - d.gbeg[JDIR]+ d.beg[JDIR];
  if( jloc >= d.beg[JDIR] && jloc < d.end[JDIR] ) {
    mask[jloc] = 1.0;
  }
  return(mask);
}


/* **************************************************************** */
void Analysis::Average(IdefixHostArray4D<real> Vin, std::vector<int> fields, std::vector<real> mask, std::vector<real> &outField)
    /*
     * compute the vertical weighted average: int dphi dz rho *infield/int dphi dz rho
     *
     **************************************************************** */
{
  std::vector<real> locField(grid.np_int[IDIR], 0.0);


  for(int k = d.beg[KDIR]; k < d.end[KDIR] ; k++) {
    for(int j = d.beg[JDIR]; j < d.end[JDIR] ; j++) {
      for(int i = d.beg[IDIR]; i < d.end[IDIR] ; i++) {
        int iglob = i - 2*d.nghost[IDIR] + d.gbeg[IDIR];
        real q=1.0;
        for(int n=0 ; n < fields.size() ; n++) {
          q = q*Vin(fields[n],k,j,i);
        }
        locField.at(iglob) += q*mask[j];
      }
    }
  }

    // Reduce
#ifdef WITH_MPI
  MPI_Allreduce(locField.data(), outField.data(), grid.np_int[IDIR], realMPI ,MPI_SUM, MPI_COMM_WORLD);
#else
  outField = locField;
#endif
  #if DIMENSIONS == 3
    // Divided by the number of points in phi
    for(int i = 0 ; i < grid.np_int[IDIR] ; i++) {
      outField[i] /= grid.np_int[KDIR];
    }
  #endif


  return ;
}

/* **************************************************************** */
void Analysis::Average(IdefixHostArray3D<real> Vin, std::vector<real> mask, std::vector<real> &outField)
    /*
     * compute the vertical weighted average: int dphi dz rho *infield/int dphi dz rho
     *
     **************************************************************** */
{
  std::vector<real> locField(grid.np_int[IDIR], 0.0);


  for(int k = d.beg[KDIR]; k < d.end[KDIR] ; k++) {
    for(int j = d.beg[JDIR]; j < d.end[JDIR] ; j++) {
      for(int i = d.beg[IDIR]; i < d.end[IDIR] ; i++) {
        int iglob = i - 2*d.nghost[IDIR] + d.gbeg[IDIR];
        locField.at(iglob) += Vin(k,j,i)*mask[j];
      }
    }
  }

    // Reduce
#ifdef WITH_MPI
  MPI_Allreduce(locField.data(), outField.data(), grid.np_int[IDIR], realMPI ,MPI_SUM, MPI_COMM_WORLD);
#else
  outField = locField;
#endif
  #if DIMENSIONS == 3
    // Divided by the number of points in phi
    for(int i = 0 ; i < grid.np_int[IDIR] ; i++) {
      outField[i] /= grid.np_int[KDIR];
    }
  #endif


  return ;
}



void Analysis::SubstractKeplerian(DataBlockHost &d) {

  for(int k = 0 ; k < d.np_tot[KDIR] ; k++) {
    for(int j = 0 ; j < d.np_tot[JDIR] ; j++) {
      for(int i = 0 ; i < d.np_tot[IDIR] ; i++) {
        real Vk = 1/sqrt(d.x[IDIR](i))/sin(d.x[JDIR](j));
        d.Vc(VX3,k,j,i) -= Vk;
      }
    }
  }
}

void Analysis::ComputeEMF(DataBlock &data) {
  ConstrainedTransport<DefaultPhysics> *emf = data.hydro->emf.get();
  // Talk to the emf object to recompute required EMFs
  IdefixArray3D<real> ex3 = emf->ez;

  // Compute Ideal EMF
    data.SetBoundaries();
  emf->CalcCornerEMF(data.t);
  Kokkos::deep_copy(this->Ex3Ideal,ex3);

  idefix_for("resetEMF", 0, data.np_tot[KDIR], 0, data.np_tot[JDIR], 0, data.np_tot[IDIR],
      KOKKOS_LAMBDA(int k, int j,int i) {
        ex3(k,j,i) = 0.0;
      });

  // We're done with the EMFs
}
/* **************************************************************** */
void Analysis::WriteProfile(real t, std::vector<real> &data, std::ofstream &file) {
/*
 * Write a global profile to a file
 *
 *
 **************************************************************** */
  if(idfx::prank==0) {
    file << std::scientific << t;
    for(int i = 0 ; i < grid.np_int[IDIR] ; i++ ) {
      file << std::scientific << "\t" << data[i];
    }
    file << std::endl;
  }
  return ;
}


void Analysis::ResetAnalysis() {
  GridHost gh(grid);
  gh.SyncFromDevice();
  std::ofstream file;
  if(idfx::prank==0) {
    file.open(filename, std::ios::trunc);
    file << "rho\t prs\t rhovr \t rhovrvphi \t rhovth_top\t rhovth_bot\t rhovthvphi_top\t rhovthvphi_bot\t BrBphi \t BthBphi_top \t BthBphi_bot\t Bth_mid \t P_mid\t rhovth\t rhovth_mid";
    #ifdef EVOLVE_VECTOR_POTENTIAL
    file << "\t A_mid \t A_top \t A_bot";
    #endif
    file << "\t Ex3Id_mid \t Ex3Id_top \t Ex3Id_bot \t Ex3Ni_mid \t Ex3Ni_top \t Ex3Ni_bot";
    file << std::endl;
    file.precision(10);
    // Put a 0 in time column
    file << std::scientific << "0.0";
    // Write radial coordinate
    for(int i = 0 ; i < grid.np_int[IDIR] ; i++ ) {
      file << std::scientific << "\t" << gh.x[IDIR](i+gh.nghost[IDIR]);
    }
    file << std::endl;
    file.close();


  }
}

void Analysis::PerformAnalysis(DataBlock &data) {
  idfx::pushRegion("Analysis::PerformAnalysis");
  if(data.t == 0) {
      this->ResetAnalysis();
  }
  d.SyncFromDevice();
  SubstractKeplerian(d);

  std::vector<real> prof(grid.np_int[IDIR]);
  std::ofstream file;

  if(idfx::prank==0) {
    file.open(filename, std::ios::app);
    file.precision(10);
  }


  // Z averages (profils along R)
  Average(d.Vc, std::vector<int> {RHO}, maskDisk, prof);
  WriteProfile(data.t, prof, file);

  Average(d.Vc, std::vector<int> {PRS}, maskDisk, prof);
  WriteProfile(data.t, prof, file);

  Average(d.Vc, std::vector<int> {RHO, VX1}, maskDisk, prof);
  WriteProfile(data.t, prof, file);

  Average(d.Vc, std::vector<int> {RHO, VX1, VX3}, maskDisk, prof);
  WriteProfile(data.t, prof, file);

  Average(d.Vc, std::vector<int> {RHO, VX2}, maskTop, prof);
  WriteProfile(data.t, prof, file);

  Average(d.Vc, std::vector<int> {RHO, VX2}, maskBottom, prof);
  WriteProfile(data.t, prof, file);

  Average(d.Vc, std::vector<int> {RHO, VX2, VX3}, maskTop, prof);
  WriteProfile(data.t, prof, file);

  Average(d.Vc, std::vector<int> {RHO, VX2, VX3}, maskBottom, prof);
  WriteProfile(data.t, prof, file);

  Average(d.Vc, std::vector<int> {BX1, BX3}, maskDisk, prof);
  WriteProfile(data.t, prof, file);

  Average(d.Vc, std::vector<int> {BX2, BX3}, maskTop, prof);
  WriteProfile(data.t, prof, file);

  Average(d.Vc, std::vector<int> {BX2, BX3}, maskBottom, prof);
  WriteProfile(data.t, prof, file);

  Average(d.Vc, std::vector<int> {BX2}, maskMid, prof);
  WriteProfile(data.t, prof, file);

  Average(d.Vc, std::vector<int> {PRS}, maskMid, prof);
  WriteProfile(data.t, prof, file);

  Average(d.Vc, std::vector<int> {RHO,VX2}, maskDisk, prof);
  WriteProfile(data.t, prof, file);

  Average(d.Vc, std::vector<int> {RHO,VX2}, maskMid, prof);
  WriteProfile(data.t, prof, file);

  #ifdef EVOLVE_VECTOR_POTENTIAL
    Average(d.Ve, std::vector<int> {AX3e}, maskMid, prof);
    WriteProfile(data.t, prof, file);

    Average(d.Ve, std::vector<int> {AX3e}, maskTop, prof);
    WriteProfile(data.t, prof, file);

    Average(d.Ve, std::vector<int> {AX3e}, maskBottom, prof);
    WriteProfile(data.t, prof, file);
  #endif

  this->ComputeEMF(data);
  Average(this->Ex3Ideal, maskMid, prof);
  WriteProfile(data.t, prof, file);
  Average(this->Ex3Ideal, maskTop, prof);
  WriteProfile(data.t, prof, file);
  Average(this->Ex3Ideal, maskBottom, prof);
  WriteProfile(data.t, prof, file);

  if(idfx::prank==0) {
    file.close();
  }
  idfx::popRegion();
}
