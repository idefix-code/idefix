// ***********************************************************************************
// Idefix MHD astrophysical code
//
// Source file test/Radiation/novae_like_3D_no_damping/analysis.hpp
//
// Last modified : 05/2025
//
// Copyright(C) by :
// - Nicolas Scepi <nicolas.scepi@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2026)
// and other code contributors
//
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef ANALYSIS_HPP_
#define ANALYSIS_HPP_

#include "idefix.hpp"
#include "input.hpp"
#include "output.hpp"
#include "grid.hpp"
#include "dataBlock.hpp"
#include "dataBlockHost.hpp"
#include <vector>
#include <string>



class Analysis {
 public:
  // Constructor from Setup arguments
  Analysis(Grid &grid, DataBlock &data, std::string filename, real);
  void ResetAnalysis();
  void PerformAnalysis(DataBlock &);
  void ComputeEMF(DataBlock &);
 private:
  void Average(IdefixHostArray4D<real> Vin, std::vector<int> fields, std::vector<real> mask, std::vector<real> &outField);
  void Average(IdefixHostArray3D<real> Vin, std::vector<real> mask, std::vector<real> &outField);
  void WriteProfile(real t, std::vector<real> &data, std::ofstream &file);
  void SubstractKeplerian(DataBlockHost&);
  std::vector<real> MakeDiskMask(real epsilon);
  std::vector<real> MakeSurfaceMask(real epsilon, int sign);
  std::vector<real> MakeMidMask();


  DataBlockHost d;
  IdefixHostArray3D<real> Ex3Ideal, Ex3Nonideal;

  Grid &grid;
  std::string filename;
  std::vector<real> maskTop;
  std::vector<real> maskBottom;
  std::vector<real> maskDisk;
  std::vector<real> maskMid;
};

#endif // ANALYSIS_HPP__
