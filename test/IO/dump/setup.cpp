// ***********************************************************************************
// Idefix MHD astrophysical code
//
// Source file test/IO/dump/setup.cpp
//
// Last modified : 03/2026
//
// Copyright(C) by :
// - Geoffroy Lesur <geoffroy.lesur@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2026)
// and other code contributors
//
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include "idefix.hpp"
#include "setup.hpp"

/*********************************************/
/**
Customized random number generator
Allow one to have consistent random numbers
generators on different architectures.
**/
/*********************************************/

Output* myOutput;
int outnum;
// =====================================================================================
// Scalars covering every data type supported by the dump I/O.
// They are registered with the Dump, written once, then randomized so that only the
// dump file can restore their original value.
// =====================================================================================
struct DumpScalars {
  int       int_      = 42;
  int64_t   int64_    = -1234567890123LL;
  float     float_    = 1.25f;
  double    double_   = 3.14159265358979;
  bool      bool_     = true;
  std::complex<float>  cfloat_  = {0.5f, -3.0f};
  std::complex<double> cdouble_ = {1.5, -2.25};
  Kokkos::complex<float>  kcfloat_  = {-7.0f, 8.5f};
  Kokkos::complex<double> kcdouble_ = {4.0, 5.5};
};
DumpScalars scalars;
// Snapshot of the original values, used for validation after the round-trip
DumpScalars scalarsRef;

// =====================================================================================
// Replicated 1D/2D arrays covering every data type supported by the dump I/O.
// Rank 0 writes them once; every rank receives an identical copy on restart.
// =====================================================================================
constexpr int N1D = 7;    // length of 1D test arrays
constexpr int N2Da = 4;   // first extent of 2D test arrays
constexpr int N2Db = 3;   // second extent of 2D test arrays

// Kokkos views cannot be constructed before Kokkos::initialize(), so they are
// allocated in a holder created from Setup::Setup (after initialisation).
struct TestArrays {
  IdefixArray1D<int>       int1D;
  IdefixArray1D<int64_t>   int641D;
  IdefixArray1D<float>     float1D;
  IdefixArray1D<double>    double1D;
  IdefixArray1D<bool>      bool1D;
  IdefixArray2D<double>    double2D;
  IdefixArray2D<bool>      bool2D;
  IdefixArray2D<int>       int2D;
  IdefixArray2D<Kokkos::complex<float>>  cfloat2D;
  IdefixArray2D<Kokkos::complex<double>> cdouble2D;

  TestArrays()
    : int1D("arrInt1D", N1D),
      int641D("arrInt641D", N1D),
      float1D("arrFloat1D", N1D),
      double1D("arrDouble1D", N1D),
      bool1D("arrBool1D", N1D),
      double2D("arrDouble2D", N2Da, N2Db),
      bool2D("arrBool2D", N2Da, N2Db),
      int2D("arrInt2D", N2Da, N2Db),
      cfloat2D("arrCFloat2D", N2Da, N2Db),
      cdouble2D("arrCDouble2D", N2Da, N2Db) {}
};

std::unique_ptr<TestArrays> testArrays;

// Fill every array with values a wrong packing would scramble
void InitArrays() {
  auto hInt1D    = Kokkos::create_mirror_view(testArrays->int1D);
  auto hInt641D  = Kokkos::create_mirror_view(testArrays->int641D);
  auto hFloat1D  = Kokkos::create_mirror_view(testArrays->float1D);
  auto hDouble1D = Kokkos::create_mirror_view(testArrays->double1D);
  auto hBool1D   = Kokkos::create_mirror_view(testArrays->bool1D);
  auto hDouble2D = Kokkos::create_mirror_view(testArrays->double2D);
  auto hBool2D   = Kokkos::create_mirror_view(testArrays->bool2D);
  auto hInt2D    = Kokkos::create_mirror_view(testArrays->int2D);
  auto hCf2D     = Kokkos::create_mirror_view(testArrays->cfloat2D);
  auto hCd2D     = Kokkos::create_mirror_view(testArrays->cdouble2D);
  for(int i = 0 ; i < N1D ; i++) {
    hInt1D(i)    = 100 + i;
    hInt641D(i)  = -9000000000000LL + i;
    hFloat1D(i)  = 0.5f * i + 1.0f;
    hDouble1D(i) = -0.25 * i + 2.0;
    hBool1D(i)   = (i % 2 == 0);
  }
  for(int i = 0 ; i < N2Da ; i++) {
    for(int j = 0 ; j < N2Db ; j++) {
      hDouble2D(i, j) = 10.0 * i + j + 0.5;
      hBool2D(i, j)   = ((i + j) % 2 == 0);
      hInt2D(i, j)    = 100 * i + j;
      hCf2D(i, j)     = Kokkos::complex<float>(i + 0.5f, -j - 1.5f);
      hCd2D(i, j)     = Kokkos::complex<double>(2.0 * i, -3.0 * j);
    }
  }
  Kokkos::deep_copy(testArrays->int1D, hInt1D);
  Kokkos::deep_copy(testArrays->int641D, hInt641D);
  Kokkos::deep_copy(testArrays->float1D, hFloat1D);
  Kokkos::deep_copy(testArrays->double1D, hDouble1D);
  Kokkos::deep_copy(testArrays->bool1D, hBool1D);
  Kokkos::deep_copy(testArrays->double2D, hDouble2D);
  Kokkos::deep_copy(testArrays->bool2D, hBool2D);
  Kokkos::deep_copy(testArrays->int2D, hInt2D);
  Kokkos::deep_copy(testArrays->cfloat2D, hCf2D);
  Kokkos::deep_copy(testArrays->cdouble2D, hCd2D);
}

// Compare the current arrays against a reference set of host copies
int CheckArraysAgainst(const IdefixHostArray1D<int> &refInt1D,
                       const IdefixHostArray1D<int64_t> &refInt641D,
                       const IdefixHostArray1D<float> &refFloat1D,
                       const IdefixHostArray1D<double> &refDouble1D,
                       const IdefixHostArray1D<bool> &refBool1D,
                       const IdefixHostArray2D<double> &refDouble2D,
                       const IdefixHostArray2D<bool> &refBool2D,
                       const IdefixHostArray2D<int> &refInt2D,
                       const IdefixHostArray2D<Kokkos::complex<float>> &refCf2D,
                       const IdefixHostArray2D<Kokkos::complex<double>> &refCd2D);

// Reference host copies, filled by SnapshotArrays() before the dump is written
struct ArrayRefs {
  IdefixHostArray1D<int>     int1D;
  IdefixHostArray1D<int64_t> int641D;
  IdefixHostArray1D<float>   float1D;
  IdefixHostArray1D<double>  double1D;
  IdefixHostArray1D<bool>    bool1D;
  IdefixHostArray2D<double>  double2D;
  IdefixHostArray2D<bool>    bool2D;
  IdefixHostArray2D<int>     int2D;
  IdefixHostArray2D<Kokkos::complex<float>>  cfloat2D;
  IdefixHostArray2D<Kokkos::complex<double>> cdouble2D;
};
ArrayRefs arrRefs;

void SnapshotArrays() {
  arrRefs.int1D    = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), testArrays->int1D);
  arrRefs.int641D  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), testArrays->int641D);
  arrRefs.float1D  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), testArrays->float1D);
  arrRefs.double1D = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), testArrays->double1D);
  arrRefs.bool1D   = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), testArrays->bool1D);
  arrRefs.double2D = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), testArrays->double2D);
  arrRefs.bool2D   = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), testArrays->bool2D);
  arrRefs.int2D    = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), testArrays->int2D);
  arrRefs.cfloat2D = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), testArrays->cfloat2D);
  arrRefs.cdouble2D = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), testArrays->cdouble2D);
}

// Overwrite every array with values that cannot match the reference snapshot
void RandomizeArrays() {
  {
    auto h = Kokkos::create_mirror_view(testArrays->int1D);
    for(int i = 0 ; i < N1D ; i++) h(i) = -1 - i;
    Kokkos::deep_copy(testArrays->int1D, h);
  }
  {
    auto h = Kokkos::create_mirror_view(testArrays->int641D);
    for(int i = 0 ; i < N1D ; i++) h(i) = 1 + i;
    Kokkos::deep_copy(testArrays->int641D, h);
  }
  {
    auto h = Kokkos::create_mirror_view(testArrays->float1D);
    for(int i = 0 ; i < N1D ; i++) h(i) = 99.0f;
    Kokkos::deep_copy(testArrays->float1D, h);
  }
  {
    auto h = Kokkos::create_mirror_view(testArrays->double1D);
    for(int i = 0 ; i < N1D ; i++) h(i) = -99.0;
    Kokkos::deep_copy(testArrays->double1D, h);
  }
  {
    auto h = Kokkos::create_mirror_view(testArrays->bool1D);
    for(int i = 0 ; i < N1D ; i++) h(i) = true;
    Kokkos::deep_copy(testArrays->bool1D, h);
  }
  {
    auto h = Kokkos::create_mirror_view(testArrays->double2D);
    for(int i = 0 ; i < N2Da ; i++) for(int j = 0 ; j < N2Db ; j++) h(i, j) = -99.0;
    Kokkos::deep_copy(testArrays->double2D, h);
  }
  {
    auto h = Kokkos::create_mirror_view(testArrays->bool2D);
    for(int i = 0 ; i < N2Da ; i++) for(int j = 0 ; j < N2Db ; j++) h(i, j) = true;
    Kokkos::deep_copy(testArrays->bool2D, h);
  }
  {
    auto h = Kokkos::create_mirror_view(testArrays->int2D);
    for(int i = 0 ; i < N2Da ; i++) for(int j = 0 ; j < N2Db ; j++) h(i, j) = -1;
    Kokkos::deep_copy(testArrays->int2D, h);
  }
  {
    auto h = Kokkos::create_mirror_view(testArrays->cfloat2D);
    for(int i = 0 ; i < N2Da ; i++) for(int j = 0 ; j < N2Db ; j++)
      h(i, j) = Kokkos::complex<float>(9.f, 9.f);
    Kokkos::deep_copy(testArrays->cfloat2D, h);
  }
  {
    auto h = Kokkos::create_mirror_view(testArrays->cdouble2D);
    for(int i = 0 ; i < N2Da ; i++) for(int j = 0 ; j < N2Db ; j++)
      h(i, j) = Kokkos::complex<double>(9., 9.);
    Kokkos::deep_copy(testArrays->cdouble2D, h);
  }
}

// Compare the current (reloaded) arrays against the snapshot
int CheckArrays() {
  return CheckArraysAgainst(arrRefs.int1D, arrRefs.int641D, arrRefs.float1D,
                            arrRefs.double1D, arrRefs.bool1D, arrRefs.double2D,
                            arrRefs.bool2D, arrRefs.int2D, arrRefs.cfloat2D,
                            arrRefs.cdouble2D);
}

// Compare the current (reloaded) arrays against the reference host copies
int CheckArraysAgainst(const IdefixHostArray1D<int> &refInt1D,
                       const IdefixHostArray1D<int64_t> &refInt641D,
                       const IdefixHostArray1D<float> &refFloat1D,
                       const IdefixHostArray1D<double> &refDouble1D,
                       const IdefixHostArray1D<bool> &refBool1D,
                       const IdefixHostArray2D<double> &refDouble2D,
                       const IdefixHostArray2D<bool> &refBool2D,
                       const IdefixHostArray2D<int> &refInt2D,
                       const IdefixHostArray2D<Kokkos::complex<float>> &refCf2D,
                       const IdefixHostArray2D<Kokkos::complex<double>> &refCd2D) {
  int errors = 0;
  auto report = [&](const char *name, bool ok) {
    if(!ok) {
      errors++;
      idfx::cout << "-----------------------------------------" << std::endl
                 << " Error: array " << name << " was not restored correctly" << std::endl;
    }
  };

  auto curInt1D    = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), testArrays->int1D);
  auto curInt641D  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), testArrays->int641D);
  auto curFloat1D  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), testArrays->float1D);
  auto curDouble1D = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), testArrays->double1D);
  auto curBool1D   = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), testArrays->bool1D);
  auto curDouble2D = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), testArrays->double2D);
  auto curBool2D   = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), testArrays->bool2D);
  auto curInt2D    = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), testArrays->int2D);
  auto curCf2D     = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), testArrays->cfloat2D);
  auto curCd2D     = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), testArrays->cdouble2D);

  bool okInt1D = true, okInt641D = true, okFloat1D = true, okDouble1D = true;
  bool okBool1D = true, okDouble2D = true, okBool2D = true, okInt2D = true;
  bool okCf2D = true, okCd2D = true;
  for(int i = 0 ; i < N1D ; i++) {
    okInt1D    = okInt1D    && (curInt1D(i)    == refInt1D(i));
    okInt641D  = okInt641D  && (curInt641D(i)  == refInt641D(i));
    okFloat1D  = okFloat1D  && (curFloat1D(i)  == refFloat1D(i));
    okDouble1D = okDouble1D && (curDouble1D(i) == refDouble1D(i));
    okBool1D   = okBool1D   && (bool(curBool1D(i)) == bool(refBool1D(i)));
  }
  for(int i = 0 ; i < N2Da ; i++) {
    for(int j = 0 ; j < N2Db ; j++) {
      okDouble2D = okDouble2D && (curDouble2D(i, j) == refDouble2D(i, j));
      okBool2D   = okBool2D   && (bool(curBool2D(i, j)) == bool(refBool2D(i, j)));
      okInt2D    = okInt2D    && (curInt2D(i, j) == refInt2D(i, j));
      okCf2D     = okCf2D     && (curCf2D(i, j) == refCf2D(i, j));
      okCd2D     = okCd2D     && (curCd2D(i, j) == refCd2D(i, j));
    }
  }
  report("testInt1D",    okInt1D);
  report("testInt641D",  okInt641D);
  report("testFloat1D",  okFloat1D);
  report("testDouble1D", okDouble1D);
  report("testBool1D",   okBool1D);
  report("testDouble2D", okDouble2D);
  report("testBool2D",   okBool2D);
  report("testInt2D",    okInt2D);
  report("testCFloat2D", okCf2D);
  report("testCDouble2D", okCd2D);
  return errors;
}

// Overwrite every scalar with a value that cannot be confused with the original
void RandomizeScalars() {
  scalars.int_      = scalars.int_ + 1;
  scalars.int64_    = -scalars.int64_ - 1;
  scalars.float_    = scalars.float_ * 2.f + 1.f;
  scalars.double_   = -scalars.double_ * 2.0 - 1.0;
  scalars.bool_     = !scalars.bool_;
  scalars.cfloat_   = std::complex<float>(9.f, 9.f);
  scalars.cdouble_  = std::complex<double>(9., 9.);
  scalars.kcfloat_  = Kokkos::complex<float>(9.f, 9.f);
  scalars.kcdouble_ = Kokkos::complex<double>(9., 9.);
}

// Check that the current scalars match the snapshot taken before the dump
int CheckScalars() {
  int errors = 0;
  auto report = [&](const char *name, bool ok) {
    if (!ok) {
      errors++;
      idfx::cout << "-----------------------------------------" << std::endl
                 << " Error: field " << name << " was not restored correctly" << std::endl;
    }
  };
  report("int",         scalars.int_      == scalarsRef.int_);
  report("int64",       scalars.int64_    == scalarsRef.int64_);
  report("float",       scalars.float_    == scalarsRef.float_);
  report("double",      scalars.double_   == scalarsRef.double_);
  report("bool",        scalars.bool_     == scalarsRef.bool_);
  report("complexFloat",  scalars.cfloat_  == scalarsRef.cfloat_);
  report("complexDouble", scalars.cdouble_ == scalarsRef.cdouble_);
  report("kcomplexFloat",  scalars.kcfloat_  == scalarsRef.kcfloat_);
  report("kcomplexDouble", scalars.kcdouble_ == scalarsRef.kcdouble_);
  return errors;
}
// Analysis function
// This analysis checks that the restart routines are performing as they should
void Analysis(DataBlock& data) {


  idfx::cout << "Analysis: Checking restart routines" << std::endl;

  // Snapshot the scalar values so they can be checked after the round-trip
  scalarsRef = scalars;

  // Fill the replicated arrays and snapshot them for the round-trip check
  InitArrays();
  SnapshotArrays();

  // Trigger dump creation
  // data.SetBoundaries();
  myOutput->ForceWriteDump(data);

    // Mirror data on Host
  DataBlockHost d(data);

  // Sync it
  d.SyncFromDevice();

  // Create local arrays to store the current physical state
  IdefixHostArray4D<real> myVc = IdefixHostArray4D<real>("myVc", d.Vc.extent(0), data.np_tot[KDIR], data.np_tot[JDIR],data.np_tot[IDIR]);
  IdefixHostArray4D<real> myVs = IdefixHostArray4D<real>("myVs", DIMENSIONS, data.np_tot[KDIR]+KOFFSET, data.np_tot[JDIR]+JOFFSET,data.np_tot[IDIR]+IOFFSET);
  #ifdef EVOLVE_VECTOR_POTENTIAL
  IdefixHostArray4D<real> myVe = IdefixHostArray4D<real>("myVe", AX3e+1, data.np_tot[KDIR]+KOFFSET, data.np_tot[JDIR]+JOFFSET,data.np_tot[IDIR]+IOFFSET);
  #endif
  // Transfer the datablock to myVc and myVs
  for(int n = 0; n < d.Vc.extent(0) ; n++) {
    for(int k = 0; k < d.np_tot[KDIR] ; k++) {
      for(int j = 0; j < d.np_tot[JDIR] ; j++) {
        for(int i = 0; i < d.np_tot[IDIR] ; i++) {
          myVc(n,k,j,i) = d.Vc(n,k,j,i);
          d.Vc(n,k,j,i) = 0.0;

        }
      }
    }
  }

  for(int n = 0; n < DIMENSIONS ; n++) {
    for(int k = 0; k < d.np_tot[KDIR] + KOFFSET; k++) {
      for(int j = 0; j < d.np_tot[JDIR] + JOFFSET; j++) {
        for(int i = 0; i < d.np_tot[IDIR] + IOFFSET; i++) {
          myVs(n,k,j,i) = d.Vs(n,k,j,i);
          d.Vs(n,k,j,i) = 0.0;
        }
      }
    }
  }
  #ifdef EVOLVE_VECTOR_POTENTIAL
    for(int n = 0; n < AX3e+1 ; n++) {
      for(int k = 0; k < d.np_tot[KDIR] + KOFFSET; k++) {
        for(int j = 0; j < d.np_tot[JDIR] + JOFFSET; j++) {
          for(int i = 0; i < d.np_tot[IDIR] + IOFFSET; i++) {
            myVe(n,k,j,i) = d.Ve(n,k,j,i);
            d.Ve(n,k,j,i) = 0.0;
          }
        }
      }
    }
  #endif

  // Push our datablockHost to erase everything
  d.SyncToDevice();
  // From this point, the dataBlock is full of zeros

  // Randomise every scalar so that only the dump can restore the original values
  RandomizeScalars();
  // Idem for the replicated arrays
  RandomizeArrays();

  // Load back the restart dump
  myOutput->RestartFromDump(data, outnum);
  data.SetBoundaries();
  #ifdef EVOLVE_VECTOR_POTENTIAL
    data.hydro->emf->ComputeMagFieldFromA(data.hydro->Ve, data.hydro->Vs);
  #endif
  d.SyncFromDevice();

  // increment outnum
  outnum++;
  int errornum;

  errornum = 0;
  idfx::cout << "Analysis: checking consistency" << std::endl;
  // Check that the save/load routines have left everything unchanged.
  for(int n = 0; n < d.Vc.extent(0) ; n++) {
    for(int k = d.beg[KDIR]; k < d.end[KDIR] ; k++) {
      for(int j = d.beg[JDIR]; j < d.end[JDIR] ; j++) {
        for(int i = d.beg[IDIR]; i < d.end[IDIR] ; i++) {
          if(myVc(n,k,j,i) != d.Vc(n,k,j,i)) {
            errornum++;
            idfx::cout << "-----------------------------------------" << std::endl
                       << " Error in Vc at (i,j,k,n) = ( " << i << ", " << j << ", " << k << ", " << n << ")" << std::endl
                       << " Coordinates (x1,x2,x3)   = ( " << d.x[IDIR](i) << ", " << d.x[JDIR](j) << ", " << d.x[KDIR](k) << ")" << std::endl
                       << "Original= " << myVc(n,k,j,i) << " New=" << d.Vc(n,k,j,i) << " diff=" << myVc(n,k,j,i)-d.Vc(n,k,j,i) << std::endl;
          }

        }
      }
    }
  }

    for(int k = d.beg[KDIR]; k < d.end[KDIR] ; k++) {
      for(int j = d.beg[JDIR]; j < d.end[JDIR] ; j++) {
        for(int i = d.beg[IDIR]; i < d.end[IDIR]+IOFFSET ; i++) {
          if(myVs(BX1s,k,j,i) != d.Vs(BX1s,k,j,i)) {
            errornum++;
            idfx::cout << "-----------------------------------------" << std::endl
                       << " Error in Vs(BX1s) at (i,j,k) = ( " << i << ", " << j << ", " << k << ")" << std::endl
                       << " Coordinates (x1,x2,x3)   = ( " << d.x[IDIR](i) << ", " << d.x[JDIR](j) << ", " << d.x[KDIR](k) << ")" << std::endl
                       << "Original= " << myVs(BX1s,k,j,i) << " New=" << d.Vs(BX1s,k,j,i) << " diff=" << myVs(BX1s,k,j,i)-d.Vs(BX1s,k,j,i) << std::endl;

          }

        }
      }
    }
    for(int k = d.beg[KDIR]; k < d.end[KDIR] ; k++) {
      for(int j = d.beg[JDIR]; j < d.end[JDIR]+JOFFSET ; j++) {
        for(int i = d.beg[IDIR]; i < d.end[IDIR] ; i++) {
          if(myVs(BX2s,k,j,i) != d.Vs(BX2s,k,j,i)) {
            errornum++;
            idfx::cout << "-----------------------------------------" << std::endl
                       << " Error in Vs(BX2s) at (i,j,k) = ( " << i << ", " << j << ", " << k << ")" << std::endl
                       << " Coordinates (x1,x2,x3)   = ( " << d.x[IDIR](i) << ", " << d.x[JDIR](j) << ", " << d.x[KDIR](k) << ")" << std::endl;
          }

        }
      }
    }
    for(int k = d.beg[KDIR]; k < d.end[KDIR]+KOFFSET ; k++) {
      for(int j = d.beg[JDIR]; j < d.end[JDIR] ; j++) {
        for(int i = d.beg[IDIR]; i < d.end[IDIR] ; i++) {
          if(myVs(BX3s,k,j,i) != d.Vs(BX3s,k,j,i)) {
            errornum++;
            idfx::cout << "-----------------------------------------" << std::endl
                       << " Error in Vs(BX3s) at (i,j,k) = ( " << i << ", " << j << ", " << k << ")" << std::endl
                       << " Coordinates (x1,x2,x3)   = ( " << d.x[IDIR](i) << ", " << d.x[JDIR](j) << ", " << d.x[KDIR](k) << ")" << std::endl
                       << "Original= " << myVs(BX3s,k,j,i) << " New=" << d.Vs(BX3s,k,j,i) << " diff=" << myVs(BX3s,k,j,i)-d.Vs(BX3s,k,j,i) << std::endl;
          }

        }
      }
    }
#ifdef EVOLVE_VECTOR_POTENTIAL
    for(int k = d.beg[KDIR]; k < d.end[KDIR]+KOFFSET ; k++) {
      for(int j = d.beg[JDIR]; j < d.end[JDIR]+JOFFSET ; j++) {
        for(int i = d.beg[IDIR]; i < d.end[IDIR] ; i++) {
          if(myVe(AX1e,k,j,i) != d.Ve(AX1e,k,j,i)) {
            errornum++;
            idfx::cout << "-----------------------------------------" << std::endl
                       << " Error in Ve(AX1e) at (i,j,k) = ( " << i << ", " << j << ", " << k << ")" << std::endl
                       << " Coordinates (x1,x2,x3)   = ( " << d.x[IDIR](i) << ", " << d.x[JDIR](j) << ", " << d.x[KDIR](k) << ")" << std::endl;

          }

        }
      }
    }
    for(int k = d.beg[KDIR]; k < d.end[KDIR]+KOFFSET ; k++) {
      for(int j = d.beg[JDIR]; j < d.end[JDIR] ; j++) {
        for(int i = d.beg[IDIR]; i < d.end[IDIR]+IOFFSET ; i++) {
          if(myVe(AX2e,k,j,i) != d.Ve(AX2e,k,j,i)) {
            errornum++;
            idfx::cout << "-----------------------------------------" << std::endl
                       << " Error in Ve(AX2e) at (i,j,k) = ( " << i << ", " << j << ", " << k << ")" << std::endl
                       << " Coordinates (x1,x2,x3)   = ( " << d.x[IDIR](i) << ", " << d.x[JDIR](j) << ", " << d.x[KDIR](k) << ")" << std::endl;
          }

        }
      }
    }
    for(int k = d.beg[KDIR]; k < d.end[KDIR] ; k++) {
      for(int j = d.beg[JDIR]; j < d.end[JDIR]+JOFFSET ; j++) {
        for(int i = d.beg[IDIR]; i < d.end[IDIR]+IOFFSET ; i++) {
          if(myVe(AX3e,k,j,i) != d.Ve(AX3e,k,j,i)) {
            errornum++;
            idfx::cout << "-----------------------------------------" << std::endl
                       << " Error in Ve(AX3e) at (i,j,k) = ( " << i << ", " << j << ", " << k << ")" << std::endl
                       << " Coordinates (x1,x2,x3)   = ( " << d.x[IDIR](i) << ", " << d.x[JDIR](j) << ", " << d.x[KDIR](k) << ")" << std::endl
                       << "Original= " << myVs(BX3s,k,j,i) << " New=" << d.Vs(BX3s,k,j,i) << " diff=" << myVs(BX3s,k,j,i)-d.Vs(BX3s,k,j,i) << std::endl;
          }

        }
      }
    }
#endif

  // Validate that every scalar survived the dump/restart round-trip
  errornum += CheckScalars();

  // Validate that every replicated array survived the dump/restart round-trip
  errornum += CheckArrays();

  idfx::cout << "Analysis: consistency check done with " << errornum << " errors " << std::endl;
  if(errornum>0) {
    IDEFIX_ERROR("Restart from dump failed validation");
  }
}

// Release the Kokkos views before Kokkos::finalize() is called
Setup::~Setup() {
  arrRefs = ArrayRefs{};
  testArrays.reset();
}

// Initialisation routine. Can be used to allocate
// Arrays or variables which are used later on
Setup::Setup(Input &input, Grid &grid, DataBlock &data, Output &output) {
   // Kokkos views must be allocated after Kokkos::initialize()
   testArrays = std::make_unique<TestArrays>();
   if(input.CheckEntry("Output","analysis")>0) {
     output.EnrollAnalysis(&Analysis);
     myOutput = &output;
     outnum=0;
   }
   // Register every supported scalar type so the dump I/O exercises them all.
   // NOTE: names are truncated to NAMESIZE-1 = 15 characters in the dump file,
   // so keep them short enough to match exactly on read-back.
   data.dump->RegisterVariable(&scalars.int_,      "testInt");
   data.dump->RegisterVariable(&scalars.int64_,    "testInt64");
   data.dump->RegisterVariable(&scalars.float_,    "testFloat");
   data.dump->RegisterVariable(&scalars.double_,   "testDouble");
   data.dump->RegisterVariable(&scalars.bool_,     "testBool");
   data.dump->RegisterVariable(&scalars.cfloat_,   "testCFloat");
   data.dump->RegisterVariable(&scalars.cdouble_,  "testCDouble");
   data.dump->RegisterVariable(&scalars.kcfloat_,  "testKCFloat");
   data.dump->RegisterVariable(&scalars.kcdouble_, "testKCDouble");

   // Replicated 1D/2D arrays: one field per data type and rank. Names must stay
   // within NAMESIZE-1 = 15 characters to match exactly on read-back.
   data.dump->RegisterVariable(testArrays->int1D,    "testInt1D");
   data.dump->RegisterVariable(testArrays->int641D,  "testInt641D");
   data.dump->RegisterVariable(testArrays->float1D,  "testFloat1D");
   data.dump->RegisterVariable(testArrays->double1D, "testDouble1D");
   data.dump->RegisterVariable(testArrays->bool1D,   "testBool1D");
   data.dump->RegisterVariable(testArrays->double2D, "testDouble2D");
   data.dump->RegisterVariable(testArrays->bool2D,   "testBool2D");
   data.dump->RegisterVariable(testArrays->int2D,    "testInt2D");
   data.dump->RegisterVariable(testArrays->cfloat2D, "testCFloat2D");
   data.dump->RegisterVariable(testArrays->cdouble2D, "testCDouble2D");
}

// This routine initialize the flow
// Note that data is on the device.
// One can therefore define locally
// a datahost and sync it, if needed
void Setup::InitFlow(DataBlock &data) {
    // Create a host copy
    DataBlockHost d(data);
    real x,y,z;
    IdefixHostArray4D<real> Ve;

    #ifndef EVOLVE_VECTOR_POTENTIAL
    Ve = IdefixHostArray4D<real>("Potential vector",3, d.np_tot[KDIR]+1, d.np_tot[JDIR]+1, d.np_tot[IDIR]+1);
    #else
    Ve = d.Ve;
    #endif

    bool haveTracer = data.hydro->haveTracer;

    real B0=1.0/sqrt(4.0*M_PI);

    for(int k = 0; k < d.np_tot[KDIR] ; k++) {
        for(int j = 0; j < d.np_tot[JDIR] ; j++) {
            for(int i = 0; i < d.np_tot[IDIR] ; i++) {
                x=d.x[IDIR](i);
                y=d.x[JDIR](j);
                z=d.x[KDIR](k);

                d.Vc(RHO,k,j,i) = 25.0/(36.0*M_PI);
                d.Vc(PRS,k,j,i) = 5.0/(12.0*M_PI);
                d.Vc(VX1,k,j,i) = -sin(2.0*M_PI*y);
                d.Vc(VX2,k,j,i) = sin(2.0*M_PI*x)+cos(2.0*M_PI*z);
                d.Vc(VX3,k,j,i) = cos(2.0*M_PI*x);

                real xl=d.xl[IDIR](i);
                real yl=d.xl[JDIR](j);
                real zl=d.xl[KDIR](k);
                Ve(IDIR,k,j,i) = B0/(2.0*M_PI)*(cos(2.0*M_PI*yl));
                Ve(JDIR,k,j,i) = B0/(2.0*M_PI)*sin(2.0*M_PI*xl);
                Ve(KDIR,k,j,i) = B0/(2.0*M_PI)*(
                                    cos(2.0*M_PI*yl) + cos(4.0*M_PI*xl)/2.0);

                if(haveTracer) {
                  d.Vc(TRG  ,k,j,i) = x>0.5?  1.0:0.0;
                  d.Vc(TRG+1,k,j,i) = z>0.5?  1.0:0.0;
                }
            }
        }
    }

    #ifndef EVOLVE_VECTOR_POTENTIAL
    d.MakeVsFromAmag(Ve);
    #endif
    // Send it all, if needed
    d.SyncToDevice();
}

// Analyse data to produce an output
void MakeAnalysis(DataBlock & data) {

}



// Do a specifically designed user step in the middle of the integration
void ComputeUserStep(DataBlock &data, real t, real dt) {

}
