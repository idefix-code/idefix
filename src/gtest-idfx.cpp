// ***********************************************************************************
// Idefix MHD astrophysical code
//
// Source file src/utils/gtest-idfx.hpp
//
// Last modified : 10/2026
//
// Copyright(C) by :
// - Sébastien Valat <sebastien.valat@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2026)
// and other code contributors
//
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

// includes
#include <string>
#if __has_include(<filesystem>)
  #include <filesystem> // NOLINT [build/c++17]
  namespace fs = std::filesystem;
#elif __has_include(<experimental/filesystem>)
  #include <experimental/filesystem>
  namespace fs = std::experimental::filesystem;
#else
  error "Missing the <filesystem> header."
#endif
#include <gtest/gtest.h>
#include <gmock/gmock.h>
#ifdef WITH_MPI
  #include <mpi.h>
#endif //WITH_MPI
#include "idefix.hpp"
#include "gtest-idfx.hpp"

/**
 * Get the full path to get access to the data/ directory on same location than the sourceFile.
 * @param sourceFile Absolute path to the source file to consider.
 * @param fname Filename to append to SRC_DIR/data/.
*/
std::string idfx_get_data_file_path(
  const std::string & sourceFile, const std::string & fname) {
  fs::path p(sourceFile);
  return p.parent_path() / "data" / fname;
}

/**
 * Configure MPI environnement to setup/tear down
**/
class MPIEnvironment : public ::testing::Environment {
 public:
  MPIEnvironment(int * argc, char *** argv): argc(argc), argv(argv) {}
  ~MPIEnvironment() override {}

  // Override this to define how to set up the environment.
  void SetUp() override {
    #ifdef WITH_MPI
      int mpiError = MPI_Init(argc, argv);
      ASSERT_EQ(mpiError, 0);
    #endif
  }

  // Override this to define how to tear down the environment.
  void TearDown() override {
    #ifdef WITH_MPI
      int mpiError = MPI_Finalize();
      ASSERT_EQ(mpiError, 0);
    #endif
  }

 private:
  int * argc{nullptr};
  char *** argv{nullptr};
};

/**
 * configure MPI environnement to setup/tear down
*/
class KokkosEnvironment : public ::testing::Environment {
 public:
  KokkosEnvironment(int argc, char ** argv): argc(argc), argv(argv) {}
  ~KokkosEnvironment() override {}

  // Override this to define how to set up the environment.
  void SetUp() override {
    Kokkos::initialize(this->argc, this->argv);
  }

  // Override this to define how to tear down the environment.
  void TearDown() override {
    Kokkos::finalize();
  }

 private:
  int argc{0};
  char ** argv{nullptr};
};

/**
 * configure MPI environnement to setup/tear down
*/
class IdefixEnvironment : public ::testing::Environment {
 public:
  ~IdefixEnvironment() override {}

  // Override this to define how to set up the environment.
  void SetUp() override {
    int status = idfx::initialize();
    ASSERT_EQ(status, 0);
  }

  // Override this to define how to tear down the environment.
  void TearDown() override {
  }
};

// setup advanded main for tests using MPI / Kokkos
int main(int argc, char* argv[]) {
  bool initKokkosBeforeMPI = false;
  // When running on GPUS with Omnipath network,
  // Kokkos needs to be initialised *before* the MPI layer
  #ifdef KOKKOS_ENABLE_CUDA
    if(std::getenv("PSM2_CUDA") != NULL) {
      initKokkosBeforeMPI = true;
    }
  #endif

  //ini gtest
  ::testing::InitGoogleTest(&argc, argv);

  //init MPI & Kokkos with right order
  if (initKokkosBeforeMPI) {
    ::testing::AddGlobalTestEnvironment(new MPIEnvironment(&argc, &argv));
    ::testing::AddGlobalTestEnvironment(new KokkosEnvironment(argc, argv));
  } else {
    ::testing::AddGlobalTestEnvironment(new KokkosEnvironment(argc, argv));
    ::testing::AddGlobalTestEnvironment(new MPIEnvironment(&argc, &argv));
  }

  //init idefix
  ::testing::AddGlobalTestEnvironment(new IdefixEnvironment);

  //run tests
  return RUN_ALL_TESTS();
}
