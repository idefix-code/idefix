// ***********************************************************************************
// Idefix MHD astrophysical code
//
// Source file src/unitTestApp.hpp
//
// Last modified : 11/2026
//
// Copyright(C) by :
// - Sébastien Valat <sebastien.valat@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2026)
// and other code contributors
//
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef UNITTESTAPP_HPP_
#define UNITTESTAPP_HPP_

//std c
#include <cstdio>
//std c++
#include <list>
#include <string>
#include <memory>
//idefix
#include "idefix.hpp"
#include "dataBlock.hpp"
#include "dataBlockHost.hpp"
#include "grid.hpp"
#include "gridHost.hpp"

namespace idfx {

/**
 * @brief Create a temporary text file.
 *
 * Create a temporary texte file, store the given data in, and delete
 * it when destroyed. This is used for unit tests.
 */
class TempConfigFile {
 public:
  explicit TempConfigFile(const std::list<std::string> & lines);
  ~TempConfigFile(void);
  const std::string & getFilePath(void) const;
 private:
  std::FILE * fp{nullptr};
  std::string fname;
};

/**
 * @brief Handle all the components of Idefix for unit tests.
 *
 * Construct all the Idefix components except physics to perform
 * unit tests on the Idefix classes.
 */
class UnitTestApp {
 public:
  UnitTestApp(size_t s1, size_t s2, size_t s3, float cell_size = 1.0);
 public:
  std::shared_ptr<Input> input;
  std::shared_ptr<Grid> grid;
  std::shared_ptr<GridHost> gridHost;
  std::shared_ptr<DataBlock> dataBlock;
  std::shared_ptr<DataBlockHost> dataBlockHost;
 private:
  void InitInput(size_t s1, size_t s2, size_t s3, float cell_size);
  void InitComponents(void);
};

} // namespace idfx

#endif //UNITTESTAPP_HPP_
