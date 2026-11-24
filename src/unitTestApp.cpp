// ***********************************************************************************
// Idefix MHD astrophysical code
//
// Source file src/unitTestApp.cpp
//
// Last modified : 11/2026
//
// Copyright(C) by :
// - Sébastien Valat <sebastien.valat@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2026)
// and other code contributors
//
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include "unitTestApp.hpp"
//unix
#include <unistd.h>
//std c
#include <cstdio>
//std c++
#include <list>
#include <memory>
#include <string>
#include <format>
//idefix
#include "units.hpp"
#include "fluid.hpp"

namespace idfx {

TempConfigFile::TempConfigFile(const std::list<std::string> & lines) {
  //tmp name
  char tmpname[] = "/tmp/idefix-ut-config-XXXXXXX";

  //open
  this->fp = fdopen(mkstemp(tmpname), "w");
  this->fname = tmpname;

  //write
  for (const auto & it : lines) {
    fputs(it.c_str(), this->fp);
    fputs("\n", this->fp);
  }

  //flush
  fflush(this->fp);
}

TempConfigFile::~TempConfigFile(void) {
  fclose(this->fp);
  unlink(this->fname.c_str());
}

const std::string & TempConfigFile::getFilePath(void) const {
  return this->fname;
}

UnitTestApp::UnitTestApp(size_t s1, size_t s2, size_t s3, float cell_size) {
  //init input and components
  this->InitInput(s1, s2, s3, cell_size);
  this->InitComponents();
}

void UnitTestApp::InitInput(size_t s1, size_t s2, size_t s3, float cell_size) {
  //create config file
  TempConfigFile configFile({
    "[Grid]",
    std::format("X1-grid    1  1.0  {}  u  {}", s1, 1 + static_cast<float>(s1)*cell_size), //R
    std::format("X2-grid    1  0.0  {}  u  3.1415926", s2, static_cast<float>(s2)*cell_size), //teta
    std::format("X3-grid    1  0.0  {}  u  {}", s3, 2 * 3.1415926), //phi
    "",
    "[Boundary]",
    "X1-beg    periodic",
    "X1-end    periodic",
    "X2-beg    periodic",
    "X2-end    periodic",
    "X3-beg    periodic",
    "X3-end    periodic",
  });

  //build args
  char* argc[] = {strdup("idefix"), strdup("-i"), strdup(configFile.getFilePath().c_str()), nullptr};

  //load it as Input object
  this->input = std::make_shared<Input>(3, argc);

  //free args
  free(argc[0]);
  free(argc[1]);
  free(argc[2]);
}

void UnitTestApp::InitComponents() {
  //ini components
  idfx::units.Init(*input);

  // Allocate the grid on device
  this->grid = std::make_shared<Grid>(*input);

  // Allocate the grid image on host
  this->gridHost = std::make_shared<GridHost>(*grid);

  // Actually make the grid on host and sync it on the device
  this->gridHost->MakeGrid(*input);
  this->gridHost->SyncToDevice();

  // instantiate required objects.
  this->dataBlock = std::make_shared<DataBlock>(*grid, *input, false);
}

} // namespace idfx
