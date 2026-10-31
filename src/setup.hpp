// ***********************************************************************************
// Idefix MHD astrophysical code
//
// Source file src/setup.hpp
//
// Last modified : 10/2023
//
// Copyright(C) by :
// - Geoffroy Lesur <geoffroy.lesur@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2020 - 2023)
// - Soufiane Baghdadi <soufiane.baghdadi@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2020)
// - Clément Robert <clement.robert@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2021)
// and other code contributors
//
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef SETUP_HPP_
#define SETUP_HPP_

#include <memory>
#include "idefix.hpp"
#include "input.hpp"
#include "grid.hpp"
#include "dataBlock.hpp"
#include "fluid.hpp"
#include "output.hpp"

// These two will likely be used in the code, so we include them here
#include "gridHost.hpp"
#include "dataBlockHost.hpp"
#include "boundary.hpp"

// Forward declarations
class Setup;

// Base class for user-defined setup with OOP style
class SetupBase {
 public:
    SetupBase(Input &input, Grid &grid, DataBlock &data, Output &output)
        : myInput(input), myGrid(grid), myData(data), myOutput(output) {}

    virtual ~SetupBase() = default;
    virtual void InitFlow(DataBlock &data) = 0;

 protected:
    Input &myInput;
    Grid &myGrid;
    DataBlock &myData;
    Output &myOutput;
};

// Default Setup class for backward compatibility (unchanged)
class Setup {
 public:
    Setup(Input &, Grid &, DataBlock &, Output&);
    ~Setup();
    void InitFlow(DataBlock &);
};

// Default Setup wrapper that inherits from SetupBase
// This allows the factory to return SetupBase* for the default case
class DefaultSetup : public SetupBase {
 public:
    DefaultSetup(Input &input, Grid &grid, DataBlock &data, Output &output)
        : SetupBase(input, grid, data, output), mySetup(input, grid, data, output) {}

    void InitFlow(DataBlock &data) override {
        mySetup.InitFlow(data);
    }

 private:
    Setup mySetup;
};

#endif // SETUP_HPP_
