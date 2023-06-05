// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2021 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef DATABLOCK_PLANETARYSYSTEM_PLANETARYSYSTEM_HPP_
#define DATABLOCK_PLANETARYSYSTEM_PLANETARYSYSTEM_HPP_

#include <vector>
#include "idefix.hpp"
#include "input.hpp"
#include "planet.hpp"

// forward class declaration
class DataBlock;



class PlanetarySystem {
 public:
    enum Integrator {RK4=1, ANALYTICAL, RK5};
    enum SmoothingFunction {PLUMMER=1, POLYNOMIAL};

    PlanetarySystem(Input&, DataBlock*);
    void EvolveSystem(DataBlock&, const real& );
    void IntegrateAnalytically(DataBlock&, const real&);
    void IntegrateRK4(DataBlock&, const real&);
    void IntegrateRK5(DataBlock&, const real&);
    void ShowConfig();
    void AddPlanetsPotential(IdefixArray3D<real> &, real);
    std::vector<PointSpeed> ComputeRHS(real&, std::vector<Planet>);

    // number of planets
    int nbp{0};
    real torqueNormalization{ONE_F};
    std::vector<Planet> planet;
    real GetSmoothingValue() const;
    real GetSmoothingExponent() const;

 protected:
    void AdvancePlanetFromDisk(DataBlock&, const real&);
    void IntegratePlanets(DataBlock&, const real&);
    friend class Planet;
    real massTaper{ZERO_F};
    real smoothingValue;
    real smoothingExponent;
    bool excludeHill{false};
    bool indirectPlanetsTerm{true};
    bool feelDisk;
    bool feelPlanets;
    bool halfdisk;
    Integrator myPlanetaryIntegrator;
    SmoothingFunction myPlanetarySmoothing;
    DataBlock *data;
};

#endif // DATABLOCK_PLANETARYSYSTEM_PLANETARYSYSTEM_HPP_
