// ********************************************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ********************************************************************************************************


#include "../idefix.hpp"
#include "hydro.hpp"




Hydro::Hydro(Input &input, Grid &grid) {
    idfx::pushRegion("Hydro::Hydro(input)");

    if(input.CheckEntry("Hydro","gamma")>0) this->gamma = input.GetReal("Hydro","gamma",0);
    else this->gamma = 5.0/3.0;

    if(input.CheckEntry("Hydro","csiso")>0) {
        real cs = input.GetReal("Hydro","csiso",0);
        this->C2Iso = cs*cs;
    }
    else this->C2Iso = 1.0;
    
    // read Solver from input file
    std::string solverString = input.GetString("Hydro","Solver",0);
    
    if (solverString.compare("tvdlf") == 0)     mySolver = TVDLF;
    else if (solverString.compare("hll") == 0)  mySolver = HLL;
    #if MHD == YES
        else if (solverString.compare("hlld") == 0) mySolver = HLLD;
    #else
        else if (solverString.compare("hllc") == 0) mySolver = HLLC;
    #endif
    else if (solverString.compare("roe") == 0)  mySolver = ROE;
    else {
        std::stringstream msg;
    #if MHD == YES
        msg << "Unknown MHD solver type " << solverString;
    #else
        msg << "Unknown HD solver type " << solverString;
    #endif
        IDEFIX_ERROR(msg);
    }
    
    // No userdefBoundary by default
    this->haveUserDefBoundary = false;
    this->haveInternalBoundary = false;

    // Source terms (always activated when non-cartesian geometry because of curvature source terms)
    #if GEOMETRY == CARTESIAN
    this->haveSourceTerms = false;
    #else
    this->haveSourceTerms = true;
    #endif
    this->haveUserSourceTerm = false;

    // Check whether we have rotation
    int rotation = input.CheckEntry("Hydro","Rotation");

    if(rotation>=0 ) {
        this->haveSourceTerms = true;
        this->haveRotation = true;
        if(rotation != 3) IDEFIX_ERROR("Rotation needs a 3 components vector in idefix.ini");
        this->OmegaX1 = input.GetReal("Hydro","Rotation",0);
        this->OmegaX2 = input.GetReal("Hydro","Rotation",1);
        this->OmegaX3 = input.GetReal("Hydro","Rotation",2);

        idfx::cout << "Hydro: Rotation enabled with Omega=(" << this->OmegaX1 << ", " << this->OmegaX2 << ", " << this->OmegaX3 << ")" << std::endl;
    }
    else {
        this->haveRotation = false;
    }

    // Check whether we have shearing box
    int shearingbox = input.CheckEntry("Hydro","ShearingBox");

    if(shearingbox>=0 ) {
        this->haveShearingBox = true;
        this->haveSourceTerms = true;
        if(shearingbox != 1) IDEFIX_ERROR("Shearing box needs a scalar value for the shear rate in idefix.ini");
        this->sbS = input.GetReal("Hydro","ShearingBox",0);

        // Get box size
        this->sbLx = grid.xend[IDIR] - grid.xbeg[IDIR];

        idfx::cout << "Hydro: ShearingBox enabled with Shear rate= " << this->sbS <<  "and Lx= " << sbLx << std::endl;
    }
    else {
        this->haveShearingBox = false;
    }

    // Gravitational potential
    this->haveGravPotential = false;
    this->gravPotentialFunc = nullptr;
    int gravPotential = input.CheckEntry("Hydro","GravPotential");
    if(gravPotential>=0) {
        std::string potentialString = input.GetString("Hydro","GravPotential",0);
        if(potentialString.compare("userdef") == 0) {
            this->haveGravPotential = true;        
            idfx::cout << "Hydro:: Enabling user-defined gravitational potential" << std::endl;
        }
        else {
            IDEFIX_ERROR("Unknown type of gravitational potential in idefix.ini. Only userdef is implemented");
        }
    }

    // Parabolic term
    haveParabolicTerms = false;

    // Nonideal MHD
    haveResistivity = Disabled;
    haveHall = Disabled;
    haveAmbipolar = Disabled;

    ohmicDiffusivityFunc = NULL;
    ambipolarDiffusivityFunc = NULL;
    hallDiffusivityFunc = NULL;

    this->needCurrent = false;

    #if MHD == YES
    if(input.CheckEntry("Hydro","Resistivity")>=0 || 
       input.CheckEntry("Hydro","Ambipolar")>=0 ||
       input.CheckEntry("Hydro","Hall")>=0 ) {
           this->needCurrent = true;
           
           if(input.CheckEntry("Hydro","Resistivity")>=0) {
               if(input.GetString("Hydro","Resistivity",0).compare("constant") == 0) {
                   idfx::cout << "Hydro: Enabling Ohmic resistivity with constant diffusivity." << std::endl;
                   this->etaO = input.GetReal("Hydro","Resistivity",1);
                   this->haveParabolicTerms = true;
                   this->haveResistivity = Constant;
               }
               else if(input.GetString("Hydro","Resistivity",0).compare("userdef") == 0) {
                   idfx::cout << "Hydro: Enabling Ohmic resistivity with user-defined diffusivity function." << std::endl;
                   this->haveParabolicTerms = true;
                   this->haveResistivity = UserDefFunction;
               }
               else {
                   IDEFIX_ERROR("Unknown resistivity definition in idefix.ini. Can only be constant or userdef.");
               }
           }
           if(input.CheckEntry("Hydro","Ambipolar")>=0) {
               if(input.GetString("Hydro","Ambipolar",0).compare("constant") == 0) {
                   idfx::cout << "Hydro: Enabling ambipolar diffusion with constant diffusivity." << std::endl;
                   this->xA = input.GetReal("Hydro","Ambipolar",1);
                   this->haveParabolicTerms = true;
                   this->haveAmbipolar = Constant;
               }
               else if(input.GetString("Hydro","Ambipolar",0).compare("userdef") == 0) {
                   idfx::cout << "Hydro: Enabling ambipolar diffusion with user-defined diffusivity function." << std::endl;
                   this->haveParabolicTerms = true;
                   this->haveAmbipolar = UserDefFunction;
               }
               else {
                   IDEFIX_ERROR("Unknown ambipolar definition in idefix.ini. Can only be constant or userdef.");
               }
           }
           if(input.CheckEntry("Hydro","Hall")>=0) {
               // Check consistency
               if(mySolver != HLL ) IDEFIX_ERROR("Hall effect is only compatible with HLL Riemann solver.");
               #if EMF_AVERAGE != ARITHMETIC
                    IDEFIX_ERROR("the Hall effect module is demonstrated stable only when using EMF_AVERAGE=ARITHMETIC");
                #endif
               if(input.GetString("Hydro","Hall",0).compare("constant") == 0) {
                   idfx::cout << "Hydro: Enabling Hall effect with constant diffusivity." << std::endl;
                   this->xH = input.GetReal("Hydro","Hall",1);
                   this->haveHall = Constant;
               }
               else if(input.GetString("Hydro","Hall",0).compare("userdef") == 0) {
                   idfx::cout << "Hydro: Enabling Hall effect with user-defined diffusivity function." << std::endl;
                   this->haveHall = UserDefFunction;
               }
               else {
                   IDEFIX_ERROR("Unknown Hall definition in idefix.ini. Can only be constant or userdef.");
               }
           }
       }
    #endif


    idfx::popRegion();
}

void Hydro::EnrollUserDefBoundary(UserDefBoundaryFunc myFunc) {
    this->userDefBoundaryFunc = myFunc;
    this->haveUserDefBoundary = true;
    idfx::cout << "Hydro: User-defined boundary condition has been enrolled" << std::endl;
}

void Hydro::EnrollInternalBoundary(InternalBoundaryFunc myFunc) {
    this->internalBoundaryFunc = myFunc;
    this->haveInternalBoundary = true;
    idfx::cout << "Hydro: User-defined internal boundary condition has been enrolled" << std::endl;
}

void Hydro::EnrollGravPotential(GravPotentialFunc myFunc) {
    if(!this->haveGravPotential) IDEFIX_ERROR("In order to enroll your gravitational potential, you need to enable it first in the .ini file.");
    this->gravPotentialFunc = myFunc;
    idfx::cout << "Hydro: User-defined gravitational potential has been enrolled" << std::endl;
}

void Hydro::EnrollUserSourceTerm(SrcTermFunc myFunc) {
    this->userSourceTerm = myFunc;
    this->haveUserSourceTerm = true;
    this->haveSourceTerms = true;
    idfx::cout << "Hydro: User-defined source term has been enrolled" << std::endl;
}

void Hydro::EnrollOhmicDiffusivity(DiffusivityFunc myFunc) {
    if(this->haveResistivity < UserDefFunction) IDEFIX_ERROR("Ohmic diffusivity enrollment requires Hydro/Resistivity to be set to userdef in .ini file");
    this->ohmicDiffusivityFunc = myFunc;
    idfx::cout << "Hydro: User-defined ohmic diffusivity has been enrolled" << std::endl;
}

void Hydro::EnrollAmbipolarDiffusivity(DiffusivityFunc myFunc) {
    if(this->haveAmbipolar < UserDefFunction) IDEFIX_ERROR("Ambipolar diffusivity enrollment requires Hydro/Ambipolar to be set to userdef in .ini file");
    this->ambipolarDiffusivityFunc = myFunc;
    idfx::cout << "Hydro: User-defined ambipolar diffusivity has been enrolled" << std::endl;
}

void Hydro::EnrollHallDiffusivity(DiffusivityFunc myFunc) {
    if(this->haveHall < UserDefFunction) IDEFIX_ERROR("Hall diffusivity enrollment requires Hydro/Hall to be set to userdef in .ini file");
    this->hallDiffusivityFunc = myFunc;
    idfx::cout << "Hydro: User-defined Hall diffusivity has been enrolled" << std::endl;
}

Hydro::Hydro() {

}

real Hydro::GetGamma() {
    return(this->gamma);
}

void Hydro::SetGamma(real newGamma) {
    this->gamma=newGamma;
}

real Hydro::GetC2iso() {
    return(this->C2Iso);
}





