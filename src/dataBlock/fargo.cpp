// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include <algorithm>
#include <string>
#include <vector>

#include "idefix.hpp"
#include "fluid.hpp"
#include "dataBlock.hpp"
#include "fargo.hpp"



Fargo::Fargo(Input &input, int nmax, DataBlock *data) {
  idfx::pushRegion("Fargo::Init");
  this->data = data;

  // A bit of arithmetic to get the sizes of the working array
  this->nghost = data->nghost;
  this->beg = data->beg;
  this->end = data->end;

  if(input.CheckBlock("Fargo")) {
    std::string opType = input.Get<std::string>("Fargo","velocity",0);
    if(opType.compare("userdef")==0) {
       #if GEOMETRY != SPHERICAL && GEOMETRY != POLAR
        IDEFIX_ERROR("Fargo+userdef is only compatible with SPHERICAL and POLAR geometries");
      #endif
      this->type=userdef;
    } else if(opType.compare("shearingbox")==0) {
      this->type=shearingbox;
      #if GEOMETRY != CARTESIAN
        // Actually, this has never really been tested, so assumes it doesn't work...
        IDEFIX_ERROR("Fargo+shearingbox is only compatible with cartesian geometry");
      #endif
    } else {
      IDEFIX_ERROR("Unknown fargo velocity in the input file. "
      "Only userdef and shearingbox are allowed");
    }
    this->maxShift = input.GetOrSet<int>("Fargo", "maxShift",0, 10);
  } else {
    // DEPRECATED: initialisation from the [Hydro] block
    if(input.CheckEntry("Hydro","fargo")>=0) {
      std::string opType = input.Get<std::string>("Hydro","fargo",0);
      if(opType.compare("userdef")==0) {
        this->type=userdef;
      } else if(opType.compare("shearingbox")==0) {
        this->type=shearingbox;
        #if GEOMETRY != CARTESIAN
          IDEFIX_ERROR("Fargo+shearingbox only compatible with cartesian geometry");
        #endif
      } else {
        IDEFIX_ERROR("Unknown fargo option in the input file. "
        "Only userdef and shearingbox are allowed");
      }
    } else {
      IDEFIX_ERROR("Fargo should be enabled in the input file with either userdef "
                    "or shearingbox options");
    }
  }

  #if GEOMETRY == CARTESIAN || GEOMETRY == POLAR
    #if DIMENSIONS < 2
      IDEFIX_ERROR("Fargo should be used with DIMENSIONS >= 2 in Cartesian or Polar geometries");
    #endif
    // Check if there is a domain decomposition in the intended fargo direction
    if(data->mygrid->nproc[JDIR]>1) {
      haveDomainDecomposition = true;
      this->nghost[JDIR] += this->maxShift;
      this->beg[JDIR] += this->maxShift;
      this->end[JDIR] += this->maxShift;
      if(data->np_int[JDIR] < this->maxShift + data->nghost[JDIR]) {
        IDEFIX_ERROR("Subdomain size < Fargo:maxShift. "
                     "Try reducting the number of processes along X2");
      }
    }
    if(this->type==userdef)
      this->meanVelocity = IdefixArray2D<real>("FargoVelocity",data->np_tot[KDIR],
                                                             data->np_tot[IDIR]);
  #elif GEOMETRY == SPHERICAL
    #if DIMENSIONS < 3
      IDEFIX_ERROR("Fargo should be used with DIMENSIONS == 3 in Spherical geometry");
    #endif
    // Check if there is a domain decomposition in the intended fargo direction
    if(data->mygrid->nproc[KDIR]>1) {
      haveDomainDecomposition = true;
      this->nghost[KDIR] += this->maxShift;
      this->beg[KDIR] += this->maxShift;
      this->end[KDIR] += this->maxShift;
      if(data->np_int[KDIR] < this->maxShift + data->nghost[KDIR]) {
        IDEFIX_ERROR("Subdomain size < Fargo:maxShift. "
                     "Try reducting the number of processes along X3");
      }
    }
    this->meanVelocity = IdefixArray2D<real>("FargoVelocity",data->np_tot[JDIR],
                                                             data->np_tot[IDIR]);
  #else
    IDEFIX_ERROR("Fargo is not compatible with the GEOMETRY you intend to use");
  #endif


  // Initialise our scratch space
  // Maximum number of variables
  int nvar = data->hydro->Vc.extent(0);
  if(data->haveDust) {
    for(int n = 0 ; n < data->dust.size() ; n++) {
      nvar = std::max(nvar,static_cast<int>(data->dust[n]->Vc.extent(0)));
    }
  }

  this->scrhUc = IdefixArray4D<real>("FargoVcScratchSpace",nvar
                                      ,end[KDIR]-beg[KDIR] + 2*nghost[KDIR]
                                      ,end[JDIR]-beg[JDIR] + 2*nghost[JDIR]
                                      ,end[IDIR]-beg[IDIR] + 2*nghost[IDIR]);

  #if MHD == YES
    if(haveDomainDecomposition) {
      this->scrhVs = IdefixArray4D<real>("FargoVsScratchSpace",DIMENSIONS
                                          ,end[KDIR]-beg[KDIR] + 2*nghost[KDIR]+KOFFSET
                                          ,end[JDIR]-beg[JDIR] + 2*nghost[JDIR]+JOFFSET
                                          ,end[IDIR]-beg[IDIR] + 2*nghost[IDIR]+IOFFSET);


    } else {
      // A separate allocation for scrhVs is only needed with domain decomposition, otherwise,
      // we just make a reference to scrhVs
      this->scrhVs = data->hydro->Vs;
    }
  #endif
  #ifdef WITH_MPI
    if(haveDomainDecomposition) {
      std::vector<int> vars;
      for(int i=0 ; i < nvar ; i++) {
        vars.push_back(i);
      }
      #if MHD == YES
        this->mpi.Init(data->mygrid, vars, this->nghost.data(), data->np_int.data(), true);
      #else
        this->mpi.Init(data->mygrid, vars, this->nghost.data(), data->np_int.data());
      #endif
    }
  #endif


  idfx::popRegion();
}

void Fargo::ShowConfig() {
  idfx::pushRegion("Fargo::ShowConfig");
  if(type==userdef) {
    idfx::cout << "Fargo: ENABLED with user-defined velocity function." << std::endl;
    if(!fargoVelocityFunc) {
      IDEFIX_ERROR("No Fargo velocity function has been enabled.");
    }
  } else if(type==shearingbox) {
    idfx::cout << "Fargo: ENABLED with shearing-box velocity function." << std::endl;
  } else {
    IDEFIX_ERROR("Something went wrong during Fargo initialisation");
  }
  #ifdef HIGH_ORDER_FARGO
    idfx::cout << "Fargo: using high order PPM advection scheme." << std::endl;
  #else
    idfx::cout << "Fargo: using standard PLM advection scheme." << std::endl;
  #endif
  if(haveDomainDecomposition) {
    idfx::cout << "Fargo: using domain decomposition along the azimuthal direction"
               << " with maxShift=" << this->maxShift << std::endl;
  }
  idfx::popRegion();
}

void Fargo::EnrollVelocity(FargoVelocityFunc myFunc) {
  if(this->type!=userdef) {
    IDEFIX_WARNING("Fargo velocity function enrollment requires Hydro/Fargo "
                 "to be set to userdef in .ini file");
  }
  this->fargoVelocityFunc = myFunc;
}

// This function fetches Fargo velocity when required
void Fargo::GetFargoVelocity(real t) {
  idfx::pushRegion("Fargo::GetFargoVelocity");
  if(type != userdef)
    IDEFIX_ERROR("Fargo::GetFargoVelocity should not be called when fargo != userdef");
  if(velocityHasBeenComputed == false) {
    if(fargoVelocityFunc== NULL) {
      IDEFIX_ERROR("No Fargo velocity function has been defined");
    }
    fargoVelocityFunc(*data, meanVelocity);
    velocityHasBeenComputed = true;
    if(this->haveDomainDecomposition) {
      CheckMaxDisplacement();
    }
  }
  idfx::popRegion();
}

// This function checks that the velocity provided does not lead to an azimuthal displacement
// larger than the one allowed
void Fargo::CheckMaxDisplacement() {
    IdefixArray2D<real> meanV = this->meanVelocity;
    int ibeg, iend, jbeg, jend;
    IdefixArray1D<real> xi;
    IdefixArray1D<real> xj;
    IdefixArray1D<real> dxk;
    [[maybe_unused]] FargoType fargoType = type;
    [[maybe_unused]] real sbS = data->hydro->sbS;
    real invDt = 0;

    // Get domain size
    #if GEOMETRY == CARTESIAN || GEOMETRY == POLAR
      ibeg = data->beg[IDIR];
      iend = data->end[IDIR];
      jbeg = data->beg[KDIR];
      jend = data->end[KDIR];
      xi = data->x[IDIR];
      xj = data->x[KDIR];
      dxk = data->dx[JDIR];
    #elif GEOMETRY == SPHERICAL
      ibeg = data->beg[IDIR];
      iend = data->end[IDIR];
      jbeg = data->beg[JDIR];
      jend = data->end[JDIR];
      xi = data->x[IDIR];
      xj = data->x[JDIR];
      dxk = data->dx[KDIR];
    #endif


    idefix_reduce("CheckMaxDt", jbeg, jend, ibeg, iend,
      KOKKOS_LAMBDA(int j, int i, real &invDtLoc) {
        real w,dphi;
        #if GEOMETRY == CARTESIAN
          if(fargoType==userdef) {
            w = meanV(j,i);
          } else if(fargoType==shearingbox) {
            w = sbS*xi(i);
          }
        #elif GEOMETRY == POLAR
          w = meanV(j,i)/xi(i);
        #elif GEOMETRY == SPHERICAL
          w = meanV(j,i)/(xi(i)*sin(xj(j)));
        #endif
        dphi = dxk(0);  // dphi is supposedly constant when using fargo.
        invDtLoc = FMAX(invDtLoc, FABS(w/dphi));
      },
      Kokkos::Max<real>(invDt));
  #ifdef WITH_MPI
    if(idfx::psize>1) {
          MPI_SAFE_CALL(MPI_Allreduce(MPI_IN_PLACE, &invDt, 1, realMPI, MPI_MAX, MPI_COMM_WORLD));
        }
  #endif
  this->dtMax = this->maxShift / invDt;
}

void Fargo::AddVelocity(const real t) {
  idfx::pushRegion("Fargo::AddVelocity");

  this->AddVelocityFluid(t, data->hydro.get());
  if(data->haveDust) {
    for(int i = 0 ; i < data->dust.size() ; i++) {
      this->AddVelocityFluid(t, data->dust[i].get());
    }
  }

  idfx::popRegion();
}

void Fargo::SubstractVelocity(const real t) {
  idfx::pushRegion("Fargo::SubstractVelocity");

  this->SubstractVelocityFluid(t, data->hydro.get());
  if(data->haveDust) {
    for(int i = 0 ; i < data->dust.size() ; i++) {
      this->SubstractVelocityFluid(t, data->dust[i].get());
    }
  }

  idfx::popRegion();
}

void Fargo::ShiftSolution(const real t, const real dt) {
  idfx::pushRegion("Fargo::ShiftFluid");

  this->ShiftFluid(t,dt,data->hydro.get());
  if(data->haveDust) {
    for(int i = 0 ; i < data->dust.size() ; i++) {
      this->ShiftFluid(t,dt,data->dust[i].get());
    }
  }

  idfx::popRegion();
}
