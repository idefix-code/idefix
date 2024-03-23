// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include <iostream>
#include <string>
#include "planet.hpp"
#include "dataBlock.hpp"
#include "planetarySystem.hpp"
#include "fluid.hpp"
/*
Planet::Planet() :
    m_vxp(state.vx),
    m_vyp(state.vy),
    m_vzp(state.vz),
    m_xp(state.x),
    m_yp(state.y),
    m_zp(state.z) {
}
*/

Planet::Planet(const Planet& p) :
    m_vxp(state.vx),
    m_vyp(state.vy),
    m_vzp(state.vz),
    m_xp(state.x),
    m_yp(state.y),
    m_zp(state.z) {
    idfx::pushRegion("Planet::Planet(Planet)");
    this->operator=(p);
    idfx::popRegion();
}

// Assignement operator (required since we internally use references)
Planet& Planet::operator=(const Planet &p) {
  idfx::pushRegion("Planet::operator=");
  this->state = p.state;
  this->m_qp = p.m_qp;
  this->m_qpIni = p.m_qpIni;
  this->m_force = p.m_force;
  this->m_tOffset = p.m_tOffset;
  this->m_isActive = p.m_isActive;
  this->m_ip = p.m_ip;
  this->pSys = p.pSys;
  this->data = p.data;
  idfx::popRegion();
  return *this;
}

Planet::Planet(int ip, Input &input, DataBlock *datain, PlanetarySystem *pSys):
m_vxp(state.vx), m_vyp(state.vy), m_vzp(state.vz), m_xp(state.x), m_yp(state.y), m_zp(state.z) {
  idfx::pushRegion("Planet::Planet");

  // Save the datablock to which we are attached from now on
  this->data = datain;
  this->pSys = pSys;

  this->m_ip = ip;

  real eccentricity{input.GetOrSet<real>("Planet","initialEccentricity",ip, ZERO_F)};
  real inclination{input.GetOrSet<real>("Planet","initialInclination",ip, ZERO_F)};
  real distIni{input.Get<real>("Planet","initialDistance",ip)};
  real massTaper = pSys->massTaper;

  std::string pintegratorString = input.Get<std::string>("Planet","integrator",0);
  if (pintegratorString.compare("analytical") == 0) {
    if (eccentricity != ZERO_F) {
      IDEFIX_ERROR(
        "planet "+std::to_string(ip)
        +": no planet eccentricity yet\n\
        if analytical orbit. Try initialEccentricity = 0.0.\n\
          You can also change integrator.");
    }
    if (inclination != ZERO_F) {
      IDEFIX_ERROR("planet "+std::to_string(ip)
        +": no planet inclination yet\n\
        if analytical orbit. Try initialInclination = 0.0.\n\
          You can also change integrator.");
    }
  }

  this->m_qpIni = input.Get<real>("Planet","planetToPrimary",ip);
  this->m_tOffset = input.GetOrSet<real>("Planet","tOffset",ip, ZERO_F);

  if (this->m_tOffset == ZERO_F) {
    this->m_isActive = true;
  } else {
    this->m_isActive = false;
  }

  // @GL: We need to discuss with GI Jonah in order to be able to add the gas indirect term if SG

  this->m_force = {
    {ZERO_F,ZERO_F,ZERO_F},
    {ZERO_F,ZERO_F,ZERO_F},
    {ZERO_F,ZERO_F,ZERO_F},
    {ZERO_F,ZERO_F,ZERO_F}
  };

  if (massTaper == ZERO_F) {
    this->m_qp = this->m_qpIni;
  } else {
    this->m_qp = ZERO_F;
  }
  // We initialize the quantities required by the planet solver
  this->m_xp = distIni*(ONE_F+eccentricity);
  this->m_yp = ZERO_F;
  this->m_zp = ZERO_F;
  this->m_vxp = ZERO_F;
  if (distIni>=ZERO_F) {
    this->m_vyp = sqrt((ONE_F+this->m_qp)/FABS(distIni))*
                  sqrt((ONE_F-eccentricity)/(ONE_F+eccentricity))*cos(inclination);
  } else {
    this->m_vyp = -sqrt((ONE_F+this->m_qp)/FABS(distIni))*
                  sqrt((ONE_F-eccentricity)/(ONE_F+eccentricity))*cos(inclination);
  }
  // Remove inertial rotating frame velocity if present
  if(this->data->hydro->haveRotation) {
    this->m_vyp += -this->data->hydro->OmegaZ * this->m_xp;
  }
  this->m_vzp = this->m_vyp*sin(inclination)/cos(inclination);

  idfx::popRegion();
}

void Planet::RegisterInDump() {
  idfx::pushRegion("Planet::RegisterInDump");
  // Register variables for dump read/write
  data->dump->RegisterVariable(&m_xp,std::string("x_p")+std::to_string(m_ip));
  data->dump->RegisterVariable(&m_yp,std::string("y_p")+std::to_string(m_ip));
  data->dump->RegisterVariable(&m_zp,std::string("z_p")+std::to_string(m_ip));
  data->dump->RegisterVariable(&m_vxp,std::string("vx_p")+std::to_string(m_ip));
  data->dump->RegisterVariable(&m_vyp,std::string("vy_p")+std::to_string(m_ip));
  data->dump->RegisterVariable(&m_vzp,std::string("vz_p")+std::to_string(m_ip));
  data->dump->RegisterVariable(&m_qp,std::string("q_p")+std::to_string(m_ip));

  idfx::popRegion();
}

void Planet::ShowConfig() {
  idfx::cout << "Planet[" << this->m_ip << "]: mass qp=" << this->m_qpIni <<std::endl;
  idfx::cout << "Planet[" << this->m_ip << "]: initial location dp="
    << sqrt(pow(this->m_xp,2)+pow(this->m_yp,2)+pow(this->m_zp,2)) <<std::endl;
  // @ GL What about the initial eccentricity/inclination?
  // --> not arguments of the class Planet...
}

void Planet::displayPlanet() const {
  idfx::cout << "Planet: Init orbital parameters. " << std::endl;
  idfx::cout << "xp = " << this->m_xp << std::endl;
  idfx::cout << "yp = " << this->m_yp << std::endl;
  idfx::cout << "zp = " << this->m_zp << std::endl;
  idfx::cout << "vxp = " << this->m_vxp << std::endl;
  idfx::cout << "vyp = " << this->m_vyp << std::endl;
  idfx::cout << "vzp = " << this->m_vzp << std::endl;
  idfx::cout << "qp = " << this->m_qp << std::endl;
}

real Planet::getMp() const {
  return this->m_qp;
}

real Planet::getXp() const {
  return this->m_xp;
}

real Planet::getYp() const {
  return this->m_yp;
}

real Planet::getZp() const {
  return this->m_zp;
}

real Planet::getVxp() const {
  return this->m_vxp;
}

real Planet::getVyp() const {
  return this->m_vyp;
}

real Planet::getVzp() const {
  return this->m_vzp;
}

bool Planet::getIsActive() const {
  return this->m_isActive;
}

int Planet::getIndex() const {
  return(this->m_ip);
}

void Planet::setMp(real qp) {
  this->m_qp = qp;
}

void Planet::setXp(real xp) {
  this->m_xp = xp;
}

void Planet::setYp(real yp) {
  this->m_yp = yp;
}

void Planet::setZp(real zp) {
  this->m_zp = zp;
}

void Planet::setVxp(real vxp) {
  this->m_vxp = vxp;
}

void Planet::setVyp(real vyp) {
  this->m_vyp = vyp;
}

void Planet::setVzp(real vzp) {
  this->m_vzp = vzp;
}


void Planet::activatePlanet(const real t) {
  if (t >= this->m_tOffset) {
    this->m_isActive = true;
  } else {
    this->m_isActive = false;
  }
}


void Planet::updateMp(const real t) {
  // idfx::cout << "tOffset = " << this->m_tOffset << std::endl;
  real mtaper = pSys->massTaper;
  real factor_mtaper = ONE_F;
  if (mtaper > ZERO_F) {
    if (t >= this->m_tOffset) {
      // changed t to t-tOffset
      factor_mtaper = (
        (t-this->m_tOffset) >= mtaper ? ONE_F : .5*(ONE_F-cos(M_PI*(t-this->m_tOffset)/mtaper))
      );
    } else {
      factor_mtaper = ZERO_F;
    }
  }
  this->m_qp = this->m_qpIni*factor_mtaper;
}

Point Planet::computeAccel(DataBlock& data, bool& isPlanet) {
  Point acceleration;
  Force &force = this->m_force;
//  Force &force = data.planet[0].force;
  computeForce(data,isPlanet);
  bool excludeHill = pSys->excludeHill;
  if (excludeHill) {
    acceleration.x = force.f_ex_inner[0]+force.f_ex_outer[0];
    acceleration.y = force.f_ex_inner[1]+force.f_ex_outer[1];
    acceleration.z = force.f_ex_inner[2]+force.f_ex_outer[2];
  } else {
    acceleration.x = force.f_inner[0]+force.f_outer[0];
    acceleration.y = force.f_inner[1]+force.f_outer[1];
    acceleration.z = force.f_inner[2]+force.f_outer[2];
  }
  return acceleration;
}

/*
Be careful: you need to substract
the azimuthally averaged density
prior to the torque evaluation (BM08 trick)
*/
void Planet::computeForce(DataBlock& data, bool& isPlanet) {
  PlanetarySystem::SmoothingFunction smoothingFunction = pSys->myPlanetarySmoothing;
  real smoothingValue = pSys->smoothingValue;
  real smoothingExponent = pSys->smoothingExponent;

  IdefixArray1D<real> x1 = data.x[IDIR];
  IdefixArray1D<real> x2 = data.x[JDIR];
  IdefixArray1D<real> x3 = data.x[KDIR];

  IdefixArray4D<real> Vc = data.hydro->Vc;
  IdefixArray3D<real> dV = data.dV;

  real xp;
  real yp;
  real zp;
  real qp;

  if(isPlanet) {
    xp = this->m_xp;
    yp = this->m_yp;
    zp = this->m_zp;
    qp = this->m_qp;
  } else {
    xp = ZERO_F;
    yp = ZERO_F;
    zp = ZERO_F;
    qp = ZERO_F;
  }
  real distPlanet = sqrt(xp*xp+yp*yp+zp*zp);
  real smoothing = smoothingValue * pow(distPlanet,ONE_F+smoothingExponent);
  real rh = pow(qp/3., 1./3.)*distPlanet;
  bool excludeHill = pSys->excludeHill;

  // since we cannot throw an error in kokkos kernel, with throw this one before the kernel.
  #if GEOMETRY == CYLINDRICAL
    IDEFIX_ERROR("Planet::ComputeForce is not compatible with the GEOMETRY you intend to use");
  #endif

  Kokkos::parallel_reduce("ComputeForce",
    Kokkos::MDRangePolicy<Kokkos::Rank<3, Kokkos::Iterate::Right, Kokkos::Iterate::Right>>
    ({data.beg[KDIR],data.beg[JDIR],data.beg[IDIR]},
      {data.end[KDIR], data.end[JDIR], data.end[IDIR]}),
    KOKKOS_LAMBDA (int k, int j, int i, Force &forceProc) {
      real cellMass = dV(k,j,i)*Vc(RHO,k,j,i);
      real xc, yc, zc;
      #if GEOMETRY == CARTESIAN
        xc = x1(i);
        yc = x2(j);
        zc = x3(k);
      #elif GEOMETRY == POLAR
        xc = x1(i)*cos(x2(j));
        yc = x1(i)*sin(x2(j));
        zc = x3(k);
      #elif GEOMETRY == SPHERICAL
        xc = x1(i)*sin(x2(j))*cos(x3(k));
        yc = x1(i)*sin(x2(j))*sin(x3(k));
        zc = x1(i)*cos(x2(j));
      #endif

      real distc = sqrt(xc*xc+yc*yc+zc*zc);
      real dist2 = ((xc-xp)*(xc-xp) + (yc-yp)*(yc-yp) + (zc-zp)*(zc-zp));
      real hillcut;

      if(excludeHill) {
        real squaredist2 = sqrt(dist2);
        if (squaredist2/rh < 0.5) {
          hillcut = ZERO_F;
        } else {
          if (squaredist2 > rh) {
            hillcut = ONE_F;
          } else {
            hillcut = pow(sin((squaredist2/rh-.5)*M_PI),2.);
          }
        }
      }

      real forceCell;
      switch(smoothingFunction) {
        case PlanetarySystem::SmoothingFunction::PLUMMER:
          {
            dist2 += smoothing*smoothing; // if default potential
            real distance = sqrt(dist2); // if default potential
            real InvDist3 = ONE_F/(dist2*distance); // if default potential
            forceCell = cellMass * InvDist3; // if default potential
            break;
          }
        case PlanetarySystem::SmoothingFunction::POLYNOMIAL:
          {
            real rmrp = sqrt(dist2); // if other potential
            if (rmrp/smoothing < 1) {
              forceCell = -cellMass*(3.0*rmrp/smoothing - 4.0)/smoothing/smoothing/smoothing;
            } else {
              forceCell = cellMass/rmrp/rmrp/rmrp;
            }
            break;
          }
        default: // do nothing
          break;
      }
      if (distc < distPlanet) {
        // INNER FORCE
        forceProc.f_inner[0] += (xc-xp)*forceCell;
        forceProc.f_inner[1] += (yc-yp)*forceCell;
        forceProc.f_inner[2] += (zc-zp)*forceCell;
        if(excludeHill) {
          forceProc.f_ex_inner[0] += (xc-xp)*forceCell*hillcut;
          forceProc.f_ex_inner[1] += (yc-yp)*forceCell*hillcut;
          forceProc.f_ex_inner[2] += (zc-zp)*forceCell*hillcut;
        } else {
          forceProc.f_ex_inner[0] += ZERO_F;
          forceProc.f_ex_inner[1] += ZERO_F;
          forceProc.f_ex_inner[2] += ZERO_F;
        }
      } else {
        // OUTER FORCE
        forceProc.f_outer[0] += (xc-xp)*forceCell;
        forceProc.f_outer[1] += (yc-yp)*forceCell;
        forceProc.f_outer[2] += (zc-zp)*forceCell;
        if(excludeHill) {
          forceProc.f_ex_outer[0] += (xc-xp)*forceCell*hillcut;
          forceProc.f_ex_outer[1] += (yc-yp)*forceCell*hillcut;
          forceProc.f_ex_outer[2] += (zc-zp)*forceCell*hillcut;
        } else {
          forceProc.f_ex_outer[0] += ZERO_F;
          forceProc.f_ex_outer[1] += ZERO_F;
          forceProc.f_ex_outer[2] += ZERO_F;
        }
      }
    }, this->m_force );

  if(pSys->halfdisk) {
    // Cancel vertical component
    m_force.f_inner[2] = 0;
    m_force.f_ex_inner[2] = 0;
    m_force.f_outer[2] = 0;
    m_force.f_ex_outer[2] = 0;

    // Multiply by 2 the remaining components
    for(int i = 0 ; i < 2 ; i++) {
      m_force.f_inner[i] *= 2;
      m_force.f_ex_inner[i] *= 2;
      m_force.f_outer[i] *= 2;
      m_force.f_ex_outer[i] *= 2;
    }
  }

  #ifdef WITH_MPI
    real forceGlob[12];
    real forceLoc[12];
    forceLoc[0] = this->m_force.f_inner[0];
    forceLoc[1] = this->m_force.f_inner[1];
    forceLoc[2] = this->m_force.f_inner[2];
    forceLoc[3] = this->m_force.f_ex_inner[0];
    forceLoc[4] = this->m_force.f_ex_inner[1];
    forceLoc[5] = this->m_force.f_ex_inner[2];
    forceLoc[6] = this->m_force.f_outer[0];
    forceLoc[7] = this->m_force.f_outer[1];
    forceLoc[8] = this->m_force.f_outer[2];
    forceLoc[9] = this->m_force.f_ex_outer[0];
    forceLoc[10] = this->m_force.f_ex_outer[1];
    forceLoc[11] = this->m_force.f_ex_outer[2];
    MPI_SAFE_CALL(MPI_Allreduce(&forceLoc, &forceGlob, 12, realMPI, MPI_SUM, MPI_COMM_WORLD));
    this->m_force.f_inner[0] = forceGlob[0];
    this->m_force.f_inner[1] = forceGlob[1];
    this->m_force.f_inner[2] = forceGlob[2];
    this->m_force.f_ex_inner[0] = forceGlob[3];
    this->m_force.f_ex_inner[1] = forceGlob[4];
    this->m_force.f_ex_inner[2] = forceGlob[5];
    this->m_force.f_outer[0] = forceGlob[6];
    this->m_force.f_outer[1] = forceGlob[7];
    this->m_force.f_outer[2] = forceGlob[8];
    this->m_force.f_ex_outer[0] = forceGlob[9];
    this->m_force.f_ex_outer[1] = forceGlob[10];
    this->m_force.f_ex_outer[2] = forceGlob[11];
  #endif
}
