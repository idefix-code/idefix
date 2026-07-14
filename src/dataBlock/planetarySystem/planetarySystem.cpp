// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include <string>
#include <vector>
#include "planetarySystem.hpp"
#include "planet.hpp"
#include "planetStructs.hpp"
#include "dataBlock.hpp"
#include "fluid.hpp"
#include "gravity.hpp"


PlanetarySystem::PlanetarySystem(Input &input, DataBlock *datain) {
  idfx::pushRegion("PlanetarySystem::Init");
  this->data = datain;
  this->feelDisk = input.Get<bool>("Planet","feelDisk",0);
  this->feelPlanets = input.Get<bool>("Planet","feelPlanets",0);
  // read Solver from input file
  std::string pintegratorString = input.GetOrSet<std::string>("Planet","integrator",0,"rk4");
  if (pintegratorString.compare("analytical") == 0) {
    this->myPlanetaryIntegrator = Integrator::ANALYTICAL;
    if ((this->feelDisk) || (this->feelPlanets)) {
      IDEFIX_ERROR("No feelDisk nor feelPlanets if analytical orbit.\n\
        Change integrator. You can also change feelDisk/feelPlanets.");
    }
  } else if (pintegratorString.compare("rk4") == 0) {
    this->myPlanetaryIntegrator = Integrator::RK4;
  } else if (pintegratorString.compare("rk5") == 0) {
    this->myPlanetaryIntegrator = Integrator::RK5;
  } else {
    std::stringstream msg;
    msg << "Unknown planet integrator type " << pintegratorString;
    IDEFIX_ERROR(msg);
  }

  std::string sfunctionString = input.Get<std::string>("Planet","smoothing",0);
  if (sfunctionString.compare("plummer") == 0) {
    this->myPlanetarySmoothing = PLUMMER;
  } else if (sfunctionString.compare("polynomial") == 0) {
    this->myPlanetarySmoothing = POLYNOMIAL;
  } else {
    std::stringstream msg;
    msg << "Unknown smoothing function type " << sfunctionString;
    IDEFIX_ERROR(msg);
  }

  this->torqueNormalization = input.GetOrSet<real>("Planet","torqueNormalization",0, 1.0);
  this->massTaper = input.GetOrSet<real>("Planet","masstaper",0, ZERO_F);
  this->excludeHill = input.GetOrSet<bool>("Planet","hillCut",0, false);
  this->indirectPlanetsTerm = input.GetOrSet<bool>("Planet","indirectPlanets",0, true);
  this->smoothingValue = input.Get<real>("Planet","smoothing",1);
  this->smoothingExponent = input.Get<real>("Planet","smoothing",2);

  // Initialize the planet object attached to this datablock
  this->nbp = input.CheckEntry("Planet","planetToPrimary");
  if (this->nbp<0) this->nbp=0;

  if ((this->nbp>1)
    && ((this->myPlanetaryIntegrator == Integrator::RK4)
      || (this->myPlanetaryIntegrator == Integrator::RK5))
    && (!(this->feelPlanets))
    && (this->indirectPlanetsTerm)) {
    IDEFIX_WARNING("Careful, the results are unphysical if Runge-Kutta with\n\
indirect term, with multiple planets that don't feel each others.");
  }

  if(this->nbp>0) {
    for(int ip = 0 ; ip < this->nbp ; ip++) {
      this->planet.emplace_back(ip, input, this->data, this);
    }
    // Now that we have initialised our planets, we register the variables for dump output
    // We cannot do this in the first loop since emplace_back might do some hidden copies
    // which then messes up the pointers used by dump i/O
    for(int ip = 0 ; ip < this->nbp ; ip++) {
      this->planet[ip].RegisterInDump();
    }
  } else {
    IDEFIX_ERROR("need to define a planet-to-primary mass ratio via planetToPrimary");
  }
#if GEOMETRY == POLAR || GEOMETRY == CARTESIAN
  #if DIMENSIONS == 3
    if ((this->data->mygrid->xbeg[KDIR] == 0)
        || (this->data->mygrid->xend[KDIR] == 0)) {
      this->halfdisk = true;
    } else {
      this->halfdisk = false;
    }
  #else // DIMENSIONS == 3
    this->halfdisk = false;
  #endif // DIMENSIONS == 3
#elif GEOMETRY == SPHERICAL
  if ( (this->data->mygrid->xbeg[JDIR] == M_PI/2.0)
        || (this->data->mygrid->xend[JDIR] == M_PI/2.0)) {
    this->halfdisk = true;
  } else {
    this->halfdisk = false;
  }
#endif

  idfx::popRegion();
}

void PlanetarySystem::ShowConfig() {
  idfx::pushRegion("PlanetarySystem::ShowConfig");
  idfx::cout << "PlanetarySystem: have " << this->nbp << " planets." << std::endl;
  switch(this->myPlanetaryIntegrator) {
    case ANALYTICAL:
      idfx::cout << "PlanetarySystem: uses analytical integration for planet location."
                 << std::endl;
      break;
    case RK4:
      idfx::cout << "PlanetarySystem: uses RK4 integration for planet location." << std::endl;
      break;
    case RK5:
      idfx::cout << "PlanetarySystem: uses RK5 integration for planet location." << std::endl;
      break;
    default:
      IDEFIX_ERROR("Unknown time integrator for planets");
  idfx::cout << "PlanetarySystem: feelDisk: " << this->feelDisk <<std::endl;
  idfx::cout << "PlanetarySystem: feelPlanets: " << this->feelPlanets <<std::endl;
  }

  switch(this->myPlanetarySmoothing) {
    case PLUMMER:
      idfx::cout << "PlanetarySystem: uses plummer expression for planet potential."
                 << std::endl;
      break;
    case POLYNOMIAL:
      idfx::cout << "PlanetarySystem: uses polynomial expression for planet potential."
                 << std::endl;
      break;
    default:
      IDEFIX_ERROR("Unknown smoothing function for planet potential");
  }

  if (this->halfdisk) {
    idfx::cout << "PlanetarySystem: half disk is detected, planet torques are computed accordingly."
               << std::endl;
  }

  // walk and show the planets
  for(Planet p : this->planet) {
    p.ShowConfig();
  }

  idfx::popRegion();
}

real PlanetarySystem::GetSmoothingValue() const {
  return this->smoothingValue;
}

real PlanetarySystem::GetSmoothingExponent() const {
  return this->smoothingExponent;
}

void PlanetarySystem::EvolveSystem(DataBlock& data, const real& dt) {
  idfx::pushRegion("PlanetarySystem::EvolveSystem");
  if(this->feelDisk) {
    this->AdvancePlanetFromDisk(data, data.dt);
  }
  this->IntegratePlanets(data, data.dt);
  idfx::popRegion();
}

void PlanetarySystem::AdvancePlanetFromDisk(DataBlock& data, const real& dt) {
  idfx::pushRegion("PlanetarySystem::AdvancePlanetFromDisk");
  for(int ip=0; ip< this->nbp ; ip++) {
    if (!(planet[ip].m_isActive)) continue;
    Point gamma;
    bool isp = true;

    real qp = planet[ip].m_qp;
    real xp = planet[ip].m_xp;
    real yp = planet[ip].m_yp;
    real zp = planet[ip].m_zp;
    real vxp = planet[ip].m_vxp;
    real vyp = planet[ip].m_vyp;
    real vzp = planet[ip].m_vzp;
    real r = sqrt(xp*xp + yp*yp + zp*zp);

    gamma = planet[ip].computeAccel (data, isp);

    planet[ip].m_vxp += dt * gamma.x*this->torqueNormalization;
    planet[ip].m_vyp += dt * gamma.y*this->torqueNormalization;
    planet[ip].m_vzp += dt * gamma.z*this->torqueNormalization;
  }
  idfx::popRegion();
}

void PlanetarySystem::IntegratePlanets(DataBlock& data, const real& dt) {
    switch(this->myPlanetaryIntegrator) {
        case ANALYTICAL:
          IntegrateAnalytically(data, dt);
          break;
        case RK4:
          {
            int subcycling = 1;
            int i;
            for (i = 0; i < subcycling; i++)
              IntegrateRK4(data, 1.0/(static_cast<double>(subcycling))*dt);
            break;
          }
        case RK5:
          {
            int subcycling = 1;
            int i;
            for (i = 0; i < subcycling; i++)
              IntegrateRK5(data, 1.0/(static_cast<double>(subcycling))*dt);
            break;
          }
        default: // do nothing
          break;
    }
}

void PlanetarySystem::IntegrateAnalytically(DataBlock& data, const real& dt) {
  idfx::pushRegion("PlanetarySystem::IntegrateAnalytically");
  for(int ip=0; ip<this->nbp ; ip++) {
    if (!(planet[ip].m_isActive)) continue;

    real qp = planet[ip].m_qp;
    real xp = planet[ip].m_xp;
    real yp = planet[ip].m_yp;

    // from cartesian to polar
    real Rp{std::sqrt(xp*xp+yp*yp)};
    real phip{std::atan2(yp,xp)};

    real omegap{std::sqrt((ONE_F+qp)/Rp/Rp/Rp)};

    real RpNew{Rp};
//    real phipNew{phip + dt*omegap};
    real phipNew;

    if(planet[ip].data->hydro->haveRotation) {
      real omegaframe = planet[ip].data->hydro->OmegaZ;
      phipNew = phip + dt*(omegap-omegaframe);
    } else {
      phipNew = phip + dt*omegap;
    }

    real vxp_tmp{ZERO_F};
    real vyp_tmp{RpNew*std::sqrt((ONE_F+qp)/RpNew/RpNew/RpNew)};

    planet[ip].m_xp = RpNew*cos(phipNew);
    planet[ip].m_yp = RpNew*sin(phipNew);
    planet[ip].m_zp = ZERO_F;
    planet[ip].m_vxp = vxp_tmp*cos(phipNew)-vyp_tmp*sin(phipNew);
    planet[ip].m_vyp = vxp_tmp*sin(phipNew)+vyp_tmp*cos(phipNew);
    planet[ip].m_vzp = ZERO_F;
  }
  idfx::popRegion();
}

void PlanetarySystem::IntegrateRK5(DataBlock& data, const real& dt) {
  idfx::pushRegion("PlanetarySystem::IntegrateRK5");
  std::vector<Planet> &ki = planet;

  real t1, t2, t3, t4;
  std::vector<Planet> k0 = ki;
  std::vector<Planet> k1 = ki;
  std::vector<Planet> k2 = ki;
  std::vector<Planet> k3 = ki;
  std::vector<Planet> k4 = ki;
  std::vector<Planet> kf = ki;

  // No need to copy since we already copied when instantiating k0

  std::vector<PointSpeed> f0 = ComputeRHS(data.t, k0);
  t1 = data.t + dt/4.0;
  for(int ip=0; ip<this->nbp; ip++) {
      k1[ip].state = ki[ip].state + dt*f0[ip]/4;
  }

  std::vector<PointSpeed> f1 = ComputeRHS(t1, k1);
  t2 = data.t + 3*dt/8.0;
  for(int ip=0; ip<this->nbp; ip++) {
      k2[ip].state = ki[ip].state + 3*dt*f0[ip]/32 + 9*dt*f1[ip]/32;
  }

  std::vector<PointSpeed> f2 = ComputeRHS(t2, k2);
  t3 = data.t + 12*dt/13.0;
  for(int ip=0; ip<this->nbp; ip++) {
      k3[ip].state = ki[ip].state + 1932*dt*f0[ip]/2197 - 7200*dt*f1[ip]/2197 + 7296*dt*f2[ip]/2197;
  }

  std::vector<PointSpeed> f3 = ComputeRHS(t3, k3);
  t4 = data.t + dt;
  for(int ip=0; ip<this->nbp; ip++) {
      k4[ip].state = ki[ip].state + 439*dt*f0[ip]/216 - 8*dt*f1[ip] +
      3680*dt*f2[ip]/513 - 845*dt*f3[ip]/4104;
  }

  std::vector<PointSpeed> f4 = ComputeRHS(t4, k4);
  for(int ip=0; ip<this->nbp; ip++) {
      kf[ip].state = ki[ip].state + dt * ( 25*f0[ip]/216 + 1408*f2[ip]/2565 +
      2197*f3[ip]/4104 - f4[ip]/5);
  }

  planet = kf;
  idfx::popRegion();
}

void PlanetarySystem::IntegrateRK4(DataBlock& data, const real& dt) {
  idfx::pushRegion("PlanetarySystem::IntegrateRK4");
  std::vector<Planet> &ki = planet;

  real t1, t2, t3;
  std::vector<Planet> k0 = ki;
  std::vector<Planet> k1 = ki;
  std::vector<Planet> k2 = ki;
  std::vector<Planet> k3 = ki;
  std::vector<Planet> kf = ki;

  // No need to copy since we already copied when instantiating k0

  std::vector<PointSpeed> f0 = ComputeRHS(data.t, k0);
  t1 = data.t + dt/2.0;
  for(int ip=0; ip<this->nbp; ip++) {
      k1[ip].state = ki[ip].state + dt*f0[ip]/2;
  }

  std::vector<PointSpeed> f1 = ComputeRHS(t1, k1);
  t2 = data.t + dt/2.0;
  for(int ip=0; ip<this->nbp; ip++) {
      k2[ip].state = ki[ip].state + dt*f1[ip]/2;
  }

  std::vector<PointSpeed> f2 = ComputeRHS(t2, k2);
  t3 = data.t + dt;
  for(int ip=0; ip<this->nbp; ip++) {
    k3[ip].state = ki[ip].state + dt*f2[ip];
  }

  std::vector<PointSpeed> f3 = ComputeRHS(t3, k3);
  for(int ip=0; ip<this->nbp; ip++) {
      kf[ip].state = ki[ip].state + dt * ( f0[ip] + 2*f1[ip] + 2*f2[ip] + f3[ip]) / 6.0;
  }

  planet = kf;
  idfx::popRegion();
}

std::vector<PointSpeed> PlanetarySystem::ComputeRHS(real& t, std::vector<Planet> planet) {
  std::vector<PointSpeed> planet_update(this->nbp);

  for(int ip=0; ip<this->nbp; ip++) {
    planet_update[ip].x = ZERO_F;
    planet_update[ip].y = ZERO_F;
    planet_update[ip].z = ZERO_F;
    planet_update[ip].vx = ZERO_F;
    planet_update[ip].vy = ZERO_F;
    planet_update[ip].vz = ZERO_F;

    if (!(planet[ip].m_isActive)) continue;

    real dist;
    real coef;
    dist = sqrt(
      planet[ip].m_xp*planet[ip].m_xp +
      planet[ip].m_yp*planet[ip].m_yp +
      planet[ip].m_zp*planet[ip].m_zp
    );
    planet_update[ip].x = planet[ip].m_vxp;
    planet_update[ip].y = planet[ip].m_vyp;
    planet_update[ip].z = planet[ip].m_vzp;
    planet_update[ip].vx = -planet[ip].m_xp/dist/dist/dist;
    planet_update[ip].vy = -planet[ip].m_yp/dist/dist/dist;
    planet_update[ip].vz = -planet[ip].m_zp/dist/dist/dist;

    for(int jp=0; jp<this->nbp; jp++) {
        if(this->indirectPlanetsTerm && planet[jp].m_isActive) {
            dist = sqrt(
              planet[jp].m_xp*planet[jp].m_xp +
              planet[jp].m_yp*planet[jp].m_yp +
              planet[jp].m_zp*planet[jp].m_zp
            );
            coef = planet[jp].m_qp/dist/dist/dist;
            planet_update[ip].vx -= coef*planet[jp].m_xp;
            planet_update[ip].vy -= coef*planet[jp].m_yp;
            planet_update[ip].vz -= coef*planet[jp].m_zp;
        }

        if((jp != ip) && this->feelPlanets && planet[jp].m_isActive) {
            dist = sqrt(
              (planet[ip].m_xp-planet[jp].m_xp)*(planet[ip].m_xp-planet[jp].m_xp) +
              (planet[ip].m_yp-planet[jp].m_yp)*(planet[ip].m_yp-planet[jp].m_yp) +
              (planet[ip].m_zp-planet[jp].m_zp)*(planet[ip].m_zp-planet[jp].m_zp)
            );
            coef = planet[jp].m_qp/dist/dist/dist;
            planet_update[ip].vx += coef*(planet[jp].m_xp-planet[ip].m_xp);
            planet_update[ip].vy += coef*(planet[jp].m_yp-planet[ip].m_yp);
            planet_update[ip].vz += coef*(planet[jp].m_zp-planet[ip].m_zp);
        }
    }
    if(planet[ip].data->hydro->haveRotation) {
      real omega = planet[ip].data->hydro->OmegaZ;


      // Add Coriolis and centrifugal forces
      planet_update[ip].vx +=  2*omega* planet[ip].m_vyp + omega*omega*planet[ip].m_xp;
      planet_update[ip].vy += -2*omega* planet[ip].m_vxp + omega*omega*planet[ip].m_yp;
    }
  }
  return planet_update;
}

void PlanetarySystem::AddPlanetsPotential(IdefixArray3D<real> &phiP, real t) {
  idfx::pushRegion("PlanetarySystem::AddPlanetsPotential");
  bool indirectPlanetsTerm = this->indirectPlanetsTerm;
  real smoothingValue = this->smoothingValue;
  real smoothingExponent = this->smoothingExponent;
  SmoothingFunction myPlanetarySmoothing = this->myPlanetarySmoothing;

  IdefixArray1D<real> x1 = this->data->x[IDIR];
  IdefixArray1D<real> x2 = this->data->x[JDIR];
  IdefixArray1D<real> x3 = this->data->x[KDIR];

  for(Planet& p : this->planet) {
    // update mass according to mass taper
    p.updateMp(t);
    p.activatePlanet(t);

    bool isActive = p.getIsActive();
    if (!(isActive)) continue;

    real qp = p.getMp();
    real xp = p.getXp();
    real yp = p.getYp();
    real zp = p.getZp();

    real distPlanet = sqrt(xp*xp+yp*yp+zp*zp);
    real smoothing = smoothingValue * pow(distPlanet,1.0+smoothingExponent);
    real Mcentral = this->data->gravity->centralMass;

    idefix_for("PlanetPotential",
      0,this->data->np_tot[KDIR],
      0, this->data->np_tot[JDIR],
      0, this->data->np_tot[IDIR],
        KOKKOS_LAMBDA (int k, int j, int i) {
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

        real dist = ((xc-xp)*(xc-xp)+
                    (yc-yp)*(yc-yp)+
                    (zc-zp)*(zc-zp));

        // term due to planet
        switch(myPlanetarySmoothing) {
            case PLUMMER:
              {
                phiP(k,j,i) += -Mcentral*qp/sqrt(dist+smoothing*smoothing);
                break;
              }
            case POLYNOMIAL:
              {
                real rmrp = sqrt(dist);
                if (rmrp/smoothing < 1) {
                  phiP(k,j,i) += -(Mcentral*qp/rmrp)*(pow(rmrp/smoothing,4.0) -
                                                     2.0*pow(rmrp/smoothing,3.0)+
                                                     2.0*rmrp/smoothing);
                } else {
                  phiP(k,j,i) += -(Mcentral*qp/rmrp);
                }
                break;
              }
            default: // do nothing
              break;
        }
        // indirect term due to planet
        if (indirectPlanetsTerm) {
          phiP(k,j,i) += Mcentral*qp*(xc*xp+yc*yp+zc*zp)/(distPlanet*distPlanet*distPlanet);
        }
    });
  }

  idfx::popRegion();
}
