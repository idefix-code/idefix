#include <algorithm>
#include <math.h>
#include "idefix.hpp"
#include "setup.hpp"
#include "planet.hpp"

real wkzMinGlob;
real wkzMaxGlob;
real wkzDampingGlob;
real sigma0Glob;
real sigmaSlopeGlob;
real h0Glob;
real flaringIndexGlob;
real alphaGlob;
real densityFloorGlob;
real masstaperGlob;
real omegaGlob;


void MySoundSpeed(DataBlock &data, const real t, IdefixArray3D<real> &cs) {
  real h0 = h0Glob;
  real flaringIndex = flaringIndexGlob;
  IdefixArray1D<real> x1=data.x[IDIR];
  IdefixArray1D<real> x2=data.x[JDIR];
  idefix_for("MySoundSpeed",0,data.np_tot[KDIR],0,data.np_tot[JDIR],0,data.np_tot[IDIR],
              KOKKOS_LAMBDA (int k, int j, int i) {
                real R = x1(i)*sin(x2(j));
                cs(k,j,i) = h0*pow(R,flaringIndex-0.5);
              });
}

void Damping(Hydro *hydro, const real t, const real dtin) {
  auto *data = hydro->data;
  IdefixArray4D<real> Vc = hydro->Vc;
  IdefixArray4D<real> Uc = hydro->Uc;
  IdefixArray1D<real> x1 = data->x[IDIR];
  IdefixArray1D<real> x2 = data->x[JDIR];

  real h0 = h0Glob;
  real sigma0 = sigma0Glob;
  real sigmaSlope = sigmaSlopeGlob;
  real flaringIndex = flaringIndexGlob;
  real omega = omegaGlob;
  real gamma = hydro->eos->GetGamma();
  real pexp = -sigmaSlope-flaringIndex-1;
  real qexp = 2*flaringIndex-1;
  real dt = dtin;
  bool isFargo = data->haveFargo;

  real rmin{data->mygrid->xbeg[0]};
  real wkzMin{wkzMinGlob};
  real rmax{data->mygrid->xend[0]};
  real wkzMax{wkzMaxGlob};
  real wkzDamping{wkzDampingGlob};
  idefix_for("Damping",
    0, data->np_tot[KDIR],
    0, data->np_tot[JDIR],
    0, data->np_tot[IDIR],
              KOKKOS_LAMBDA (int k, int j, int i) {
                real r = x1(i);
                real th = x2(j);
                real R = r*sin(th);
                real cs2 = h0*h0*pow(R,2*flaringIndex-1);
                real Rmin = rmin*sin(th);
                real Rmax = rmax*sin(th);
                real wkzRmin = wkzMin*sin(th);
                real wkzRmax = wkzMax*sin(th);
                real Vk = 1.0/sqrt(R);

                /*
                // Cooling function
                real cs2 = h0*h0*pow(R,2*flaringIndex-1.0);

                // cooling time in local Keplerian units
                real tau = 0.1*pow(R,1.5);
                real Ptarget = cs2*Vc(RHO,k,j,i);

                Uc(ENG,k,j,i) += -dt*(Vc(PRS,k,j,i)-Ptarget)/(tau*(gamma-1.0));
                */

                real lambda = 0.0;

                // Damping function for poloidal velocity field in damping zones
                // Damp whatever is at R<wkzRmin and R>wkzRmax
                if (R<wkzRmin) {
                  lambda = 1.0/(wkzDamping*2.0*M_PI*pow(Rmin,1.5))*(1.0 - pow(sin(M_PI*( (R-Rmin) / (wkzRmin-Rmin) )/2.0),2.0));
                }
                if (R>wkzRmax) {
                  lambda = 1.0/(wkzDamping*2.0*M_PI*pow(Rmax,1.5))*pow(sin(M_PI*( (R-wkzRmax) / (Rmax-wkzRmax) )/2.0),2.0);
                }

                real rhoTarget = sigma0/sqrt(2.0*M_PI)/(h0*pow(R,flaringIndex+1))*pow(R,-sigmaSlope) * exp(1.0/ (cs2) * (1.0/r-1.0/R)) ;
                real vx3Target = 0.0;
                if(!isFargo) {
                  vx3Target = Vk*sqrt((pexp+qexp)*h0*h0*pow(R,2*flaringIndex)+(1+qexp)-qexp*R/r) - omega * R;
                }

                // relaxation
                real drho = lambda*(Vc(RHO,k,j,i)-rhoTarget);
                real dvx1 = lambda*Vc(RHO,k,j,i)*Vc(VX1,k,j,i);
                real dvx2 = lambda*Vc(RHO,k,j,i)*Vc(VX2,k,j,i);
                real dvx3 = lambda*Vc(RHO,k,j,i)*(Vc(VX3,k,j,i)-vx3Target);

                real dmx1 = dvx1 + Vc(VX1,k,j,i) * drho;
                real dmx2 = dvx2 + Vc(VX2,k,j,i) * drho;
                real dmx3 = dvx3 + Vc(VX3,k,j,i) * drho;

                // Kinetic energy fluctuations due to above damping
                // must be compensated in total energy conservation
                /*
                real deng = 0.5*(Vc(VX1,k,j,i)*Vc(VX1,k,j,i)
                                +Vc(VX2,k,j,i)*Vc(VX2,k,j,i))*drho
                                +Vc(VX3,k,j,i)*Vc(VX3,k,j,i))*drho
                            + Vc(RHO,k,j,i) * (
                                  Vc(VX1,k,j,i)*dvx1
                                + Vc(VX2,k,j,i)*dvx2
                                + Vc(VX3,k,j,i)*dvx3
                            );
                */

                Uc(RHO,k,j,i) += -drho*dt;
                Uc(MX1,k,j,i) += -dmx1*dt;
                Uc(MX2,k,j,i) += -dmx2*dt;
                Uc(MX3,k,j,i) += -dmx3*dt;
                // Uc(ENG,k,j,i) += -deng*dt;

});

}

// User-defined boundaries
void UserdefBoundary(Hydro *hydro, int dir, BoundarySide side, real t) {
  auto *data = hydro->data;
  IdefixArray4D<real> Vc = hydro->Vc;
  IdefixArray1D<real> x1 = data->x[IDIR];
  IdefixArray1D<real> x2 = data->x[JDIR];
  real h0=h0Glob;
  real flaringIndex=flaringIndexGlob;
  real sigmaSlope=sigmaSlopeGlob;
  real omega = omegaGlob;
  real pexp = -sigmaSlope-flaringIndex-1;
  real qexp = 2*flaringIndex-1;

    if(dir==IDIR) {
        int ighost,ibeg,iend;
        if(side == left) {
            ighost = data->beg[IDIR];
            ibeg = 0;
            iend = data->beg[IDIR];
            //return;
        }
        else if(side==right) {
            ighost = data->end[IDIR]-1;
            ibeg=data->end[IDIR];
            iend=data->np_tot[IDIR];
        }


        idefix_for("UserDefBoundary",
          0, data->np_tot[KDIR],
          0, data->np_tot[JDIR],
          ibeg, iend,
                    KOKKOS_LAMBDA (int k, int j, int i) {
                        real r = x1(i);
                        real rghost = x1(ighost);
                        real R = x1(i)*sin(x2(j));
                        real Rghost = x1(ighost)*sin(x2(j));
                        real Vk = 1.0/sqrt(R);
                        real cs2 = h0*h0*pow(R,2*flaringIndex-1);

                        Vc(RHO,k,j,i) = Vc(RHO,k,j,ighost) * pow(R/Rghost,-sigmaSlope-flaringIndex-1) * exp(1.0/cs2 * (1.0/r - 1.0/rghost + 1.0/Rghost - 1.0/R)) ;
                        Vc(VX1,k,j,i) = ZERO_F;//- Vc(VX1,k,j,2*ighost-i+sign);
                        Vc(VX2,k,j,i) = ZERO_F;//Vc(VX2,k,j,2*ighost-i+sign);
                        Vc(VX3,k,j,i) = Vk*sqrt((pexp+qexp)*h0*h0*pow(R,2*flaringIndex)+(1+qexp)-qexp*R/r) - omega * R;

//                         Vc(RHO,k,j,i) = sigma0/sqrt(2.0*M_PI)/(h0)*pow(R,-sigmaSlope-1.5) * exp(1.0/cs2 * (1.0/r - 1.0/R)) ;
//                         Vc(VX1,k,j,i) = - Vc(VX1,k,j,2*ighost-i+1);
//                         Vc(VX2,k,j,i) = ZERO_F;
//                         Vc(VX3,k,j,i) = Vk*sqrt(1 - (1.5+sigmaSlope)*cs2*R);

                    });
    }

    if( dir==JDIR) {
        IdefixArray4D<real> Vc = hydro->Vc;
        int jghost;
        int jbeg,jend;
        // UPPER LAYER
        if(side == left) {
            jghost = data->beg[JDIR];
            jbeg = 0;
            jend = data->beg[JDIR];
            //return;
            idefix_for("UserDefBoundary",
              0, data->np_tot[KDIR],
              jbeg, jend,
              0, data->np_tot[IDIR],
                        KOKKOS_LAMBDA (int k, int j, int i) {
                            real r = x1(i);
                            real R = x1(i)*sin(x2(j));
                            real Rghost = x1(i)*sin(x2(jghost));
                            real Vk = 1.0/sqrt(R);
                            real cs2 = h0*h0*pow(R,2*flaringIndex-1);

                            Vc(RHO,k,j,i) = Vc(RHO,k,jghost,i) * pow(R/Rghost,-sigmaSlope-flaringIndex-1) * exp(1.0/(cs2)*(1.0/Rghost - 1.0/R)) ;
                            // Vc(RHO,k,j,i) = sigma0/sqrt(2.0*M_PI)/(h0)*pow(R,-sigmaSlope-1.5) * exp(1.0/(cs2) * (1.0/r - 1.0/R)) ;
                            Vc(VX1,k,j,i) = ZERO_F; //Vc(VX1,k,2*jghost-j-1,i);
                            Vc(VX2,k,j,i) = ZERO_F;
                            Vc(VX3,k,j,i) = Vk*sqrt((pexp+qexp)*h0*h0*pow(R,2*flaringIndex)+(1+qexp)-qexp*R/r) - omega * R;
                        });
        }
        // MIDPLANE
        else if(side==right) {
            jghost = data->end[JDIR]-1;
            jbeg = data->end[JDIR];
            jend = data->np_tot[JDIR];
            idefix_for("UserDefBoundary",
              0, data->np_tot[KDIR],
              jbeg, jend,
              0, data->np_tot[IDIR],
                        KOKKOS_LAMBDA (int k, int j, int i) {
                          real r = x1(i);
                          real R = x1(i)*sin(x2(j));
                          real Vk = 1.0/sqrt(R);

                          Vc(RHO,k,j,i) = Vc(RHO,k,2*jghost-j+1,i);
                          Vc(VX1,k,j,i) = Vc(VX1,k,2*jghost-j+1,i);
                          Vc(VX2,k,j,i) = - Vc(VX2,k,2*jghost-j+1,i);
                          Vc(VX3,k,j,i) = Vc(VX3,k,2*jghost-j+1,i);
                        });
        }
    }

}

void FargoVelocity(DataBlock &data, IdefixArray2D<real> &Vphi) {
  IdefixArray1D<real> x1 = data.x[IDIR];
  IdefixArray1D<real> x2 = data.x[JDIR];
  real h0=h0Glob;
  real flaringIndex=flaringIndexGlob;
  real sigmaSlope=sigmaSlopeGlob;
  real pexp = -sigmaSlope-flaringIndex-1;
  real qexp = 2*flaringIndex-1;
  real omega = omegaGlob;
  idefix_for("FargoVphi",0,data.np_tot[JDIR], 0, data.np_tot[IDIR],
      KOKKOS_LAMBDA (int j, int i) {
        real r = x1(i);
        real R = x1(i)*sin(x2(j));
        real Vk = 1.0/sqrt(R);
        Vphi(j,i) = Vk*sqrt((pexp+qexp)*h0*h0*pow(R,2*flaringIndex)+(1+qexp)-qexp*R/r) - omega * R;
  });
}

// Analyse data to produce an ascii output
void Analysis(DataBlock & data) {
  // Mirror data on Host
  DataBlockHost d(data);

  // Sync it
  d.SyncFromDevice();

  for(int ip=0; ip < data.planetarySystem->nbp ; ip++) {
    // Get the orbital parameters
    real timeStep = data.dt;
    real xp = data.planetarySystem->planet[ip].getXp();
    real yp = data.planetarySystem->planet[ip].getYp();
    real zp = data.planetarySystem->planet[ip].getZp();
    real vxp = data.planetarySystem->planet[ip].getVxp();
    real vyp = data.planetarySystem->planet[ip].getVyp();
    real vzp = data.planetarySystem->planet[ip].getVzp();
    real qp = data.planetarySystem->planet[ip].getMp();
    real time = data.t;

    std::string planetName, tqwkName;
    std::stringstream pName, tName;
    pName << "planet" << ip << ".dat";
    tName << "tqwk" << ip << ".dat";
    planetName = pName.str();
    tqwkName = tName.str();
    // Write the data in ascii to our file
    if(idfx::prank==0) {
      std::ofstream f;
      f.open(planetName,std::ios::app);
      f.precision(10);
      f << std::scientific << timeStep << "    " << xp << "    " << yp << "    " << zp << "    " << vxp << "    " << vyp << "    " << vzp << "    " << qp << "    " << time << std::endl;
      f.close();
    }

    // Force force = {{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}};
    Force &force = data.planetarySystem->planet[ip].m_force;
    bool isp = true;
    data.planetarySystem->planet[ip].computeForce(data,isp);

    // Get the torque and work
    // real fxi = force.f_inner[0];
    // real fyi = force.f_inner[1];
    // real fzi = force.f_inner[2];
    // real fxo = force.f_outer[0];
    // real fyo = force.f_outer[1];
    // real fzo = force.f_outer[2];
    // real fxhi = force.f_ex_inner[0];
    // real fyhi = force.f_ex_inner[1];
    // real fzhi = force.f_ex_inner[2];
    // real fxho = force.f_ex_outer[0];
    // real fyho = force.f_ex_outer[1];
    // real fzho = force.f_ex_outer[2];

    real tq_inner = xp*force.f_inner[1]-yp*force.f_inner[0];
    real tq_outer = xp*force.f_outer[1]-yp*force.f_outer[0];
    real tq_ex_inner = xp*force.f_ex_inner[1]-yp*force.f_ex_inner[0];
    real tq_ex_outer = xp*force.f_ex_outer[1]-yp*force.f_ex_outer[0];
    real wk_inner = vxp*force.f_inner[0]+vyp*force.f_inner[1];
    real wk_outer = vxp*force.f_outer[0]+vyp*force.f_outer[1];
    real wk_ex_inner = vxp*force.f_ex_inner[0]+vyp*force.f_ex_inner[1];
    real wk_ex_outer = vxp*force.f_ex_outer[0]+vyp*force.f_ex_outer[1];

    // Write the data in ascii to our file
    if(idfx::prank==0) {
      std::ofstream ft;
      ft.open(tqwkName,std::ios::app);
      ft.precision(10);
      // ft << std::scientific << timeStep << "    " << fxi << "    " << fyi << "    " << fzi << "    " << fxo << "    " << fyo << "    " << fzo << "    " << fxhi << "    " << fyhi << "    " << fzhi << "    " << fxho << "    " << fyho << "    " << fzho << "    " << time << std::endl;
      ft << std::scientific << timeStep << "    " << tq_inner << "    " << tq_outer << "    " << tq_ex_inner << "    " << tq_ex_outer << "    " << wk_inner << "    " << wk_outer << "    " << wk_ex_inner << "    " << wk_ex_outer << "    " << time << std::endl;
      ft.close();
    }
  }
}

// Default constructor
// Initialisation routine. Can be used to allocate
// Arrays or variables which are used later on
Setup::Setup(Input &input, Grid &grid, DataBlock &data, Output &output)// : m_planet(0)//, Planet &planet)
{
  // Set the function for userdefboundary
  data.hydro->EnrollUserDefBoundary(&UserdefBoundary);
  data.hydro->EnrollUserSourceTerm(&Damping);
  data.hydro->EnrollIsoSoundSpeed(&MySoundSpeed);
  if(data.haveFargo)
    data.fargo->EnrollVelocity(&FargoVelocity);
  if(data.hydro->haveRotation) {
    omegaGlob = data.hydro->OmegaZ;
  } else {
    omegaGlob = 0.0;
  }
  // Enroll the analysis function
  output.EnrollAnalysis(&Analysis);
  wkzMinGlob = input.Get<real>("Setup","wkzMin",0);
  wkzMaxGlob = input.Get<real>("Setup","wkzMax",0);
  wkzDampingGlob = input.Get<real>("Setup","wkzDamping",0);
  sigma0Glob = input.Get<real>("Setup","sigma0",0);
  sigmaSlopeGlob = input.Get<real>("Setup","sigmaSlope",0);
  h0Glob = input.Get<real>("Setup","h0",0);
  flaringIndexGlob = input.Get<real>("Setup","flaringIndex",0);
  densityFloorGlob = input.Get<real>("Setup","densityFloor",0);
  masstaperGlob = input.Get<real>("Planet","masstaper",0);
  // delete file planet0.dat at initialization if we do not restart the simulation.
  for(int ip=0; ip < data.planetarySystem->nbp ; ip++) {
    std::string planetName, tqwkName;
    std::stringstream pName, tName;
    pName << "planet" << ip << ".dat";
    tName << "tqwk" << ip << ".dat";
    planetName = pName.str();
    tqwkName = tName.str();
    if (!(input.restartRequested)) {
      if(idfx::prank==0) {
        std::ofstream f;
        f.open(planetName,std::ios::out);
        f.close();
        std::ofstream ft;
        ft.open(tqwkName,std::ios::out);
        ft.close();
      }
    }
    /*
    else {
      if(idfx::prank==0) {
        char ch;
        int nline = input.restartFileNumber*input.Get<real>("Output","dmp",0)/input.Get<real>("Output","analysis",0);
        int line = 1;

        // planet0.dat
        std::ifstream fin(planetName);
        std::ofstream fout;
        fout.open("temp.dat", std::ios::out);
        while(fin.get(ch))
        {
          if(ch == '\n')
            line++;

          if(line <= nline)      // content not to be deleted
            fout<<ch;
        }
        fout << std::endl;
        fout.close();
        fin.close();
        std::remove(planetName.c_str());
        std::rename("temp.dat",planetName.c_str());

        // tqwk0.dat
        std::ifstream fint(tqwkName);
        std::ofstream foutt;
        foutt.open("tempt.dat", std::ios::out);
        nline = input.restartFileNumber*input.Get<real>("Output","dmp",0)/input.Get<real>("Output","analysis",0);
        line = 1;
        while(fint.get(ch))
        {
          if(ch == '\n')
            line++;

          if(line <= nline)      // content not to be deleted
            foutt<<ch;
        }

        foutt << std::endl;
        foutt.close();
        fint.close();
        std::remove(tqwkName.c_str());
        std::rename("tempt.dat",tqwkName.c_str());
      }
    }
    */
  }
}


// This routine initialize the flow
// Note that data is on the device.
// One can therefore define locally
// a datahost and sync it, if needed
void Setup::InitFlow(DataBlock &data) {
    // Create a host copy
    DataBlockHost d(data);
    real h0=h0Glob;
    real flaringIndex=flaringIndexGlob;
    real sigma0=sigma0Glob;
    real sigmaSlope=sigmaSlopeGlob;
    real pexp = -sigmaSlope-flaringIndex-1;
    real qexp = 2*flaringIndex-1;

    for(int k = 0; k < d.np_tot[KDIR] ; k++) {
        for(int j = 0; j < d.np_tot[JDIR] ; j++) {
            for(int i = 0; i < d.np_tot[IDIR] ; i++) {
                real r = d.x[IDIR](i);
                real R = d.x[IDIR](i)*sin(d.x[JDIR](j));
                real Vk = 1.0/sqrt(R);

                real cs2 = h0*h0*pow(R,2*flaringIndex-1);

                d.Vc(RHO,k,j,i) = sigma0/sqrt(2.0*M_PI)/(h0)*pow(R,-sigmaSlope-flaringIndex-1) * exp(1.0/cs2 * (1.0/r - 1.0/R)) ;
                d.Vc(VX1,k,j,i) = ZERO_F;
                d.Vc(VX2,k,j,i) = ZERO_F;
                d.Vc(VX3,k,j,i) = Vk*sqrt((pexp+qexp)*h0*h0*pow(R,2*flaringIndex)+(1+qexp)-qexp*R/r) - omegaGlob * R;
            }
        }
    }

    // Send it all, if needed
    d.SyncToDevice();
}
