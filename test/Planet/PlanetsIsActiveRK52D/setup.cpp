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
real omegaGlob;


void MySoundSpeed(DataBlock &data, const real t, IdefixArray3D<real> &cs) {
  real h0 = h0Glob;
  real flaringIndex = flaringIndexGlob;
  IdefixArray1D<real> x1=data.x[IDIR];
  idefix_for("MySoundSpeed",0,data.np_tot[KDIR],0,data.np_tot[JDIR],0,data.np_tot[IDIR],
              KOKKOS_LAMBDA (int k, int j, int i) {
                real R = x1(i);
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
  real dt = dtin;
  bool isFargo = data->haveFargo;

  real rmin{data->mygrid->xbeg[0]};
  real wkzMin{wkzMinGlob};
  real rmax{data->mygrid->xend[0]};
  real wkzMax{wkzMaxGlob};
  real wkzDamping{wkzDampingGlob};
  idefix_for("MySourceTerm",
    0, data->np_tot[KDIR],
    0, data->np_tot[JDIR],
    0, data->np_tot[IDIR],
              KOKKOS_LAMBDA (int k, int j, int i) {
                real R = x1(i);
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
                // Damp whatever is at R<wkzMin and R>wkzMax
                if (R<wkzMin) {
                  lambda = 1.0/(wkzDamping*2.0*M_PI*pow(rmin,1.5))*(1.0 - pow(sin(M_PI*( (R-rmin) / (wkzMin-rmin) )/2.0),2.0));
                }
                if (R>wkzMax) {
                  lambda = 1.0/(wkzDamping*2.0*M_PI*pow(rmax,1.5))*pow(sin(M_PI*( (R-wkzMax) / (rmax-wkzMax) )/2.0),2.0);
                }

                real rhoTarget = sigma0*pow(R,-sigmaSlope) ;
                real vx2Target = 0.0;
                if(!isFargo) {
                  vx2Target = Vk*sqrt(1.0-(1.0+sigmaSlope-2*flaringIndex)*h0*h0*pow(R,2*flaringIndex)) - omega * R;
                }

                // relaxation
                real drho = lambda*(Vc(RHO,k,j,i)-rhoTarget);
                real dvx1 = lambda*Vc(RHO,k,j,i)*Vc(VX1,k,j,i);
                real dvx2 = lambda*Vc(RHO,k,j,i)*(Vc(VX2,k,j,i)-vx2Target);
                real dvx3 = lambda*Vc(RHO,k,j,i)*Vc(VX3,k,j,i);

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
  DataBlock &data = *hydro->data;
  IdefixArray4D<real> Vc = hydro->Vc;
  IdefixArray1D<real> x1 = data.x[IDIR];
  real sigmaSlope=sigmaSlopeGlob;
  real omega = omegaGlob;

    if(dir==IDIR) {
        int ighost,ibeg,iend,sign;
        if(side == left) {
            ighost = data.beg[IDIR];
            ibeg = 0;
            iend = data.beg[IDIR];
            sign=-1;
            //return;
        }
        else if(side==right) {
            ighost = data.end[IDIR]-1;
            ibeg=data.end[IDIR];
            iend=data.np_tot[IDIR];
            sign=1;
        }


        idefix_for("UserDefBoundary",0,data.np_tot[KDIR],0,data.np_tot[JDIR],ibeg,iend,
                    KOKKOS_LAMBDA (int k, int j, int i) {
                        real R=x1(i);
                        real R0=x1(ighost);
                        real Vk = 1.0/sqrt(R);

                        Vc(RHO,k,j,i) = Vc(RHO,k,j,ighost)*pow(R/R0,-sigmaSlope) ;
                        Vc(VX1,k,j,i) = - Vc(VX1,k,j,2*ighost-i+sign);
                        Vc(VX2,k,j,i) = Vc(VX2,k,j,ighost)*pow(R/R0,-0.5) - omega*R ;
                        Vc(VX3,k,j,i) = Vc(VX3,k,j,2*ighost-i+sign);
                    });
    }
}

void FargoVelocity(DataBlock &data, IdefixArray2D<real> &Vphi) {
  IdefixArray1D<real> x1 = data.x[IDIR];
  real omega = omegaGlob;

  idefix_for("FargoVphi",0,data.np_tot[KDIR], 0, data.np_tot[IDIR],
      KOKKOS_LAMBDA (int k, int i) {
        real R = x1(i);
        Vphi(k,i) = 1.0/sqrt(R) - omega*R;
  });
}

// Analyse data to produce an ascii output
void Analysis(DataBlock & data) {
  // Mirror data on Host
  DataBlockHost d(data);

  // Sync it
  d.SyncFromDevice();

//  data.DumpToFile("totou");

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
    bool isActive = data.planetarySystem->planet[ip].getIsActive();
    // idfx::cout << "isactive: " << isActive << std::endl;

    std::string isActiveName;
    std::stringstream iaName;
    iaName << "isactive" << ip << ".dat";
    isActiveName = iaName.str();
    // Write the data in ascii to our file
    if(idfx::prank==0) {
      std::ofstream fi;
      fi.open(isActiveName,std::ios::app);
      fi.precision(10);
      fi << std::scientific << time << "    " << isActive << std::endl;
      fi.close();
    }

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


// Compute user variables which will be written in vtk files
void ComputeUserVars(DataBlock & data, UserDefVariablesContainer &variables) {
  // Mirror data on Host
  real h0=h0Glob;
  real flaringIndex=flaringIndexGlob;
  DataBlockHost d(data);

  // Sync it
  d.SyncFromDevice();

  // Make references to the user-defined arrays (variables is a container of IdefixHostArray3D)
  // Note that the labels should match the variable names in the input file
  IdefixHostArray3D<real> press = variables["PRS"];

  for(int k = 0; k < d.np_tot[KDIR] ; k++) {
    for(int j = 0; j < d.np_tot[JDIR] ; j++) {
      for(int i = 0; i < d.np_tot[IDIR] ; i++) {
        real R = d.x[IDIR](i);
        real cs2 = h0*h0*pow(R,2*flaringIndex-1.);
        press(k,j,i) = cs2*d.Vc(RHO,k,j,i);
      }
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
  // Enroll our user-defined variables
  output.EnrollUserDefVariables(&ComputeUserVars);
  wkzMinGlob = input.Get<real>("Setup","wkzMin",0);
  wkzMaxGlob = input.Get<real>("Setup","wkzMax",0);
  wkzDampingGlob = input.Get<real>("Setup","wkzDamping",0);
  sigma0Glob = input.Get<real>("Setup","sigma0",0);
  sigmaSlopeGlob = input.Get<real>("Setup","sigmaSlope",0);
  h0Glob = input.Get<real>("Setup","h0",0);
  flaringIndexGlob = input.Get<real>("Setup","flaringIndex",0);
  densityFloorGlob = input.Get<real>("Setup","densityFloor",0);
  // delete file planet0.dat at initialization if we do not restart the simulation.
  for(int ip=0; ip < data.planetarySystem->nbp ; ip++) {
    std::string planetName, tqwkName, isActiveName;
    std::stringstream pName, tName, iaName;
    pName << "planet" << ip << ".dat";
    tName << "tqwk" << ip << ".dat";
    iaName << "isactive" << ip << ".dat";
    planetName = pName.str();
    tqwkName = tName.str();
    isActiveName = iaName.str();
    if (!(input.restartRequested)) {
      if(idfx::prank==0) {
        std::ofstream f;
        f.open(planetName,std::ios::out);
        f.close();
        std::ofstream ft;
        ft.open(tqwkName,std::ios::out);
        ft.close();
        std::ofstream fia;
        fia.open(isActiveName,std::ios::out);
        fia.close();
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
    }*/
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
    real sigma0=sigma0Glob;
    real sigmaSlope=sigmaSlopeGlob;
    real flaringIndex=flaringIndexGlob;
    real omega=omegaGlob;

    for(int k = 0; k < d.np_tot[KDIR] ; k++) {
        for(int j = 0; j < d.np_tot[JDIR] ; j++) {
            for(int i = 0; i < d.np_tot[IDIR] ; i++) {
                real R = d.x[IDIR](i);
                real z = d.x[KDIR](k);
                real Vk = 1.0/sqrt(R);

                d.Vc(RHO,k,j,i) = sigma0*pow(R,-sigmaSlope) ;
                d.Vc(VX1,k,j,i) = 0.0;
                d.Vc(VX2,k,j,i) = Vk*sqrt(1.0-(1.0+sigmaSlope-2*flaringIndex)*h0*h0*pow(R,2*flaringIndex)) - omega * R;
                d.Vc(VX3,k,j,i) = 0.0;

                if(d.Vc(RHO,k,j,i) < densityFloorGlob) {
                  d.Vc(RHO,k,j,i) = densityFloorGlob;
                }
            }
        }
    }

    // Send it all, if needed
    d.SyncToDevice();
}

// Analyse data to produce an output
void MakeAnalysis(DataBlock & data) {

}
