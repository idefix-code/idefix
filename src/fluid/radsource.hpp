// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************
#ifndef FLUID_RADSOURCE_HPP_
#define FLUID_RADSOURCE_HPP_

#include <string>
#include "idefix.hpp"
#include "input.hpp"
#include "fluid_defs.hpp"
#include "eos.hpp"
#include "units.hpp"
#include "lookupTable.hpp"
#include "column.hpp"

class RadSource {
 public:
  // Type of opacity
  enum class Opacity{constant,kramers,usertable,userfunc};
  // Type of irradiation flux
  enum class IradiationFlux{constant,usertable,userfunc,usergeometry};
  // Type of implicit solver for radiation source terms
  enum class ImplicitSolver{full_implicit,fixed_point_rad,fixed_point_gas};

  // RadSource constructor
  template <typename Phys>
  RadSource(Input &, Fluid<Phys> *);

  void ShowConfig();                    // print configuration
  void AddRadSource(const real);        // Effectively add the radiation source terms

  // Type of implicit solver for radiation source terms
  void SourceFullImplicit(const real);
  void SourceFixedPointRad(const real);
  void SourceFixedPointGas(const real);

  // Add Relativistic Corrections in an explicit way as in Melon & Fuksman & Klahr 2022
  void RelativistCorrection(const real);

  // Compute divergence of external irradiation flux
  void IrrFlux(IdefixArray3D<real>);

  // Compute viscous heating
  void Qviscous(IdefixArray3D<real>);

  // Arrays containing the userfunc opacities (copied from Fluid class)
  IdefixArray3D<real> xiArr;
  IdefixArray3D<real> kappapArr;
  IdefixArray3D<real> kapparArr;
  IdefixArray3D<real> kappairrArr;

  // Arrays containing the irradiation field (copied from Fluid class)
  IdefixArray3D<real> irrArr;

  // Array containing kappa*rho for userfunc irradiation flux
  IdefixArray3D<real> kapparhoArr;

  IdefixArray4D<real> UcRad;  // Radiation conservative quantities
  IdefixArray4D<real> UcGas;  // Gas conservative quantities
  IdefixArray4D<real> VcRad;  // Radiation primitive quantities
  IdefixArray4D<real> VcGas;  // Gas primitive quantities
  IdefixArray3D<real> InvDt;  // The InvDt of current radiation multigroup

  // Data related to current instance of the Rad object
  int instanceNumber;

  // Compute limiting diffusion speed for Riemann solver in opt. thick media
  KOKKOS_INLINE_FUNCTION real LimitSpeedsRad(int i, int j, int k, real dx) const {
    auto VcGas = this->VcGas;

    real kappa,xi;
    real mu = eos.GetMu(VcGas(PRS,k,j,i),VcGas(RHO,k,j,i));

    if (kappa_type == Opacity::constant) {
      kappa = this->kappar_0;
    } else if (kappa_type == Opacity::kramers) {
      real T = VcGas(PRS,k,j,i)/(VcGas(RHO,k,j,i))*this->unit_Kelvin*mu;
      kappa = this->kappar_0*VcGas(RHO,k,j,i)*this->unit_density/this->rho_0;
      kappa *= std::pow(T/this->T_0,-3.5);
    } else if (kappa_type == Opacity::usertable) {
      real T = VcGas(PRS,k,j,i)/(VcGas(RHO,k,j,i))*this->unit_Kelvin*mu;
      real logT = std::log10(T);
      real logrho = std::log10(VcGas(RHO,k,j,i)*this->unit_density);
      if (this->kappa_ndim == 1) {
        kappa = this->kappa_ross_1D.Get(&logT);
      } else if (this->kappa_ndim == 2) {
        real x[2];
        x[1] = FMIN(FMAX(logT,2.5),5.98);
        x[0] = FMIN(-4.05,FMAX(-14.,logrho));
        kappa = std::pow(10.,this->kappa_ross_2D.Get(x));
      }
    } else if (kappa_type == Opacity::userfunc) {
      kappa = this->kapparArr(k,j,i);
    }

    if (xi_type == Opacity::constant) {
      xi = this->xi_0;
    } else if (xi_type == Opacity::usertable) {
      real T = VcGas(PRS,k,j,i)/(VcGas(RHO,k,j,i))*this->unit_Kelvin*mu;
      real logT = std::log10(T);
      real logrho = std::log10(VcGas(RHO,k,j,i)*this->unit_density);
      if (this->xi_ndim == 1) {
        xi = this->xi_1D.Get(&logT);
      } else if (this->xi_ndim == 2) {
        real x[2];
        x[1] = logT;
        x[0] = logrho;
        xi = std::pow(10.,this->xi_2D.Get(x));
      }
    } else if (xi_type == Opacity::userfunc) {
      xi = this->xiArr(k,j,i);
    }


    // Compute optical depth across one cell
    real tau = VcGas(RHO,k,j,i)*this->unit_density*(kappa+xi)*dx*this->unit_length;

    // return characteristic velocity of radiative diffusion
    return 4./(3.*tau)*this->reduced_c;
  }

 private:
  DataBlock* data;
  real kappap_0;
  real kappar_0;
  real kappap_es;
  real kappar_es;
  real xi_0;
  real kappa_irr;
  real rho_0;
  real T_0;
  real rs;
  real Ts;
  real kappa_star;
  real reduced_c;
  real gamma;
  real mu;
  int count_max;

  // EOS
  EquationOfState eos;

  // Dimension of Planck and Rosseland opacities tables
  int kappa_ndim;
  int xi_ndim;

  // User-defined opacity tables
  LookupTable<1> kappa_planck_1D;
  LookupTable<1> kappa_ross_1D;
  LookupTable<1> xi_1D;
  LookupTable<2> kappa_planck_2D;
  LookupTable<2> kappa_ross_2D;
  LookupTable<2> xi_2D;

  // Have irradiation or not
  bool haveIrradiation{false};

  // Have relativist correction or not
  bool haveRelativistCorrection{false};

  // Dimension of irradiation flux table
  int irr_ndim;

  // Irradiation flux table
  LookupTable<1> irr_1D;

  Column *column_rho;        // Column density
  IdefixArray3D<real> divF;  // Divergence of irradiation flux
  IdefixArray3D<real> Qvisc;  // Viscous heating

  ImplicitSolver source_solver;    // Type of implicit solver for radiation source terms
  Opacity kappa_type;          // Type of absorption opacity definition
  Opacity xi_type;             // Type of scattering opacity definition
  IradiationFlux irr_type;             // Type of irradiation flux definition

  //Units
  real unit_density = idfx::units.GetDensity();
  real unit_length = idfx::units.GetLength();
  real unit_Kelvin = idfx::units.GetKelvin();
};

#include "fluid.hpp"

template<typename Phys>
RadSource::RadSource(Input &input, Fluid<Phys> *hydroin):
                      UcRad{hydroin->Uc},
                      UcGas{hydroin->data->hydro->Uc},
                      VcRad{hydroin->Vc},
                      VcGas{hydroin->data->hydro->Vc},
                      InvDt{hydroin->InvDt},
                      eos{*(hydroin->data->hydro->eos.get())} {
  idfx::pushRegion("RadSource::RadSource");

  //do not yet support ORDER=3 with Radiation support
  if (ORDER == 3) {
    IDEFIX_ERROR("ORDER=3 is not yet valid when used in conjunction with Radiation support !");
  }

  // Create our own prefix
  std::string prefix = std::string(Phys::prefix);

  // Save the parent hydro object
  this->data = hydroin->data;

  // Check in which block we should fetch our information
  std::string BlockName;
  if(Phys::radiation) {
    BlockName = "Rad";
  } else {
    IDEFIX_ERROR("Fluid is not radiative");
  }

  // Reduced velocity of light
  this->reduced_c =  hydroin->reduced_c;

  // Information on scattering opacity coefficient
  if(input.CheckEntry(BlockName,"xi")>=0) {
    // Fetch the opacity coefficient for the current radiation group.
    const int n = hydroin->instanceNumber;

    std::string xiType = input.Get<std::string>(BlockName,"xi",0);
    if(xiType.compare("constant") == 0) {
      this->xi_type = Opacity::constant;
      this->xi_0 = input.Get<real>(BlockName,"xi",n+1);
    } else if(xiType.compare("usertable") == 0) {
      this->xi_type = Opacity::usertable;
      this->xi_ndim = input.Get<int>(BlockName,"xi",n+1);
      std::string xi_file = input.Get<std::string>(BlockName,"xi",n+2);
      if (this->xi_ndim == 1) {
        this->xi_1D = LookupTable<1>(xi_file,',');
      } else if (this->xi_ndim == 2) {
        this->xi_2D = LookupTable<2>(xi_file,',');
      } else {
        std::stringstream msg;
        msg << "Only 1 or 2 dimension for scattering opacity tables"
               "are currently accepted." << std::endl;
        IDEFIX_ERROR(msg);
      }
    } else if (xiType.compare("userfunc") == 0) {
      this->xi_type = Opacity::userfunc;
      this->xiArr = hydroin->xiArr;
    } else {
      std::stringstream msg;
      msg << "Unknown xi type \"" <<  xiType
          << "\" in your input file." << std::endl
          << "Allowed values are: constant, usertable." << std::endl;
      IDEFIX_ERROR(msg);
    }
  } else {
    IDEFIX_ERROR("A *xi* line in your [Rad] block is required"
                 "in your input file to define the scattering opacity.");
  }

  // Information on absorption opacity coefficient
  if(input.CheckEntry(BlockName,"kappa")>=0) {
    // Fetch the opacity coefficient for the current radiation group.
    const int n = hydroin->instanceNumber;

    std::string kappaType = input.Get<std::string>(BlockName,"kappa",0);
    if(kappaType.compare("constant") == 0) {
      this->kappa_type = Opacity::constant;
      this->kappap_0 = input.Get<real>(BlockName,"kappa",n+1);
      this->kappar_0 = input.Get<real>(BlockName,"kappa",n+2);
    } else if(kappaType.compare("kramers") == 0) {
      this->kappap_0 = input.Get<real>(BlockName,"kappa",n+1);
      this->kappar_0 = input.Get<real>(BlockName,"kappa",n+2);
      this->kappa_type = Opacity::kramers;
      this->rho_0 = input.Get<real>(BlockName,"kappa",n+3);
      this->T_0 = input.Get<real>(BlockName,"kappa",n+4);
      this->kappap_es = input.Get<real>(BlockName,"kappa",n+5);
      this->kappar_es = input.Get<real>(BlockName,"kappa",n+6);
    } else if(kappaType.compare("usertable") == 0) {
      this->kappa_type = Opacity::usertable;
      this->kappa_ndim = input.Get<int>(BlockName,"kappa",n+1);
      std::string kappap_file = input.Get<std::string>(BlockName,"kappa",n+2);
      std::string kappar_file = input.Get<std::string>(BlockName,"kappa",n+3);
      if (this->kappa_ndim == 1) {
        this->kappa_planck_1D = LookupTable<1>(kappap_file,',');
        this->kappa_ross_1D = LookupTable<1>(kappar_file,',');
      } else if (this->kappa_ndim == 2) {
        this->kappa_planck_2D = LookupTable<2>(kappap_file,',');
        this->kappa_ross_2D = LookupTable<2>(kappar_file,',');
      } else {
        std::stringstream msg;
        msg << "Only 1 or 2 dimensions for absorption opacity tables"
               "are currently accepted." << std::endl;
        IDEFIX_ERROR(msg);
      }
    } else if (kappaType.compare("userfunc") == 0) {
      this->kappa_type = Opacity::userfunc;
      this->kappapArr = hydroin->kappapArr;
      this->kapparArr = hydroin->kapparArr;
    } else {
      std::stringstream msg;
      msg << "Unknown kappa type \"" <<  kappaType
          << "\" in your input file." << std::endl
          << "Allowed values are: constant, kramers, usertable, userfunc" << std::endl;

      IDEFIX_ERROR(msg);
    }
  } else {
    IDEFIX_ERROR("A *kappa* line in your [Rad] block is required"
                 "in your input file to define the absorption opacity.");
  }

  // Information on solver for source terms
  if(input.CheckEntry(BlockName,"source")>=0) {
    // Fetch the opacity coefficient for the current radiation group.
    const int n = hydroin->instanceNumber;

    std::string sourceType = input.Get<std::string>(BlockName,"source",0);
    if(sourceType.compare("full_implicit") == 0) {
      this->source_solver = ImplicitSolver::full_implicit;
    } else if(sourceType.compare("fixed_point_rad") == 0) {
      this->source_solver = ImplicitSolver::fixed_point_rad;
    } else if(sourceType.compare("fixed_point_gas") == 0) {
      this->source_solver = ImplicitSolver::fixed_point_gas;
    } else {
      std::stringstream msg;
      msg << "Unknown solver for source terms \"" <<  sourceType
          << "\" in your input file." << std::endl
          << "Allowed values are: full_implicit, fixed_point_rad, fixed_point_gas." << std::endl;

      IDEFIX_ERROR(msg);
    }

  } else {
    IDEFIX_ERROR("A *source* line in your [Rad] block is required"
                 "in your input file to define the solver for the radiation source terms.");
  }

  // Information on relativist correction
  haveRelativistCorrection = input.GetOrSet<bool>(BlockName,"relativist_correction",0,false);

  // Information on irradiation source term
  if(input.CheckEntry(BlockName,"irr")>=0) {
    haveIrradiation = true;
    // Fetch the opacity coefficient for the current radiation group.
    const int n = hydroin->instanceNumber;

    this->column_rho = new Column(IDIR,1,data);
    this->divF = IdefixArray3D<real>(prefix+"_divF",data->np_tot[KDIR],
                                     data->np_tot[JDIR], data->np_tot[IDIR]);
    this->Qvisc = IdefixArray3D<real>(prefix+"_Qvisc",data->np_tot[KDIR],
                                      data->np_tot[JDIR], data->np_tot[IDIR]);


    std::string irrType = input.Get<std::string>(BlockName,"irr",0);

    if(irrType.compare("constant") == 0) {
      this->irr_type = IradiationFlux::constant;
      this->rs = input.Get<real>(BlockName,"irr",n+1);
      this->Ts = input.Get<real>(BlockName,"irr",n+2);
      this->kappa_irr = input.Get<real>(BlockName,"irr",n+3);
    } else if(irrType.compare("usertable") == 0) {
      this->irr_type = IradiationFlux::usertable;
      this->rs = input.Get<real>(BlockName,"irr",n+1);
      this->Ts = input.Get<real>(BlockName,"irr",n+2);
      this->irr_ndim = input.Get<int>(BlockName,"irr",n+3);
      std::string irr_file = input.Get<std::string>(BlockName,"irr",n+4);
      if (input.Get<int>(BlockName,"irr",n+3) == 1) {
        this->irr_1D = LookupTable<1>(irr_file,',');
      } else {
        std::stringstream msg;
        msg << "Only 1 dimension for irradiation flux tables are currently accepted." << std::endl;
        IDEFIX_ERROR(msg);
      }
    } else if(irrType.compare("userfunc") == 0) {
      this->irr_type = IradiationFlux::userfunc;
      this->rs = input.Get<real>(BlockName,"irr",n+1);
      this->Ts = input.Get<real>(BlockName,"irr",n+2);
      this->kappa_star = input.Get<real>(BlockName,"irr",n+3);
      this->kapparhoArr = IdefixArray3D<real>("kapparrhoArray",data->np_tot[KDIR],
                                                 data->np_tot[JDIR],
                                                 data->np_tot[IDIR]);
      this->kappairrArr = hydroin->kappairrArr;
    } else if(irrType.compare("usergeometry") == 0) {
      this->irr_type = IradiationFlux::usergeometry;
      this->irrArr = hydroin->irrArr;
    } else {
      std::stringstream msg;
      msg << "Unknown irr type \"" <<  irrType
          << "\" in your input file." << std::endl
          << "Allowed values are: constant, usertable, userfunc, usergeometry." << std::endl;

      IDEFIX_ERROR(msg);
    }
  }



  idfx::popRegion();
}

#endif // FLUID_RADSOURCE_HPP_
