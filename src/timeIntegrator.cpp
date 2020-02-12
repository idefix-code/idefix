#include <cstdio>
#include "idefix.hpp"
#include "timeIntegrator.hpp"

TimeIntegrator::TimeIntegrator(Input & input, DataBlock &datain) {
    this->data=datain;
    nstages=input.nstages;
    if(nstages>1) {
        // Temporary array to store initial field in the RK2-3 loops
        V0 = IdefixArray4D<real>("TimeIntegrator_V0", NVAR, data.np_tot[KDIR], data.np_tot[JDIR], data.np_tot[IDIR]);
    }
    if(nstages==2) {
        wc[0] = 0.5;
        w0[0] = 0.5;
    }
    if(nstages==3) {
        wc[0] = 0.25;
        w0[0] = 0.5;
        wc[1] = 2.0/3.0;
        w0[1] = 1.0/3.0;
    }
    InvDtHyp = IdefixArray3D<real>("TimeIntegrator_InvDtHyp", data.np_tot[KDIR], data.np_tot[JDIR], data.np_tot[IDIR]);
    InvDtPar = IdefixArray3D<real>("TimeIntegrator_InvDtPar", data.np_tot[KDIR], data.np_tot[JDIR], data.np_tot[IDIR]);
}

// Compute one Stage of the time Integrator
void TimeIntegrator::Stage() {}


// Compute one full cycle of the time Integrator
void TimeIntegrator::Cycle() {
    // Do one cycle
    IdefixArray4D<real> Vc = data.Vc;
    IdefixArray4D<real> V0 = this->V0;
    IdefixArray3D<real> InvDtHypLoc=this->InvDtHyp;
    IdefixArray3D<real> InvDtParLoc=this->InvDtPar;

    // Store initial stage for multi-stage time integrators
    if(nstages>1) Kokkos::deep_copy(V0,Vc);

    for(int stage=0; stage < nstages ; stage++) {
        // Update Vc
        Stage();
        // Is this the last stage?
        if(stage<nstages-1) {
            // No!
            real wcs=wc[stage];
            real w0s=w0[stage];

            idefix_for("Cycle-update",0,NVAR,0,data.np_tot[KDIR],0,data.np_tot[JDIR],0,data.np_tot[IDIR],
                        KOKKOS_LAMBDA (int n, int k, int j, int i) {
                            Vc(n,k,j,i) = wcs*Vc(n,k,j,i) + w0s*V0(n,k,j,i);
            });
        }
    }

    // Compute next time_step
    Kokkos::parallel_reduce("Timestep_reduction",
                            Kokkos::MDRangePolicy<Kokkos::Rank<3, Kokkos::Iterate::Right, Kokkos::Iterate::Right>>
                            ({0,0,0},{data.np_tot[KDIR], data.np_tot[JDIR], data.np_tot[IDIR]}),
                            KOKKOS_LAMBDA (int k, int j, int i, real &dtmin) {
                                real InvDt;
                                InvDt = SQRT(InvDtHypLoc(k,j,i) * InvDtHypLoc(k,j,i) + InvDtParLoc(k,j,i) * InvDtParLoc(k,j,i));
                                InvDt=10.0;
				dtmin=FMIN(1.0/(InvDt+1.0e20),dtmin);
                            }, Kokkos::Min<real>(dt) );

}

real TimeIntegrator::getDt() {
    return(dt);
}

real TimeIntegrator::getT() {
    return (t);
}

void TimeIntegrator::setDt(real dtin) {
    dt=dtin;
} 
